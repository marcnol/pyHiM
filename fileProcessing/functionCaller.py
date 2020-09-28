#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 15:26:00 2020

@author: marcnol
"""


import os
import argparse
from datetime import datetime

from dask.distributed import Client, LocalCluster, get_client, as_completed, fire_and_forget

from fileProcessing.fileManagement import daskCluster, writeString2File, log, session, retrieveNumberUniqueBarcodesRootFolder

from imageProcessing.alignImages import alignImages, appliesRegistrations
from imageProcessing.makeProjections import makeProjections
from imageProcessing.segmentMasks import segmentMasks
from imageProcessing.localDriftCorrection import localDriftCorrection
from imageProcessing.projectsBarcodes import projectsBarcodes
from matrixOperations.alignBarcodesMasks import processesPWDmatrices
from imageProcessing.refitBarcodes3D import refitBarcodesClass


class HiMfunctionCaller:
    def __init__(self, runParameters, sessionName="HiM_analysis"):
        self.runParameters=runParameters
        self.rootFolder=runParameters["rootFolder"]
        self.parallel=runParameters["parallel"]
        self.sessionName=sessionName
        
        self.log1 = log(rootFolder = self.rootFolder,parallel=self.parallel)
        
        self.labels2Process = [
            {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
            {"label": "barcode", "parameterFile": "infoList_barcode.json"},
            {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
            {"label": "RNA", "parameterFile": "infoList_RNA.json"},
        ]
        
        self.session1 = session(self.rootFolder, self.sessionName)
          
    def initialize(self):

        print("\n--------------------------------------------------------------------------")

        print("parameters> rootFolder: {}".format(self.rootFolder))

        begin_time = datetime.now()

        # setup logs
        # log1 = log(rootFolder = self.rootFolder,parallel=self.parallel)
        self.log1.addSimpleText("\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format(self.sessionName))
        if self.log1.fileNameMD=='.md':
            self.log1.fileNameMD='HiM_report.md'
            
        self.log1.report("Hi-M analysis MD: {}".format(self.log1.fileNameMD))
        writeString2File(
            self.log1.fileNameMD, "# Hi-M analysis {}".format(begin_time.strftime("%Y/%m/%d %H:%M:%S")), "w",
        )  # initialises MD file
        
        
    def lauchDaskScheduler(self):
        if self.parallel:
            parametersFile = self.rootFolder + os.sep + self.labels2Process[0]["parameterFile"]
            numberUniqueCycles = retrieveNumberUniqueBarcodesRootFolder(self.rootFolder,parametersFile)
            print("Found {} unique cycles in rootFolder".format(numberUniqueCycles))
            self.daskClusterInstance = daskCluster(numberUniqueCycles)
            print("Go to http://localhost:8787/status for information on progress...")
            
            self.cluster = LocalCluster(n_workers=self.daskClusterInstance.nThreads,
                                # processes=True,
                                # threads_per_worker=1,
                                # memory_limit='2GB',
                                # ip='tcp://localhost:8787',
                                ) 
            self.client = Client(self.cluster)
            
    def makeProjections(self,param):
        if not self.runParameters["parallel"]:
            makeProjections(param, self.log1, self.session1)
        else:
            result = self.client.submit(makeProjections,param, self.log1, self.session1)
            _ = self.client.gather(result)
        
    def alignImages(self, param, ilabel):
        if self.getLabel(ilabel) == "fiducial" and param.param["acquisition"]["label"] == "fiducial":
            self.log1.report("Making image registrations, ilabel: {}, label: {}".format(ilabel, self.getLabel(ilabel)), "info")
            if not self.parallel:
                alignImages(param, self.log1, self.session1)        
            else:
                 result = self.client.submit(alignImages,param, self.log1, self.session1)
                 _ = self.client.gather(result)

    def appliesRegistrations(self, param, ilabel):
        if self.getLabel(ilabel) != "fiducial" and param.param["acquisition"]["label"] != "fiducial":
            self.log1.report("Applying image registrations, ilabel: {}, label: {}".format(ilabel, self.getLabel(ilabel)), "info")

            if not self.parallel:
                appliesRegistrations(param, self.log1, self.session1)
            else:
                result = self.client.submit(appliesRegistrations,param, self.log1, self.session1)
                _ = self.client.gather(result)

    def segmentMasks(self, param, ilabel):
        if (self.getLabel(ilabel)!= "fiducial" and \
            param.param["acquisition"]["label"] != "fiducial" and \
            self.getLabel(ilabel)!= "RNA" and \
            param.param["acquisition"]["label"] != "RNA"):
            if not self.parallel:
                segmentMasks(param, self.log1, self.session1)
            else:
                result = self.client.submit(segmentMasks,param, self.log1, self.session1)
                _ = self.client.gather(result)

    def projectsBarcodes(self, param, ilabel):
        if self.getLabel(ilabel)== "barcode":
            if not self.parallel:
                projectsBarcodes(param, self.log1, self.session1)
            else:
                result = self.client.submit(projectsBarcodes,param, self.log1, self.session1)
                _ = self.client.gather(result)

                                
    def refitBarcodes(self, param, ilabel):
        if self.getLabel(ilabel) == "barcode" and self.runParameters["refit"]:
            fittingSession = refitBarcodesClass(param, self.log1, self.session1,parallel=self.parallel)
            if not self.parallel:
                fittingSession.refitFolders()            
            else:
                result = self.client.submit(fittingSession.refitFolders)
                _ = self.client.gather(result)
                
    def localDriftCorrection(self, param, ilabel):
        if self.getLabel(ilabel) == "DAPI" and self.runParameters["localAlignment"]:

            if not self.parallel:
                errorCode, _, _ = localDriftCorrection(param, self.log1, self.session1)                
            else:
                result = self.client.submit(localDriftCorrection,param, self.log1, self.session1)
                errorCode, _, _ = self.client.gather(result)

    def processesPWDmatrices(self, param, ilabel):
        if self.getLabel(ilabel) == "DAPI":
            if not self.parallel:
                processesPWDmatrices(param, self.log1, self.session1)
            else:
                result = self.client.submit(processesPWDmatrices,param, self.log1, self.session1)
                a = self.client.gather(result)
                
    def getLabel(self, ilabel):
        return self.labels2Process[ilabel]["label"]

def HiM_parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("--parallel", help="Runs in parallel mode", action="store_true")
    parser.add_argument("--localAlignment", help="Runs localAlignment function", action="store_true")
    parser.add_argument("--refit", help="Refits barcode spots using a Gaussian axial fitting function.", action="store_true")
    
    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")
    runParameters={}
    if args.rootFolder:
        runParameters["rootFolder"] = args.rootFolder
    else:
        runParameters["rootFolder"] = '.' # os.getcwd()
       
    if args.parallel:
        runParameters["parallel"] = args.parallel
    else:
        runParameters["parallel"] = False

    if args.localAlignment:
        runParameters["localAlignment"] = args.localAlignment
    else:
        runParameters["localAlignment"] = False

    if args.refit:
        runParameters["refit"] = args.refit
    else:
        runParameters["refit"] = False

    return runParameters
