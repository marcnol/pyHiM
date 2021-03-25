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

from fileProcessing.fileManagement import (
    daskCluster,
    writeString2File,
    log,
    session,
    retrieveNumberUniqueBarcodesRootFolder,
)

from imageProcessing.alignImages import alignImages, appliesRegistrations
from imageProcessing.makeProjections import makeProjections
from imageProcessing.segmentMasks import segmentMasks
from imageProcessing.localDriftCorrection import localDriftCorrection
from imageProcessing.projectsBarcodes import projectsBarcodes
from matrixOperations.alignBarcodesMasks import processesPWDmatrices
from imageProcessing.refitBarcodes3D import refitBarcodesClass
from imageProcessing.alignImages3D import drift3D
from imageProcessing.segmentSources3D import segmentSources3D

class HiMfunctionCaller:
    def __init__(self, runParameters, sessionName="HiM_analysis"):
        self.runParameters = runParameters
        self.rootFolder = runParameters["rootFolder"]
        self.parallel = runParameters["parallel"]
        self.sessionName = sessionName

        self.log1 = log(rootFolder=self.rootFolder, parallel=self.parallel)

        self.labels2Process = [
            {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
            {"label": "barcode", "parameterFile": "infoList_barcode.json"},
            {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
            {"label": "RNA", "parameterFile": "infoList_RNA.json"},
        ]

        self.session1 = session(self.rootFolder, self.sessionName)

    def initialize(self):

        print("\n--------------------------------------------------------------------------")

        print("$ rootFolder: {}".format(self.rootFolder))

        begin_time = datetime.now()

        # setup logs
        # log1 = log(rootFolder = self.rootFolder,parallel=self.parallel)
        self.log1.addSimpleText(
            "\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format(self.sessionName)
        )
        if self.log1.fileNameMD == ".md":
            self.log1.fileNameMD = "HiM_report.md"

        self.log1.addSimpleText("$ Hi-M analysis will be written tos: {}".format(self.log1.fileNameMD))
        writeString2File(
            self.log1.fileNameMD, "# Hi-M analysis {}".format(begin_time.strftime("%Y/%m/%d %H:%M:%S")), "w",
        )  # initialises MD file

    def lauchDaskScheduler(self,threadsRequested=25,maximumLoad=0.8):
        if self.parallel:
            print("$ Requested {} threads".format(threadsRequested))

            daskClusterInstance = daskCluster(threadsRequested,maximumLoad=maximumLoad)

            daskClusterInstance.createDistributedClient()
            self.client = daskClusterInstance.client
            self.cluster = daskClusterInstance.cluster

    def makeProjections(self, param):
        if not self.runParameters["parallel"]:
            makeProjections(param, self.log1, self.session1)
        else:
            result = self.client.submit(makeProjections, param, self.log1, self.session1)
            _ = self.client.gather(result)

    def alignImages(self, param, ilabel):
        if self.getLabel(ilabel) == "fiducial" and param.param["acquisition"]["label"] == "fiducial":
            self.log1.addSimpleText(
                "> Making image registrations, ilabel: {}, label: {}".format(ilabel, self.getLabel(ilabel)))
            if not self.parallel:
                alignImages(param, self.log1, self.session1)
            else:
                result = self.client.submit(alignImages, param, self.log1, self.session1)
                _ = self.client.gather(result)

    def alignImages3D(self, param, ilabel):
        if self.getLabel(ilabel) == "fiducial" and "block3D" in param.param["alignImages"]["localAlignment"]:
            self.log1.addSimpleText(
                "> Making 3D image registrations, ilabel: {}, label: {}".format(ilabel, self.getLabel(ilabel)))
            _drift3D = drift3D(param, self.log1, self.session1, parallel=self.parallel)
            # if not self.parallel:
            _drift3D.alignFiducials3D()
            # else:
            #     result = self.client.submit(_drift3D.alignFiducials3D)
            #     _ = self.client.gather(result)

    def appliesRegistrations(self, param, ilabel):
        if self.getLabel(ilabel) != "fiducial" and param.param["acquisition"]["label"] != "fiducial":
            self.log1.addSimpleText(
                "> Applying image registrations, ilabel: {}, label: {}".format(ilabel, self.getLabel(ilabel)))

            if not self.parallel:
                appliesRegistrations(param, self.log1, self.session1)
            else:
                result = self.client.submit(appliesRegistrations, param, self.log1, self.session1)
                _ = self.client.gather(result)

    def segmentMasks(self, param, ilabel):

        if "segmentedObjects" in param.param.keys():
            operation = param.param["segmentedObjects"]["operation"]
        else:
            operation = [""]

        if (
            # self.getLabel(ilabel) != "fiducial"
            # and param.param["acquisition"]["label"] != "fiducial"
            self.getLabel(ilabel) != "RNA"
            and param.param["acquisition"]["label"] != "RNA"
            and "2D" in operation
        ):
            if not self.parallel:
                segmentMasks(param, self.log1, self.session1)
            else:
                result = self.client.submit(segmentMasks, param, self.log1, self.session1)
                _ = self.client.gather(result)


    def segmentSources3D(self, param, ilabel):
        if (
            self.getLabel(ilabel) == "barcode"
            and "3D" in param.param["segmentedObjects"]["operation"]
        ):
            self.log1.report(
                "Making 3D image segmentations, ilabel: {}, label: {}".format(ilabel, self.getLabel(ilabel)), "info"
            )
            _segmentSources3D = segmentSources3D(param, self.log1, self.session1, parallel=self.parallel)
            # if not self.parallel:
            _segmentSources3D.segmentSources3D()
            # else:
            #     result = self.client.submit(_segmentSources3D.segmentSources3D)
            #     _ = self.client.gather(result)


    def projectsBarcodes(self, param, ilabel):
        if self.getLabel(ilabel) == "barcode":
            if not self.parallel:
                projectsBarcodes(param, self.log1, self.session1)
            else:
                result = self.client.submit(projectsBarcodes, param, self.log1, self.session1)
                _ = self.client.gather(result)

    def refitBarcodes(self, param, ilabel):
        if self.getLabel(ilabel) == "barcode" and self.runParameters["refit"]:
            fittingSession = refitBarcodesClass(param, self.log1, self.session1, parallel=self.parallel)
            if not self.parallel:
                fittingSession.refitFolders()
            else:
                result = self.client.submit(fittingSession.refitFolders)
                _ = self.client.gather(result)

    def localDriftCorrection(self, param, ilabel):

        # runs mask 2D aligment
        if self.getLabel(ilabel) == "DAPI" and ("mask2D" in param.param["alignImages"]["localAlignment"]):

            if not self.parallel:
                errorCode, _, _ = localDriftCorrection(param, self.log1, self.session1)
            else:
                result = self.client.submit(localDriftCorrection, param, self.log1, self.session1)
                errorCode, _, _ = self.client.gather(result)

    def processesPWDmatrices(self, param, ilabel):
        if self.getLabel(ilabel) == "DAPI":
            if not self.parallel:
                processesPWDmatrices(param, self.log1, self.session1)
            else:
                result = self.client.submit(processesPWDmatrices, param, self.log1, self.session1)
                a = self.client.gather(result)

    def getLabel(self, ilabel):
        return self.labels2Process[ilabel]["label"]


def makeListCommands():
    return ["makeProjections", "appliesRegistrations","alignImages","alignImages3D", "segmentMasks",\
                "segmentSources3D","refitBarcodes3D","localDriftCorrection","projectBarcodes","buildHiMmatrix"]

def HiM_parseArguments():
    parser = argparse.ArgumentParser()

    available_commands = makeListCommands()

    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("-C", "--cmd", help="Comma-separated list of routines to run (order matters !): makeProjections, appliesRegistrations,\
                        alignImages, alignImages3D,  segmentMasks,\
                        segmentSources3D, refitBarcodes3D,\
                        localDriftCorrection, projectBarcodes, buildHiMmatrix")
    parser.add_argument("--threads", help="Number of threads to run in parallel mode. If none, then it will run with one thread.")
    # parser.add_argument("--localAlignment", help="Runs localAlignment function", action="store_true")
    # parser.add_argument("--refit", help="Refits barcode spots using a Gaussian axial fitting function.", action="store_true")

    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")
    runParameters = {}
    if args.rootFolder:
        runParameters["rootFolder"] = args.rootFolder
    else:
        if "docker" in os.environ.keys():
            # runParameters["rootFolder"] = os.environ["HiMdata"] #os.getenv("PWD")
            runParameters["rootFolder"] = "/data"
            print("\n\n$ Running in docker, HiMdata: {}".format(runParameters["rootFolder"]))
        else:
            print("\n\n# HiMdata: NOT FOUND")
            runParameters["rootFolder"] = os.getenv("PWD")  # os.getcwd()

    if args.threads:
        runParameters["threads"] = int(args.threads)
        runParameters["parallel"] = True
    else:
        runParameters["threads"] = 1
        runParameters["parallel"] = False

    if args.cmd:
        runParameters["cmd"] = args.cmd.split(",")
    else:
        runParameters["cmd"] = available_commands

    for cmd in runParameters["cmd"]:
        if cmd not in available_commands:
            print("\n\n# ERROR: {} not found in list of available commands: {}\n".format(cmd,available_commands))
            raise SystemExit

    # if args.localAlignment:
    #     runParameters["localAlignment"] = args.localAlignment
    # else:
    #     runParameters["localAlignment"] = False

    # if args.refit:
    #     runParameters["refit"] = args.refit
    # else:
    #     runParameters["refit"] = False

    print("\n$ Parameters loaded:")
    for key in runParameters.keys():
        if len(key)>7:
            print("\t{}\t{}".format(key,runParameters[key]))
        else:
            print("\t{}\t\t{}".format(key,runParameters[key]))
    print("\n")

    return runParameters
