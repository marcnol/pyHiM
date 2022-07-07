#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 15:26:00 2020

@author: marcnol
"""


import os
import argparse
from datetime import datetime
import logging

from dask.distributed import Client, LocalCluster, get_client, as_completed, fire_and_forget

from fileProcessing.fileManagement import (
    daskCluster,
    writeString2File,
    log,
    session,
    printDict,
    retrieveNumberUniqueBarcodesRootFolder,
    printLog,
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
from imageProcessing.segmentMasks3D import segmentMasks3D
from matrixOperations.filter_localizations import FilterLocalizations
from matrixOperations.register_localizations import RegisterLocalizations
from matrixOperations.build_traces import BuildTraces
from matrixOperations.build_matrix import BuildMatrix

class HiMfunctionCaller:
    def __init__(self, runParameters, sessionName="HiM_analysis"):
        self.runParameters = runParameters
        self.rootFolder = runParameters["rootFolder"]
        self.parallel = runParameters["parallel"]
        self.sessionName = sessionName

        self.log1 = log(rootFolder=self.rootFolder, parallel=self.parallel)

        self.session1 = session(self.rootFolder, self.sessionName)

    def initialize(self):

        printLog("\n--------------------------------------------------------------------------")

        printLog("$ rootFolder: {}".format(self.rootFolder))

        begin_time = datetime.now()

        #####################
        # setup markdown file
        #####################
        printLog("\n======================{}======================\n".format(self.sessionName))
        now = datetime.now()
        dateTime = now.strftime("%d%m%Y_%H%M%S")

        fileNameRoot="HiM_analysis"
        self.logFile = self.rootFolder + os.sep + fileNameRoot + dateTime + ".log"
        self.fileNameMD = self.logFile.split(".")[0] + ".md"

        printLog("$ Hi-M analysis will be written tos: {}".format(self.fileNameMD))
        writeString2File(
            self.fileNameMD, "# Hi-M analysis {}".format(begin_time.strftime("%Y/%m/%d %H:%M:%S")), "w",
        )  # initialises MD file

        ##############
        # setupLogger
        ##############

        # creates output formats for terminal and log file
        formatter1 = logging.Formatter("%(asctime)s: %(levelname)s: %(message)s")
        formatter2 = logging.Formatter("%(message)s")

        # clears up any existing logger
        logger = logging.getLogger()
        logger.handlers = []
        for hdlr in logger.handlers[:]:
            if isinstance(hdlr,logging.FileHandler):
                logger.removeHandler(hdlr)

        # initializes handlers for terminal and file
        filehandler = logging.FileHandler(self.logFile, 'w')
        ch = logging.StreamHandler()

        filehandler.setLevel(logging.INFO)
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)

        logger.addHandler(ch)
        logger.addHandler(filehandler)

        filehandler.setFormatter(formatter1)
        ch.setFormatter(formatter2)

    def lauchDaskScheduler(self,threadsRequested=25,maximumLoad=0.8):
        if self.parallel:
            printLog("$ Requested {} threads".format(threadsRequested))

            daskClusterInstance = daskCluster(threadsRequested,maximumLoad=maximumLoad)

            daskClusterInstance.createDistributedClient()
            self.client = daskClusterInstance.client
            self.cluster = daskClusterInstance.cluster

    def makeProjections(self, param):
        if not self.runParameters["parallel"]:
            makeProjections(param, self.session1)
        else:
            result = self.client.submit(makeProjections, param, self.session1)
            _ = self.client.gather(result)

    def alignImages(self, param, label):
        if label == "fiducial" and param.param["acquisition"]["label"] == "fiducial":
            printLog(
                "> Making image registrations for label: {}".format(label))
            if not self.parallel:
                alignImages(param, self.session1)
            else:
                result = self.client.submit(alignImages, param, self.session1)
                _ = self.client.gather(result)

    def alignImages3D(self, param, label):
        if label == "fiducial" and "block3D" in param.param["alignImages"]["localAlignment"]:
            printLog(
                "> Making 3D image registrations label: {}".format(label))
            _drift3D = drift3D(param, self.session1, parallel=self.parallel)
            _drift3D.alignFiducials3D()

    def appliesRegistrations(self, param, label):
        if label != "fiducial" and param.param["acquisition"]["label"] != "fiducial":
            printLog(
                "> Applying image registrations for label: {}".format(label))

            if not self.parallel:
                appliesRegistrations(param, self.session1)
            else:
                result = self.client.submit(appliesRegistrations, param, self.session1)
                _ = self.client.gather(result)

    def segmentMasks(self, param, label):
        if "segmentedObjects" in param.param.keys():
            operation = param.param["segmentedObjects"]["operation"]
        else:
            operation = [""]

        if (
            label != "RNA"
            and param.param["acquisition"]["label"] != "RNA"
            and "2D" in operation
        ):
            if not self.parallel:
                segmentMasks(param, self.session1)
            else:
                result = self.client.submit(segmentMasks, param, self.session1)
                _ = self.client.gather(result)



    def segmentMasks3D(self, param, label):
        if (
            (label == "DAPI" or label == 'mask')
            and "3D" in param.param["segmentedObjects"]["operation"]
        ):
            printLog("Making 3D image segmentations for label: {}".format(label))
            printLog(">>>>>>Label in functionCaller:{}".format(label))

            _segmentSources3D = segmentMasks3D(param, self.session1, parallel=self.parallel)
            _segmentSources3D.segmentMasks3D()


    def segmentSources3D(self, param, label):
        if (
            label == "barcode"
            and "3D" in param.param["segmentedObjects"]["operation"]
        ):
            printLog("Making 3D image segmentations for label: {}".format(label))
            printLog(">>>>>>Label in functionCaller:{}".format(label))

            _segmentSources3D = segmentSources3D(param, self.session1, parallel=self.parallel)
            _segmentSources3D.segmentSources3D()


    # This function will be removed in new release
    def projectsBarcodes(self, param, label):
        if label == "barcode":
            if not self.parallel:
                projectsBarcodes(param, self.log1, self.session1)
            else:
                result = self.client.submit(projectsBarcodes, param, self.log1, self.session1)
                _ = self.client.gather(result)

    # This function will be removed in new release
    def refitBarcodes(self, param, label):
        if label == "barcode":
            fittingSession = refitBarcodesClass(param, self.log1, self.session1, parallel=self.parallel)
            if not self.parallel:
                fittingSession.refitFolders()
            else:
                result = self.client.submit(fittingSession.refitFolders)
                _ = self.client.gather(result)

    # filters barcode localization table
    def filter_localizations(self, param, label):
        if label == "barcode":
            filter_localizations_instance = FilterLocalizations(param)
            filter_localizations_instance.filter_folder()

    # filters barcode localization table
    def register_localizations(self, param, label):
        if label == "barcode":
            register_localizations_instance = RegisterLocalizations(param)
            register_localizations_instance.register()

    # build traces
    def build_traces(self, param, label):
        if label == "barcode":
            build_traces_instance = BuildTraces(param)
            build_traces_instance.run()

    # build matrices
    def build_matrix(self, param, label):
        if label == "barcode":
            build_matrix_instance = BuildMatrix(param)
            build_matrix_instance.run()

    # This function will be removed in new release
    def localDriftCorrection(self, param, label):

        # runs mask 2D aligment
        if label == "DAPI" and ("mask2D" in param.param["alignImages"]["localAlignment"]):

            if not self.parallel:
                errorCode, _, _ = localDriftCorrection(param, self.log1, self.session1)
            else:
                result = self.client.submit(localDriftCorrection, param, self.log1, self.session1)
                errorCode, _, _ = self.client.gather(result)

    def processesPWDmatrices(self, param, label):
        if (label == "DAPI" or label == 'mask'):
            if not self.parallel:
                processesPWDmatrices(param, self.session1)
            else:
                result = self.client.submit(processesPWDmatrices, param, self.session1)
                _ = self.client.gather(result)

    def getLabel(self, ilabel):
        return self.labels2Process[ilabel]["label"]


def availableListCommands():
    return ["makeProjections", "appliesRegistrations","alignImages","alignImages3D", "segmentMasks",\
                "segmentMasks3D","segmentSources3D","refitBarcodes3D","localDriftCorrection",\
                "projectBarcodes","filter_localizations","register_localizations","build_traces","build_matrix","buildHiMmatrix"]


def defaultListCommands():
    return ["makeProjections", "appliesRegistrations","alignImages","alignImages3D", "segmentMasks",\
                "segmentMasks3D", "segmentSources3D","buildHiMmatrix"]

def HiM_parseArguments():
    parser = argparse.ArgumentParser()

    availableCommands=availableListCommands()
    defaultCommands=defaultListCommands()

    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("-C", "--cmd", help="Comma-separated list of routines to run (order matters !): makeProjections alignImages \
                        appliesRegistrations alignImages3D segmentMasks \
                        segmentMasks3D segmentSources3D buildHiMmatrix \
                        optional: [ filter_localizations register_localizations build_traces build_matrix]")
                        # to be removed: refitBarcodes3D localDriftCorrection projectBarcodes

    parser.add_argument("--threads", help="Number of threads to run in parallel mode. If none, then it will run with one thread.")
    args = parser.parse_args()

    printLog("\n--------------------------------------------------------------------------")
    runParameters = {}
    if args.rootFolder:
        runParameters["rootFolder"] = args.rootFolder
    else:
        if "docker" in os.environ.keys():
            # runParameters["rootFolder"] = os.environ["HiMdata"] #os.getenv("PWD")
            runParameters["rootFolder"] = "/data"
            printLog("\n\n$ Running in docker, HiMdata: {}".format(runParameters["rootFolder"]))
        else:
            printLog("\n\n# HiMdata: NOT FOUND")
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
        runParameters["cmd"] = defaultCommands

    for cmd in runParameters["cmd"]:
        if cmd not in availableCommands:
            printLog("\n\n# ERROR: {} not found in list of available commands: {}\n".format(cmd,availableCommands),status='WARN')
            raise SystemExit

    printDict(runParameters)

    return runParameters
