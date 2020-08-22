#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:48:55 2020

@author: marcnol
"""

""
# =============================================================================
# IMPORTS
# =============================================================================

import glob, os
import argparse
from datetime import datetime

from fileProcessing.fileManagement import folders, session, log, Parameters

from imageProcessing.localDriftCorrection import localDriftCorrection

# =============================================================================
# Local functions
# =============================================================================

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("-x", "--fileName", nargs='+', help="fileName to analyze")
    parser.add_argument("--parallel", help="Runs in parallel mode", action="store_true")

    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")
    runParameters={}
    if args.rootFolder:
        runParameters["rootFolder"] = args.rootFolder
    else:
        runParameters["rootFolder"] = os.getcwd()

    if args.fileName:
        runParameters["fileName"] = args.fileName
    else:
        runParameters["fileName"] = None
        
    if args.parallel:
        runParameters["parallel"] = args.parallel
    else:
        runParameters["parallel"] = False

    return runParameters


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    runParameters=parseArguments()    

    # rootFolder = "/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0"

    print("parameters> rootFolder: {}".format(runParameters["rootFolder"]))
    sessionName = "localDriftCorrection"

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
    ]

    # session
    session1 = session(runParameters["rootFolder"], sessionName)

    # setup logs
    log1 = log(runParameters["rootFolder"])
    # labels2Process indeces: 0 fiducial, 1:
    labelParameterFile = labels2Process[2]["parameterFile"]
    param = Parameters(runParameters["rootFolder"], labelParameterFile)

    dataFolder = folders(param.param["rootFolder"])

    for currentFolder in dataFolder.listFolders:
        filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
        dataFolder.createsFolders(currentFolder, param)

        # generates lists of files to process
        param.files2Process(filesFolder)
        if runParameters["parallel"]:
            param.param['parallel']=True
        else:
            param.param['parallel']=False
            
        for fileName in param.fileList2Process:
            session1.add(fileName, sessionName)

        errorCode, dictShift, alignmentResultsTable = localDriftCorrection(param, log1, session1)

    if errorCode != 0:
        print("Error code reported: {}".format(errorCode))
    else:
        print("normal termination")

    # for fileName in param.fileList2Process:
    #     session1.add(fileName, sessionName)
    print("Elapsed time: {}".format(datetime.now() - begin_time))
