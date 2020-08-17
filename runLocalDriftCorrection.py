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

from fileProcessing.fileManagement import folders, session, log, Parameters

from imageProcessing.localDriftCorrection import localDriftCorrection

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        # rootFolder = "/home/marcnol/data/Experiment_20/Embryo_1"
        # rootFolder='/home/marcnol/data/Experiment_15/Embryo_006_ROI18'
        rootFolder = "/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0"
        # rootFolder='/home/marcnol/data/Embryo_debug_dataset/rawImages'

    print("parameters> rootFolder: {}".format(rootFolder))
    sessionName = "localDriftCorrection"

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
    ]

    # session
    session1 = session(rootFolder, sessionName)

    # setup logs
    log1 = log(rootFolder)
    # labels2Process indeces: 0 fiducial, 1:
    labelParameterFile = labels2Process[2]["parameterFile"]
    param = Parameters(rootFolder, labelParameterFile)

    dataFolder = folders(param.param["rootFolder"])

    for currentFolder in dataFolder.listFolders:
        filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
        dataFolder.createsFolders(currentFolder, param)

        # generates lists of files to process
        param.files2Process(filesFolder)

        for fileName in param.fileList2Process:
            session1.add(fileName, sessionName)

    errorCode, dictShift, alignmentResultsTable = localDriftCorrection(param, log1, session1)

    if errorCode != 0:
        print("Error code reported: {}".format(errorCode))
    else:
        print("normal termination")

    # for fileName in param.fileList2Process:
    #     session1.add(fileName, sessionName)
