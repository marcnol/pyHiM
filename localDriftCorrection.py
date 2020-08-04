#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:45:32 2020

@author: marcnol


Purpose: Loads masks and performs a local drift correction in the region surrounding the mask.

steps for each ROI: 
    - load reference fiducial 2D projection
    - iterate over cycles
    - load fiducial for cycle <i>
    - check if masks are available
    - load mask file 
    - iterate over masks <j>
    - obtain subvolume from reference and cycle <i> fiducial images for mask <j>
    - cross-correlate and obtain a second-order correction value
    - determine if we keep or not
    - store in database.


"""


import glob, os
import matplotlib.pylab as plt
import numpy as np
import uuid
import argparse

from imageProcessing import Image, saveImage2Dcmd
from fileManagement import folders,session, log, Parameters
from fileManagement import writeString2File,ROI2ReferenceFiducialFileName


def loadsFiducial(param,fileName):
    
    # finds name of reference fiducial file
    positionROIinformation = param.param["acquisition"]["positionROIinformation"]
    ROI = os.path.basename(fileName).split("_")[positionROIinformation]
    referenceBarcode = param.param["alignImages"]["referenceFiducial"]
    fiducialFilename = ROI2ReferenceFiducialFileName(param, fileName, referenceBarcode, positionROIinformation)
    
    # loads reference fiducial file
    if len(fiducialFilename) < 1:
        print('Error, no reference candidate found for ROI:{} | filename:{}\n'.format(ROI, os.path.basename(fileName)))
        return -1
    elif len(fiducialFilename) > 1:
        print('Error, too many reference candidates found for ROI:{} | filename:{}\n'.format(ROI, os.path.basename(fileName)))
        return -1
    else:
        print('Using reference fiducial> {}\n'.format(os.path.basename(fiducialFilename[0])))
        fullFiducialFilename = currentFolder+os.sep+fiducialFilename[0]
    
    imReference = Image()
    imReference.loadImage2D(fullFiducialFilename , log1, dataFolder.outputFolders["zProject"])
    imReference.imageShow(show=True)
    
    return ROI, imReference 


def localDriftCorrection(param, log1, session1):
    
    # processes folders and files
    dataFolder = folders(param.param["rootFolder"])
    log1.addSimpleText(
        "\n===================={}:{}====================\n".format(sessionName, param.param["acquisition"]["label"])
    )
    log1.report("folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(
        log1.fileNameMD, "## {}: {}\n".format(session1.name, param.param["acquisition"]["label"]), "a",
    )
    
    for currentFolder in dataFolder.listFolders:
        filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
        dataFolder.createsFolders(currentFolder, param)

        # generates lists of files to process
        param.files2Process(filesFolder)
        log1.report("-------> Processing Folder: {}".format(currentFolder))
        log1.report("About to read {} files\n".format(len(param.fileList2Process)))

        # iterates over ROIs
        for fileName in param.fileList2Process:
            
            label = param.param["acquisition"]["label"]

            if label == "DAPI":
                # appends filename to session and defines ROI
                session1.add(fileName, sessionName)
                print('processing> {}\n'.format(os.path.basename(fileName)))

                # loads reference fiducial file
                ROI, imReference = loadsFiducial(param,fileName)
                
                # - iterate over cycles

                # - load fiducial for cycle <i>

                # - check if masks are available

                # - load mask file 

                # - iterate over masks <j>

                # - obtain subvolume from reference and cycle <i> fiducial images for mask <j>

                # - cross-correlate and obtain a second-order correction value

                # - determine if we keep or not

                # - store in database.                

    return 0
                
                
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
        # rootFolder='/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0'
        rootFolder='/home/marcnol/data/Embryo_debug_dataset/rawImages'
        
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


    errorCode=localDriftCorrection(param, log1, session1)
            
    if errorCode!=0:
        print("Error code reported: {}".format(errorCode))
    else:
        print("normal termination")
            
    # for fileName in param.fileList2Process:
    #     session1.add(fileName, sessionName)


