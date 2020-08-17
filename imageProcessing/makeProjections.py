#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 23:17:58 2020

@author: marcnol

This file contains functions to project 3D images to 2D

Operation will be defined in the parameters file. Options are:
    - user-defined range
    - all z range
    - optimal range based on detection of focal plane and use of user defined window around it
    

"""
# =============================================================================
# IMPORTS
# =============================================================================


import glob, os

# import matplotlib.pylab as plt
import numpy as np
import argparse
from datetime import datetime

# import cv2
import matplotlib.pyplot as plt
from imageProcessing.imageProcessing import Image

from fileProcessing.fileManagement import (
    folders,session, writeString2File, folders, session, log, Parameters)


# =============================================================================
# FUNCTIONS
# =============================================================================


def makes2DProjectionsFile(fileName, param, log1, session1, dataFolder):

    if fileName in session1.data and param.param["zProject"]["operation"] != "overwrite":
        # creates image object
        Im = Image()
        Im.loadImage2D(fileName, log1, dataFolder.outputFolders["zProject"])
        if param.param["zProject"]["display"]:
            Im.imageShow()
        log1.report("File already projected: {}".format(os.path.basename(fileName)))
    else:

        log1.report("Analysing file: {}".format(os.path.basename(fileName)))

        # creates image object
        Im = Image()

        # loads image
        Im.loadImage(fileName)

        # makes actual 2d projection
        Im.zProjectionRange(param, log1)

        # outputs information from file
        if Im.fileName:
            Im.printImageProperties()

        # saves output 2d zProjection as png
        if param.param["zProject"]["display"]:
            pngFileName = dataFolder.outputFolders["zProject"] + os.sep + os.path.basename(fileName) + "_2d.png"
            Im.imageShow(save=param.param["zProject"]["saveImage"], outputName=pngFileName)
            writeString2File(
                log1.fileNameMD, "{}\n ![]({})\n".format(os.path.basename(fileName), pngFileName), "a",
            )  # initialises MD file

        # saves output 2d zProjection as matrix
        Im.saveImage2D(log1, dataFolder.outputFolders["zProject"])

        del Im


def makeProjections(param, log1, session1,fileName=None):
    sessionName = "makesProjections"

    # processes folders and files
    dataFolder = folders(param.param["rootFolder"])
    log1.addSimpleText("\n===================={}====================\n".format(sessionName))
    log1.report("folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(
        log1.fileNameMD, "## {}: {}\n".format(sessionName, param.param["acquisition"]["label"]), "a",
    )  # initialises MD file

    for currentFolder in dataFolder.listFolders:
        filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
        dataFolder.createsFolders(currentFolder, param)

        # generates lists of files to process
        param.files2Process(filesFolder)
        log1.report("-------> Processing Folder: {}".format(currentFolder))
        log1.report("About to read {} files\n".format(len(param.fileList2Process)))

        for fileName2Process in param.fileList2Process:
            # print("Looking for {} in {}".format(fileName,fileName2Process))

            if fileName==None:
                makes2DProjectionsFile(fileName2Process, param, log1, session1, dataFolder)
                session1.add(fileName2Process, sessionName)
            elif fileName!=None and (os.path.basename(fileName2Process) in [os.path.basename(x) for x in fileName]):
                makes2DProjectionsFile(fileName2Process, param, log1, session1, dataFolder)
                session1.add(fileName2Process, sessionName)
                # print("******File {} processed!!!".format(fileName2Process))
            else:
                pass
                # print("File {} not found in list".format(fileName2Process))


# =============================================================================
# MAIN
# =============================================================================

# if __name__ == "__main__":

#     begin_time = datetime.now()

#     parser = argparse.ArgumentParser()
#     parser.add_argument("-F", "--rootFolder", help="Folder with images")
#     parser.add_argument("-x", "--fileName", help="fileName to analyze")
#     args = parser.parse_args()

#     print("\n--------------------------------------------------------------------------")

#     if args.rootFolder:
#         rootFolder = args.rootFolder
#     else:
#         rootFolder = os.getcwd()

#     if args.fileName:
#         fileName = args.fileName
#     else:
#         fileName = None
        
#         print("parameters> rootFolder: {}".format(rootFolder))
#     now = datetime.now()

#     labels2Process = [
#         {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
#         {"label": "barcode", "parameterFile": "infoList_barcode.json"},
#         {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
#         {"label": "RNA", "parameterFile": "infoList_RNA.json"},
#     ]

#     # session
#     sessionName = "makesProjections"
#     session1 = session(rootFolder, sessionName)

#     # setup logs
#     log1 = log(rootFolder)
#     log1.addSimpleText("\n-------------------------{}-------------------------\n".format(sessionName))
#     log1.report("Hi-M analysis MD: {}".format(log1.fileNameMD))
#     writeString2File(
#         log1.fileNameMD, "# Hi-M analysis {}".format(now.strftime("%Y/%m/%d %H:%M:%S")), "w",
#     )  # initialises MD file

#     for ilabel in range(len(labels2Process)):
#         label = labels2Process[ilabel]["label"]
#         labelParameterFile = labels2Process[ilabel]["parameterFile"]
#         log1.addSimpleText("**Analyzing label: {}**".format(label))

#         # sets parameters
#         param = Parameters(rootFolder, labelParameterFile)

#         # [projects 3D images in 2d]
#         makeProjections(param, log1, session1, fileName)

#         print("\n")
#         del param
#     # exits
#     session1.save(log1)
#     log1.addSimpleText("\n===================={}====================\n".format("Normal termination"))

#     del log1, session1
#     print("Elapsed time: {}".format(datetime.now() - begin_time))
