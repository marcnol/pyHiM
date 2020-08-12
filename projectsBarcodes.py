#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 20:48:39 2020

This function will project all barcodes into a single image
This will be used (maybe) to then segment barcodes belonging to the same chromosome
or cluster. This could be used alternatively or in complement to DAPI segmentation

@author: marcnol
"""

# =============================================================================
# IMPORTS
# =============================================================================
import glob, os
import argparse

# import matplotlib.pylab as plt
# import numpy as np
# from astropy.stats import sigma_clipped_stats,SigmaClip,gaussian_fwhm_to_sigma
# from astropy.convolution import Gaussian2DKernel
# from astropy.visualization import SqrtStretch,simple_norm
# from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.table import Table, vstack, Column

# from photutils import DAOStarFinder,CircularAperture,detect_sources
# from photutils import detect_threshold,deblend_sources
# from photutils import Background2D, MedianBackground
from imageProcessing import Image, saveImage2Dcmd, imageAdjust
from fileManagement import folders, session, log, Parameters
from fileManagement import writeString2File, FileHandling


# =============================================================================
# FUNCTIONS
# =============================================================================


def gets2DcorrectedImage(fileName, param, log1, session1, dataFolder):
    """
    Reads new 2D projected, drift-corrected image from a barcode and returns it

    Parameters
    ----------
    fileName : TYPE
        DESCRIPTION.
    param : TYPE
        DESCRIPTION.
    log1 : TYPE
        DESCRIPTION.
    session1 : TYPE
        DESCRIPTION.
    dataFolder : TYPE
        DESCRIPTION.

    Returns
    -------
    im=image;
    errorCode=0, normal termination; -1 image not found.

    """

    rootFileName = os.path.basename(fileName).split(".")[0]
    fileName_2d_aligned = dataFolder.outputFolders["alignImages"] + os.sep + rootFileName + "_2d_registered.npy"

    if os.path.exists(fileName_2d_aligned):  # file exists
        # loading registered 2D projection
        Im = Image()
        Im.loadImage2D(
            fileName, log1, dataFolder.outputFolders["alignImages"], tag="_2d_registered",
        )
        im = Im.data_2D
        del Im
        return im, 0

    else:
        return [], -1


def projectsBarcodes(param, log1, session1):

    if param.param["projectsBarcodes"]["operation"] == "overwrite":
        sessionName = "projectsBarcodes"

        # processes folders and files
        dataFolder = folders(param.param["rootFolder"])
        log1.addSimpleText(
            "\n===================={}:{}====================\n".format(sessionName, param.param["acquisition"]["label"])
        )
        log1.report("folders read: {}".format(len(dataFolder.listFolders)))
        writeString2File(
            log1.fileNameMD, "## {}: {}\n".format(sessionName, param.param["acquisition"]["label"]), "a",
        )
        # barcodesCoordinates=Table()

        for currentFolder in dataFolder.listFolders:
            # currentFolder=dataFolder.listFolders[0]
            filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
            dataFolder.createsFolders(currentFolder, param)
            log1.report("-------> Processing Folder: {}".format(currentFolder))

            # generates lists of files to process
            param.files2Process(filesFolder)
            log1.report("About to read {} files\n".format(len(param.fileList2Process)))

            # should find out how many ROIs are there
            imageStack = {}

            for fileName in param.fileList2Process:
                # gets ROI
                newFile = FileHandling(fileName)
                ROI = newFile.getROI()

                newImage, errorCode = gets2DcorrectedImage(fileName, param, log1, session1, dataFolder)

                if errorCode == 0:
                    # adjusts image levels before stacking
                    (newImage, hist1_before, hist1, lower_cutoff, higher_cutoff,) = imageAdjust(
                        newImage, lower_threshold=0.3, higher_threshold=0.9999
                    )

                    # convolve with a gaussian kernel to homogeneize region occupied by barcodes

                    # accumulates new image to stack by summing
                    if ROI in imageStack:
                        imageStack[ROI] += newImage  # accumulates images
                    else:
                        imageStack[ROI] = newImage  # starts dict entry
                    log1.report(
                        "File {} accumulated to stack of ROI {}.".format(fileName, ROI), "info",
                    )

                else:
                    log1.report(
                        "No 2d corrected image for file {} could be found --> not accumulated".format(fileName),
                        "warning",
                    )

                # saves projected image into file
                session1.add(fileName, sessionName)

            # [saves imageStack to file before starting a new currentFolder ]
            for ROI in imageStack:
                imageFileNameOutput = dataFolder.outputFiles["projectsBarcodes"] + "_" + ROI + ".npy"
                saveImage2Dcmd(imageStack[ROI], imageFileNameOutput, log1)

                ImtoSave = Image()
                ImtoSave.data_2D = imageStack[ROI]
                outputName = dataFolder.outputFiles["projectsBarcodes"] + "_" + ROI + ".png"
                ImtoSave.imageShow(outputName=outputName, normalization="simple")

                log1.report("Output image File {}".format(outputName), "info")
                del ImtoSave


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
        rootFolder = "/home/marcnol/data/Experiment_20/Embryo_1"
        # rootFolder='/home/marcnol/data/Experiment_15/Embryo_006_ROI18'
        # rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset'

    print("parameters> rootFolder: {}".format(rootFolder))

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
    ]

    # session
    session1 = session(rootFolder, "processingPipeline")

    # setup logs
    log1 = log(rootFolder)

    labelParameterFile = labels2Process[1]["parameterFile"]
    param = Parameters(rootFolder, labelParameterFile)

    projectsBarcodes(param, log1, session1)
