#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 18:16:10 2020

@author: marcnol

"""


# Imports
"""import logging
logging.basicConfig(level=logging.INFO
                #,format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
                #,datefmt='%d/%m/%Y %I:%M:%S %p'
                #,filename='testImageProcessing.log'
                )  
"""

import glob, os, sys

# from copy import deepcopy
import matplotlib.pylab as plt
import numpy as np
import cv2
import matplotlib.pyplot as plt
from matplotlib import cm
from skimage import io
from imageProcessing import Image
import parameters as parameters
from fileManagement import log
from fileManagement import folders

# import cPickle as pickle
# from scipy.spatial import ConvexHull
# from scipy.spatial.distance import cdist
# from tqdm import tqdm_notebook as tqdm

# Internal packages
# import Fitting_v4 as ft
# import workers_cells_v3 as wkc


def processImage(fileName, param, log1):
    log1.report("Analysing file: {}\n".format(fileName))

    # creates image object
    Im = Image()

    # loads image
    Im.loadImage(fileName)

    Im.zProjectionRange(param, log1)

    if Im.fileName:
        Im.printImageProperties()

    Im.imageShow(save=False)

    del Im


if __name__ == "__main__":

    # defines folders
    rootFolder = "/home/marcnol/Documents/Images"
    logFile = "testImageProcess.log"

    # sets parameters
    param = parameters.Parameters()
    param.initializeStandardParameters()
    paramFile = rootFolder + os.sep + "infoList.inf"
    param.loadParametersFile(paramFile)

    # setup logs
    logFileName = rootFolder + os.sep + logFile
    log1 = log(logFileName)
    log1.eraseFile()

    # processes folders and files
    dataFolder = folders(rootFolder)
    dataFolder.setsFolders()
    log1.report("folders read: {}".format(len(dataFolder.listFolders)))
    filesFolder = glob.glob(dataFolder.listFolders[0] + os.sep + "*.tif")
    log1.report("About to read {} files\n".format(len(filesFolder)))

    # Processes all DAPI masks
    for fileName in filesFolder:
        if fileName.split("_")[-1].split(".")[0] == "ch00" and "DAPI" in fileName.split("_"):
            processImage(fileName, param, log1)

    # Processes all DAPI masks
    for fileName in filesFolder:
        res = [i for i in fileName.split("_") if "RT" in i]
        if len(res) > 0 and fileName.split("_")[-1].split(".")[0] == "ch00":
            barcodeNumber = res[0].split("RT")[1]
            log1.report("Processing barcode: {}".format(barcodeNumber))
            processImage(fileName, param, log1)

    # exits
    log1.report("Normal exit.")
