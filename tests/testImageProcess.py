#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 18:16:10 2020

@author: marcnol

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
from fileManagement import Log
from fileManagement import Folders

# import cPickle as pickle
# from scipy.spatial import ConvexHull
# from scipy.spatial.distance import cdist
# from tqdm import tqdm_notebook as tqdm

# Internal packages
# import Fitting_v4 as ft
# import workers_cells_v3 as wkc


def processImage(file_name, current_param, current_log):
    current_log.report("Analysing file: {}\n".format(file_name))

    # creates image object
    im_obj = Image()

    # loads image
    im_obj.load_image(file_name)

    im_obj.z_projection_range(current_param, current_log)

    if im_obj.file_name:
        im_obj.print_image_properties()

    im_obj.show_image(save=False)

    del im_obj


if __name__ == "__main__":

    # defines folders
    root_folder = "/home/marcnol/Documents/Images"
    log_file = "testImageProcess.log"

    # sets parameters
    current_param = parameters.Parameters()
    current_param.initialize_standard_parameters()
    param_file = root_folder + os.sep + "infoList.inf"
    current_param.load_parameters_file(param_file)

    # setup logs
    logFileName = root_folder + os.sep + log_file
    current_log = Log(logFileName)
    current_log.erase_file()

    # processes folders and files
    data_folder = Folders(root_folder)
    data_folder.set_folders()
    current_log.report("folders read: {}".format(len(data_folder.list_folders)))
    files_folder = glob.glob(data_folder.list_folders[0] + os.sep + "*.tif")
    current_log.report("About to read {} files\n".format(len(files_folder)))

    # Processes all DAPI masks
    for file_name in files_folder:
        if file_name.split("_")[-1].split(".")[0] == "ch00" and "DAPI" in file_name.split("_"):
            processImage(file_name, current_param, current_log)

    # Processes all DAPI masks
    for file_name in files_folder:
        res = [i for i in file_name.split("_") if "RT" in i]
        if len(res) > 0 and file_name.split("_")[-1].split(".")[0] == "ch00":
            barcodeNumber = res[0].split("RT")[1]
            current_log.report("Processing barcode: {}".format(barcodeNumber))
            processImage(file_name, current_param, current_log)

    # exits
    current_log.report("Normal exit.")
