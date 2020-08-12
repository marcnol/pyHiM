#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 14:28:21 2020

@author: marcnol
"""

import os
from glob import glob

from tifffile import imread

import scipy.io as spio

rootFolder = "/mnt/PALM_dataserv/DATA/JB/JB/Sara/Deep_Learning/Training_data/Embryo/Marcelo_embryo_data/DAPI_nuclei_segmentation/stage_14"

folderMasks = "Labeled_images"
folderImages = "Original_images"

ListMasks, ListImages = [], []
ListMasks = sorted(glob(rootFolder + os.sep + folderMasks + os.sep + "*.tif"))

baseNameMasks = [os.path.basename(basename) for basename in ListMasks]


for target in baseNameMasks:

    expectedFolder = rootFolder + os.sep + folderImages + os.sep + target.split(".tif")[0]
    expectedFolder45 = rootFolder + os.sep + folderImages + os.sep + target.split("_45.tif")[0]

    if os.path.exists(expectedFolder):
        fileName = expectedFolder + os.sep + "00_Raw_Embryo_segmentation.mat"
        if os.path.exists(fileName):
            ListImages.append(fileName)

    elif os.path.exists(expectedFolder45):
        fileName = expectedFolder45 + os.sep + "00_Raw_Embryo_segmentation_45.mat"
        print("blabla")
        if os.path.exists(fileName):
            ListImages.append(fileName)

print("Number of Masks: {}".format(len(ListMasks)))
print("Number of Images: {}".format(len(ListImages)))

if len(ListMasks) == len(ListImages):

    Y = list(map(imread, ListMasks))
    Xmat = list(map(spio.loadmat, ListImages))
    # for y in Ymat:
    #     if 'im_raw' in y.keys():
    #         Y.append(y['im_raw'])
    #     elif 'im_raw_45' in y.keys():
    #         Y.append(y['im_raw_45'])

    X = [x[list(x.keys())[-1]] for x in Xmat]

else:
    print("Warning, something is wrong...")
# mat = spio.loadmat(rootFolder+fileName, squeeze_me=True)
