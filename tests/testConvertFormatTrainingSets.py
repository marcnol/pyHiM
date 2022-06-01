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

root_folder = "/mnt/PALM_dataserv/DATA/JB/JB/Sara/Deep_Learning/Training_data/Embryo/Marcelo_embryo_data/DAPI_nuclei_segmentation/stage_14"

folder_masks = "Labeled_images"
folder_images = "Original_images"

list_masks, list_images = [], []
list_masks = sorted(glob(root_folder + os.sep + folder_masks + os.sep + "*.tif"))

base_name_masks = [os.path.basename(basename) for basename in list_masks]


for target in base_name_masks:

    expected_folder = root_folder + os.sep + folder_images + os.sep + target.split(".tif")[0]
    expected_folder_45 = root_folder + os.sep + folder_images + os.sep + target.split("_45.tif")[0]

    if os.path.exists(expected_folder):
        file_name = expected_folder + os.sep + "00_Raw_Embryo_segmentation.mat"
        if os.path.exists(file_name):
            list_images.append(file_name)

    elif os.path.exists(expected_folder_45):
        file_name = expected_folder_45 + os.sep + "00_Raw_Embryo_segmentation_45.mat"
        print("blabla")
        if os.path.exists(file_name):
            list_images.append(file_name)

print("Number of Masks: {}".format(len(list_masks)))
print("Number of Images: {}".format(len(list_images)))

if len(list_masks) == len(list_images):

    Y = list(map(imread, list_masks))
    x_mat = list(map(spio.loadmat, list_images))
    # for y in Ymat:
    #     if 'im_raw' in y.keys():
    #         Y.append(y['im_raw'])
    #     elif 'im_raw_45' in y.keys():
    #         Y.append(y['im_raw_45'])

    X = [x[list(x.keys())[-1]] for x in x_mat]

else:
    print("Warning, something is wrong...")
# mat = spio.loadmat(root_folder+file_name, squeeze_me=True)
