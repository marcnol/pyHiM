#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 16:04:23 2020

@author: marcnol
"""

import os, glob
from fileProcessing.fileManagement import Folders, Parameters

# function to get unique values
def unique(list1):

    # intilize a null list
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)

    return unique_list


root_folder = "./"
all_files_in_root_folder = glob.glob(root_folder + "*tif")
print("Number of files: {}".format(len(all_files_in_root_folder)))

parameter_file = "infoList_DAPI.json"

# data_folder = Folders(root_folder)
current_param = Parameters(root_folder, root_folder + parameter_file)

rois, rts = [], []
for x in all_files_in_root_folder:
    file_parts = current_param.decode_file_parts(x)
    rois.append(file_parts["roi"])
    rts.append(file_parts["cycle"])

number_unique_cycles = len(unique(rts))

print("Number Unique cycles: {}".format(number_unique_cycles))
