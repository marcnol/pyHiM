#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 16:04:23 2020

@author: marcnol
"""

import os, glob
from fileProcessing.fileManagement import folders, Parameters

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


rootFolder = "./"
allFilesinRootFolder = glob.glob(rootFolder + "*tif")
print("Number of files: {}".format(len(allFilesinRootFolder)))

parameterFile = "infoList_DAPI.json"

# dataFolder = folders(rootFolder)
param = Parameters(rootFolder, rootFolder + parameterFile)

ROIs, RTs = [], []
for x in allFilesinRootFolder:
    fileParts = param.decodesFileParts(x)
    ROIs.append(fileParts["roi"])
    RTs.append(fileParts["cycle"])

numberUniqueCycles = len(unique(RTs))

print("Number Unique cycles: {}".format(numberUniqueCycles))
