#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 15:45:56 2020

@author: marcnol

testLoadsROIs

"""

import os
import numpy as np
from matplotlib import pyplot as plt

from skimage import measure
from imageProcessing import Image

os.chdir("/home/marcnol/Repositories/pyHiM")
masks = np.load("masks.npy")


# Create image
rootFolder = "/home/marcnol/data/Experiment_4/0_Embryo/alignImages/"
fileNameRNA = rootFolder + "scan_002_DAPI_001_ROI_converted_decon_ch01_2d_registered.npy"
img = np.load(fileNameRNA).squeeze()
Im = Image()
Im.data_2D = img
ax = Im.imageShow(show=True)

# fig, ax = plt.subplots()
allMasks = np.sum(masks, axis=2)

ax.imshow(allMasks, cmap=plt.cm.gray)

for imask in range(masks.shape[2]):
    # ax=plt.imshow(masks[:,:,imask], interpolation='nearest', cmap="Greys")
    print("shown mask {}".format(imask))

    contours = measure.find_contours(masks[:, :, imask], 0.8)

    for n, contour in enumerate(contours):
        ax.plot(contour[:, 1], contour[:, 0], linewidth=2)

plt.show()


#%%

from astropy.table import Table, Column, vstack


ROI = 1
fileNameDAPImask = "./segmentedObjects/scan_002_DAPI_001_ROI_converted_decon_ch00_Masks.npy"
fileName = "./segmentedObjects/scan_002_DAPI_001_ROI_converted_decon_ch01_2d_registered_SNDmask_sna.npy"
# load DAPI mask
maskDAPI = Image()
maskDAPI.data_2D = np.load(fileNameDAPImask).squeeze()

# load SND mask
maskSND = Image()
maskSND.data_2D = np.load(fileName).squeeze()

# matches cells and SND masks
newMatrix = maskDAPI.data_2D * maskSND.data_2D
cellsWithinMask = np.unique(newMatrix)
allResults = Table()

if len(cellsWithinMask) > 0:
    newTable = Table()
    newTable["ROI #"] = int(ROI) * np.ones(len(cellsWithinMask))
    newTable["CellID #"] = cellsWithinMask
    allResults = vstack([allResults, newTable])
