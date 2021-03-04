#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 17:10:17 2020

@author: marcnol
"""

import numpy as np
import os

from mayavi.mlab import contour3d
from matplotlib.pylab import plt

from imageProcessing.imageProcessing import Image
from fileProcessing.fileManagement import Parameters, log
from skimage import io

#%% loads data

rootFolder = "/home/marcnol/Downloads"
fileName = "scan_006_DAPI_001_ROI_converted_decon_ch00.tif"

rootFolder = "/home/marcnol/grey/rawData_2020/Exp_Combinatory_3_Tof/ROIs/ROI001"
fileName = "scan_007_DAPI_001_ROI_converted_decon_ch00.tif"
fullFileName = rootFolder + os.sep + fileName

fileNameF = "scan_007_DAPI_001_ROI_converted_decon_ch01.tif"
fullFileNameF = rootFolder + os.sep + fileNameF

data = io.imread(fullFileName).squeeze()
dataF = io.imread(fullFileNameF).squeeze()

subdata = data[22:45, 0:2000, 0:2000]


#%% Displays overlays between two channels


# fig, ax=plt.subplots()
# ax.imshow(data[:,500,:],alpha=0.7,cmap='Reds')
# ax.imshow(dataF[:,500,:]*1, alpha=.7, cmap='Blues')

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.imshow(data[:, :, 500])
ax2.imshow(dataF[:, :, 500] * 10)


#%% Displays XY and YZ projections of the 3D volume


#%% Displays XY and YZ abd XZ slices of the 3D volume


#%% Displays 3D level reconstruction using mayavi

contour3d(subdata, vmin=5000, vmax=40000)

# contour3d(data)
