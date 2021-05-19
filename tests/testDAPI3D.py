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

# from imageProcessing.imageProcessing import Image
# from fileProcessing.fileManagement import Parameters, log
from skimage import io

#%% loads data

rootFolder = "/home/marcnol/Downloads"

rootFolder = "/home/marcnol/grey/rawData_2020/Exp_Combinatory_3_Tof/ROIs/ROI001"
fileName = "scan_007_DAPI_001_ROI_converted_decon_ch00.tif"
fileNameF = "scan_007_DAPI_001_ROI_converted_decon_ch01.tif"

rootFolder = "/home/marcnol/data/Embryo_debug_dataset/test_dataset"
fileName = "scan_001_RT27_001_ROI_converted_decon_ch01.tif"
fileNameF="scan_006_DAPI_001_ROI_converted_decon_ch00.tif"

rootFolder = "/home/marcnol/data/Embryo_debug_dataset/test_dataset"
fileName = "scan_001_RT27_001_ROI_converted_decon_ch01.tif"
fileNameF="scan_006_DAPI_001_ROI_converted_decon_ch00.tif"

rootFolder = "/home/marcnol/grey/users/marcnol/test_HiM/run_zBinning2"
fileName = "scan_001_RT41_003_ROI_converted_decon_ch01.tif"
fileNameF = "scan_006_DAPI_003_ROI_converted_decon_ch00.tif"

fullFileName = rootFolder + os.sep + fileName
fullFileNameF = rootFolder + os.sep + fileNameF

data = io.imread(fullFileName).squeeze()
dataF = io.imread(fullFileNameF).squeeze()


#%% Displays overlays between two channels

subdata = data[10:50, 0:2000, 0:2000]
subdataF = dataF[10:50, 0:2000, 0:2000]

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.imshow(np.sum(data[:,450:550,:]*10,axis=1),alpha=0.7,cmap='Reds')
ax1.imshow(np.sum(dataF[:,450:550,:]*1,axis=1), alpha=.7, cmap='Blues')

ax2.imshow(np.sum(data[:,:,:]*10,axis=0),alpha=0.7,cmap='Reds')
ax2.imshow(np.sum(dataF[:,:,:]*1,axis=0), alpha=.7, cmap='Blues')

# fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
# ax1.imshow(data[:, :, 500]*100)
# ax2.imshow(dataF[:, :, 500] * 10)
# ax3.imshow(dataF[35,:,:])

#%% Displays XY and YZ projections of the 3D volume


#%% Displays XY and YZ abd XZ slices of the 3D volume


#%% Displays 3D level reconstruction of barcode using mayavi

contour3d(subdata, vmin=0, vmax=25000)

# Displays 3D level reconstruction of DAPI using mayavi

contour3d(subdataF, vmin=0, vmax=25000)
