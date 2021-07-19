#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 17:10:17 2020

@author: marcnol
"""

import numpy as np
import os
import glob
from skimage import exposure,color

from mayavi.mlab import contour3d,pipeline

from matplotlib.pylab import plt

# from imageProcessing.imageProcessing import Image
# from fileProcessing.fileManagement import Parameters, log
from skimage import io

#%% loads data for only one barcode and one DAPI

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

fullFileName = rootFolder + os.sep + fileName
fullFileNameF = rootFolder + os.sep + fileNameF

data = [io.imread(fullFileName).squeeze()]
data_DAPI = [io.imread(fullFileNameF).squeeze()]

#%% selects all files in a folder within a given ROI
rootFolder = "/home/marcnol/grey/users/marcnol/test_HiM/run_zBinning2"
rootFolder = "/mnt/PALM_dataserv/DATA/Olivier/Thesis/Insulators_Project/Paper_Insulator/HiM_analysis/Doc_Locus/Embryo_001"

ROI="003"

fileNames_RT = [x for x in glob.glob(rootFolder+os.sep+"*tif") if ROI in os.path.basename(x).split("_")[3]\
             and "DAPI" not in x and "ch01" in x]
fileNames_DAPI = [x for x in glob.glob(rootFolder+os.sep+"*tif") if ROI in os.path.basename(x).split("_")[3]\
             and "DAPI" in x and "ch00" in x]

print("\nRT found = {}".format(fileNames_RT))
print("\nDAPI found = {}".format(fileNames_DAPI))

data = [io.imread(x).squeeze() for x in fileNames_RT]
data_DAPI = [io.imread(x).squeeze() for x in fileNames_DAPI]

#%% preprocesses

xyRange = (800,1200)
zRange = (10,50)
zShift = 8

data = [exposure.rescale_intensity(x, out_range=(0, 1)) for x in data]
data_DAPI = [exposure.rescale_intensity(x, out_range=(0, 1)) for x in data_DAPI]

subdata = [x[zShift+zRange[0]:zShift+zRange[1], xyRange[0]:xyRange[1], xyRange[0]:xyRange[1]] for x in data]
subdata_DAPI = [x[zRange[0]:zRange[1], xyRange[0]:xyRange[1], xyRange[0]:xyRange[1]] for x in data_DAPI]

#%% Displays overlays between two channels
preFactor = 1

f0ig, (ax1, ax2) = plt.subplots(2, 1)
# ax1.imshow(np.sum(subdata_DAPI[0][:,450:550,:]*1,axis=1), alpha=.7, cmap='Blues')
ax2.imshow(np.sum(subdata_DAPI[0][:,:,:]*1,axis=0), alpha=.7, cmap='Blues')

sum_RT_zx = 0.0*np.sum(subdata[0][:,450:550,:],axis=1)
sum_RT_xy = 0.0*np.sum(subdata[0][:,:,:]*preFactor,axis=0)
for x in subdata:
    sum_RT_zx = sum_RT_zx + np.sum(x[:,450:550,:],axis=1)
    sum_RT_xy = sum_RT_xy + np.sum(x[:,:,:]*preFactor,axis=0)

cmap_RT = 'RdBu'
ax1.imshow(sum_RT_zx,alpha=0.5,cmap=cmap_RT)
ax2.imshow(sum_RT_xy,alpha=0.5,cmap=cmap_RT)

# fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
# ax1.imshow(data[:, :, 500]*100)
# ax2.imshow(dataF[:, :, 500] * 10)
# ax3.imshow(dataF[35,:,:])

#%% Displays 3D level reconstruction of DAPI using mayavi

# contour3d(subdata_DAPI[0], vmin=0.2, vmax=1, line_width = 0, colormap = 'Pastel2')
sf = pipeline.scalar_field(subdata_DAPI[0])
iso = pipeline.iso_surface(sf,contours=[0.1],line_width = 0, colormap = 'spectral',opacity=.9)
iso.actor.property.representation = 'wireframe'

# Displays 3D level reconstruction of barcode using mayavi
cmaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
barcodes2plot=[0,1,2,3,4,5,6,7,8,9,10]
subdata2Plot = [x for i,x in enumerate(subdata) if i in barcodes2plot]
for x,cmap in zip(subdata2Plot,cmaps[0:len(subdata2Plot)]):
    contour3d(x, vmin=0, vmax=1,colormap=cmap,contours=[0.05])


