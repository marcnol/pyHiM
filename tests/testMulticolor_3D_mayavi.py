#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 17:10:17 2020

@author: marcnol

this test script is intended to plot two or three 3D images using mayavi. The file names are provided manually at the beginning of the script
"""

import numpy as np
import os
import glob
from skimage import exposure,color

from mayavi.mlab import contour3d,pipeline

from matplotlib.pylab import plt

from skimage import io

#%% loads data for only one barcode and one DAPI

rootFolder = "/mnt/grey/DATA/rawData_2019/Experiment_17/2019-08-08/Embryos_deconvolved/Embryo_001"

fileNames = ["scan_002_DAPI_000_ROI_converted_decon_ch00.tif","scan_002_RT1_000_ROI_converted_decon_ch00.tif"]

fullFileNames = [rootFolder + os.sep + fileName for fileName in fileNames]

data = [io.imread(x).squeeze() for x in fullFileNames]

#%% rescales images
data = [exposure.rescale_intensity(x, out_range=(0, 1)) for x in data]

#%%
# takes ROI
xyRange = (900,1300) 
zRange = (10,50)
zShift = 8

subdata = [x[zShift+zRange[0]:zShift+zRange[1], xyRange[0]:xyRange[1], xyRange[0]:xyRange[1]] for x in data]

#%% Displays overlays between two channels
preFactor = 1

fig, (ax1, ax2) = plt.subplots(2, 1)
ax2.imshow(np.sum(subdata[0][:,:,:]*1,axis=0), alpha=.7, cmap='terrain')

sum_RT_zx = 0.0*np.sum(subdata[0][:,450:550,:],axis=1)
sum_RT_xy = 0.0*np.sum(subdata[0][:,:,:]*preFactor,axis=0)
for x in subdata[1:]:
    sum_RT_zx = sum_RT_zx + np.sum(x[:,450:550,:],axis=1)
    sum_RT_xy = sum_RT_xy + np.sum(x[:,:,:]*preFactor,axis=0)
    
cmap_RT = 'RdBu'
ax1.imshow(sum_RT_zx,alpha=0.5,cmap=cmap_RT)
ax2.imshow(sum_RT_xy,alpha=0.5,cmap=cmap_RT)


#%% Displays 3D level reconstruction of DAPI using mayavi

# contour3d(subdata_DAPI[0], vmin=0.2, vmax=1, line_width = 0, colormap = 'Pastel2')
sf = pipeline.scalar_field(subdata[0])
iso = pipeline.iso_surface(sf,contours=[0.005],line_width = 0, colormap = 'spectral',opacity=.9)
iso.actor.property.representation = 'wireframe'

# Displays 3D level reconstruction of barcode using mayavi
cmaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
barcodes2plot=[1]
subdata2Plot = [x for i,x in enumerate(subdata) if i in barcodes2plot]
for x,cmap in zip(subdata2Plot,cmaps[0:len(subdata2Plot)]):
    contour3d(x, vmin=0, vmax=1,colormap=cmap,contours=[0.05])


