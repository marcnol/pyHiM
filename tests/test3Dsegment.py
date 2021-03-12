#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:05:02 2021

@author: marcnol
"""

import os, argparse, sys, glob
from datetime import datetime
import matplotlib.pylab as plt
from skimage import io
import numpy as np
from tifffile import imsave
from tqdm import tqdm, trange
from skimage import exposure,color
from imageProcessing.imageProcessing  import (
    _removesInhomogeneousBackground2D,
    imageAdjust,
    _segments3DvolumesByThresholding,
    savesImageAsBlocks,
    )
from photutils import detect_sources
from astropy.convolution import Gaussian2DKernel
from skimage import measure
from astropy.visualization import SqrtStretch, simple_norm
from skimage.util.shape import view_as_blocks

def display3D(image3D = None,labels=None, z=40, rangeXY=1000, norm=True):


    if image3D is not None:
        images = list()        
        images.append(image3D[z,:,:])
        images.append(image3D[:,rangeXY,:])
        images.append(image3D[:,:,rangeXY])
    else:
        images=[1,1,1]
        
    if labels is not None:
        segmented = list()        
        segmented.append(labels[z,:,:])
        segmented.append(labels[:,rangeXY,:])
        segmented.append(labels[:,:,rangeXY])
    else:
        segmented=[1,1,1]
    percent=99.5
    
    fig, axes = plt.subplots(1, len(images))
    fig.set_size_inches(len(images) * 50, 50)
    ax = axes.ravel()
        
    for image,segm,axis in zip(images,segmented,ax):
        if image3D is not None:
            if norm:
                norm = simple_norm(image, "sqrt", percent=percent)
                axis.imshow(image, cmap="Greys", origin="lower", norm=norm)
            else:
                axis.imshow(image, cmap="Greys", origin="lower")
        if labels is not None:
            axis.imshow(color.label2rgb(segm, bg_label=0),alpha=.3)
    
#%% loads and segments a file

rootFolder="/mnt/grey/DATA/users/marcnol/models/StarDist3D/training3Dbarcodes/dataset1/"
file = rootFolder+'scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0.tif'

threshold_over_std,sigma ,boxSize, filter_size=1, 3, (32, 32),(3, 3)
nlevels=64
contrast=0.001
image3D = io.imread(file).squeeze()

binary, segmentedImage3D = _segments3DvolumesByThresholding(image3D,
                                                       threshold_over_std=threshold_over_std, 
                                                       sigma = 3, 
                                                       boxSize=(32, 32),
                                                       filter_size=(3, 3),
                                                       nlevels=nlevels,
                                                       contrast=contrast,
                                                       deblend3D=True)


display3D(image3D=image3D,labels=segmentedImage3D,z=40, rangeXY=1000)


#%% decompose output labeled image in blocks

fileName = rootFolder + "blockTest"
savesImageAsBlocks(segmentedImage3D,fileName,blockSizeXY=256)
    
#%% Deblend image in 3D by watersheding
    
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi

from skimage.segmentation import watershed
from skimage.feature import peak_local_max

binary0=binary>0

# Now we want to separate the two objects in image
# Generate the markers as local maxima of the distance to the background
print("Constructing distance matrix from 3D binary mask...")
distance = ndi.distance_transform_edt(binary0)
coords = peak_local_max(distance, footprint=np.ones((10, 10, 25)), labels=binary0)
mask = np.zeros(distance.shape, dtype=bool)
mask[tuple(coords.T)] = True
markers, _ = ndi.label(mask)
print("Deblending sources in 3D by watersheding...")
labels = watershed(-distance, markers, mask=binary0)

display3D(image3D=image3D,labels=labels,z=40, rangeXY=1000)


#%% loads and display output of a segmentation run

folder = "/home/marcnol/grey/users/marcnol/models/StarDist3D/training3Dbarcodes/dataset1/"
file_image3D=folder+os.sep+"scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0.tif"

file_segmented=folder+os.sep+"scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0_segmented.tif"
image3D = io.imread(file_image3D).squeeze()
labeled = io.imread(file_segmented).squeeze()

display3D(image3D=image3D,labels=labeled,z=40, rangeXY=1000)

#%%

display3D(image3D=image3D,z=40, rangeXY=1000,norm=False)