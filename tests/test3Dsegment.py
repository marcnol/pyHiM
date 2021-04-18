#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:05:02 2021

@author: marcnol

- tests several algorithms for 3D segmentation:
    - based on ASTROPY (plane by plane object-based segmentation and deblending)

"""

import os, argparse, sys, glob
from datetime import datetime
import matplotlib.pylab as plt
from skimage import io
import numpy as np
from tifffile import imsave
from tqdm import tqdm, trange
from skimage import exposure
from imageProcessing.imageProcessing  import (
    _removesInhomogeneousBackground2D,
    _removesInhomogeneousBackground,
    imageAdjust,
    _segments3DvolumesByThresholding,
    savesImageAsBlocks,
    display3D,
    combinesBlocksImageByReprojection,
    display3D_assembled,
    appliesXYshift3Dimages,
    )
from photutils import detect_sources
from astropy.convolution import Gaussian2DKernel
from skimage import measure
from astropy.visualization import SqrtStretch, simple_norm
from skimage.util.shape import view_as_blocks

from scipy import ndimage as ndi

from skimage.segmentation import watershed
from skimage.feature import peak_local_max


#%% This first test will load a pre-processed 3D image and will run _segments3DvolumesByThresholding to
# find the objects in 3D.

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

# decompose output labeled image in blocks

fileName = rootFolder + "blockTest"
savesImageAsBlocks(segmentedImage3D,fileName,blockSizeXY=256)

# Deblend image in 3D by watersheding

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

#%% loads and display output of a previous segmentation run

folder = "/home/marcnol/grey/users/marcnol/models/StarDist3D/training3Dbarcodes/dataset1/"
file_image3D=folder+os.sep+"scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0.tif"

file_segmented=folder+os.sep+"scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0_segmented.tif"
image3D = io.imread(file_image3D).squeeze()
labeled = io.imread(file_segmented).squeeze()

display3D(image3D=image3D,labels=labeled,z=40, rangeXY=1000)

