#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 10:38:40 2021

@author: marcnol

- tests several algorithms for 3D segmentation:
    - based on ASTROPY (plane by plane object-based segmentation and deblending)
    - blob_log from skimage


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


#%% Test 3D localization of centroids using skimage: this is what starFISH currently do for 3D gaussian fitting
# it does not work very well. It is biased strongly by the shape of the spots.

from skimage import data, feature, exposure
from skimage import io

z = 40

rootFolder="/mnt/grey/DATA/users/marcnol/models/StarDist3D/training3Dbarcodes/dataset1/"
file = rootFolder+'scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0.tif'

print("loading image {}".format(os.path.basename(file)))
image3D = io.imread(file).squeeze()

print("Normalizing exposures")
image3D = exposure.equalize_hist(image3D)  # improves detection

print("Calling blob_log to detect in 3D")
# localizationTable = feature.blob_log(image3D, threshold = .3)
localizationTable = feature.blob_dog(image3D[z-1:z+2,:,:], threshold = .3)

display3D(image3D=image3D,localizations = localizationTable, z=40, rangeXY=1000,norm=False)

