#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 09:17:56 2021

@author: marcnol
"""


from skimage import io
import os
import numpy as np
import matplotlib.pylab as plt
from imageProcessing.imageProcessing import _reinterpolatesFocalPlane, imageShowWithValues, findsFocusFromBlocks

# from astropy.stats import SigmaClip
from scipy.stats import sigmaclip

rootFolder = "/home/marcnol/data/Embryo_debug_dataset/Experiment_18"
filename = "scan_001_RT27_001_ROI_converted_decon_ch00.tif"
fullfilename = rootFolder + os.sep + filename

data = io.imread(fullfilename).squeeze()

blockSizeXY = 128

output, focalPlaneMatrix, zRange, focusPlane, LaplacianMeans = _reinterpolatesFocalPlane(data, blockSizeXY)

focus, filteredMatrix = findsFocusFromBlocks(focalPlaneMatrix, LaplacianMeans, threshold=0.1)
print("Consensus focal plane: {}".format(focus))
imageShowWithValues([focalPlaneMatrix, filteredMatrix], verbose=True, title="focal plane = " + "{:.2f}".format(focus))
#%%
