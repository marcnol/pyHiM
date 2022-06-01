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
from imageProcessing.imageProcessing import _reinterpolate_focal_plane, image_show_with_values, findsFocusFromBlocks

# from astropy.stats import SigmaClip
from scipy.stats import sigmaclip

root_folder = "/home/marcnol/data/Embryo_debug_dataset/Experiment_18"
filename = "scan_001_RT27_001_ROI_converted_decon_ch00.tif"
fullfilename = root_folder + os.sep + filename

data = io.imread(fullfilename).squeeze()

block_size_xy = 128

output, focal_plane_matrix, z_range, focus_plane, LaplacianMeans = _reinterpolate_focal_plane(data, block_size_xy)

focus, filteredMatrix = findsFocusFromBlocks(focal_plane_matrix, LaplacianMeans, threshold=0.1)
print("Consensus focal plane: {}".format(focus))
image_show_with_values([focal_plane_matrix, filteredMatrix], verbose=True, title="focal plane = " + "{:.2f}".format(focus))
#%%
