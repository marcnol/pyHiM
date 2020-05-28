#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 21:06:41 2020

@author: marcnol
"""


import numpy as np
import matplotlib.pyplot as plt
import os
from skimage import io

from imageProcessing import Image, align2ImagesCrossCorrelation, save2imagesRGB
from fileManagement import folders, Parameters, session
from fileManagement import log
from alignImages import align2Files
from scipy.ndimage import shift as shiftImage

import astroalign as aa
import imreg

rootFolder = "/home/marcnol/data/Experiment_15/Embryo_006_ROI18/rawData/zProject/"
fileName = rootFolder + "scan_001_DAPI_018_ROI_converted_decon_ch01_2d.npy"
fileNameReference = rootFolder + "scan_004_RT18_018_ROI_converted_decon_ch00_2d.npy"
"""
#rootFolder='/home/marcnol/Documents/Images/Experiment15_embryo001/rawImages'
#fileName1='scan_001_RT14_002_ROI_converted_decon_ch01.tif'
#fileName2='scan_001_RT1_002_ROI_converted_decon_ch01.tif'
rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset/rawImages'
fileName1='scan_001_DAPI_017_ROI_converted_decon_ch01.tif'
#fileName1='scan_004_RT18_017_ROI_converted_decon_ch00.tif'
fileName2='scan_004_RT20_017_ROI_converted_decon_ch00.tif'
"""
outputFileName = "./test/" + os.path.basename(fileName).split(".")[0]

# loads images
image1 = np.load(fileNameReference).squeeze()
image2 = np.load(fileName).squeeze()

image1_uncorrected = image1 / image1.max()
image2_uncorrected = image2 / image2.max()

#%% This is my method from MATLAB. Rescales images and performs CC and then refits max of CC peak
(
    shift,
    error,
    diffphase,
    lower_threshold,
    I_histogram,
    image2_corrected,
    image1_adjusted,
    image2_adjusted,
) = align2ImagesCrossCorrelation(image1_uncorrected, image2_uncorrected)

image2_corrected_raw = shiftImage(image2_uncorrected, shift)
image2_corrected_raw[image2_corrected_raw < 0] = 0

# thresholds corrected images for better display and saves
image1_corrected = image1_adjusted > 0.1
image2_corrected = image2_corrected > 0.1
save2imagesRGB(image1_uncorrected, image2_corrected_raw, outputFileName + "_overlay_corrected.png")

#%% image_registration package from astropy
""" 
this seems to give the same result as the routine I am currently using
the only nice thing is that it provides errors of alignment
"""

from image_registration import chi2_shift

# import image_registration

# on raw images. There is a 1px difference with respect to non-adjusted images
dx, dy, edx, edy = chi2_shift(image1_uncorrected, image2_uncorrected, upsample_factor="auto")

shift1 = -np.array((dy, dx))
image2_corrected = shiftImage(image2_adjusted, shift1)
image2_corrected[image2_corrected < 0] = 0
image2_corrected /= image2_corrected.max()
save2imagesRGB(image1_uncorrected, image2_corrected, outputFileName + "_overlay_corrected_chi2.png")

# on adjusted images
dx2, dy2, edx2, edy2 = chi2_shift(image1_adjusted, image2_adjusted, upsample_factor="auto")

shift2 = -np.array((dy2, dx2))


#%% uses imreg_dft
"""
very slow. Does rotations and scaling. Results are not what I get using other methods. Documentation
is crappy so not clear what it is doing and how. Typically from Enrico...
"""
import imreg_dft as ird

result = ird.similarity(image1_uncorrected, image2_uncorrected, numiter=5)
ird.imshow(image1_uncorrected, image2_uncorrected, result["timg"])
plt.show()
image2_corrected = result["timg"]
image2_corrected[image2_corrected < 0] = 0
image2_corrected /= image2_corrected.max()
# thresholds corrected images for better display and saves
save2imagesRGB(
    image1_uncorrected, image2_corrected, outputFileName + "_overlay_corrected_imregdft.png",
)

#%% uses imreg
import imreg

im2, scale, angle, (t0, t1) = imreg.similarity(image1_uncorrected, image2_uncorrected)
t0, t1 = imreg.translation(image1_uncorrected, image2_uncorrected)

plt.imshow(image1_uncorrected, image2_uncorrected, im2)

#%% uses astroalign
import astroalign as aa

p, (pos_img, pos_img_rot) = aa.find_transform(image1_adjusted, image2_adjusted)

#%%
img_aligned, footprint = aa.register(image1_uncorrected, image2_uncorrected)

fig, axes = plt.subplots(2, 2, figsize=(10, 10))
axes[0, 0].imshow(img, cmap="gray", interpolation="none", origin="lower")
axes[0, 0].axis("off")
axes[0, 0].set_title("Source Image")

axes[0, 1].imshow(img_rotated, cmap="gray", interpolation="none", origin="lower")
axes[0, 1].axis("off")
axes[0, 1].set_title("Target Image")

axes[1, 1].imshow(img_aligned, cmap="gray", interpolation="none", origin="lower")
axes[1, 1].axis("off")
axes[1, 1].set_title("Source Image aligned with Target")

axes[1, 0].axis("off")

plt.tight_layout()
plt.show()
