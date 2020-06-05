#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 21:06:41 2020

@author: marcnol
"""


import numpy as np
import matplotlib.pyplot as plt
import os

from skimage import data, filters
from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft
from skimage import exposure
from scipy.ndimage import fourier_shift
from scipy.ndimage import shift as shiftImage
from imageProcessing import Image
from fileManagement import folders, Parameters
from fileManagement import log
import cv2


def displaysEquializationHistograms(
    hist1_before, hist1_after, hist2_before, hist2_after, min_threshold, vebose=False, fileName="test",
):
    fig = plt.figure(figsize=(6, 3))
    ax1 = plt.subplot(2, 2, 1)
    ax2 = plt.subplot(2, 2, 2)
    ax3 = plt.subplot(2, 2, 3)
    ax4 = plt.subplot(2, 2, 4)

    ax1.plot(hist1_before[1], hist1_before[0])
    ax2.plot(hist2_before[1], hist2_before[0])
    ax3.plot(hist1_before[1], hist1_before[0])
    ax4.plot(hist2_before[1], hist2_before[0])
    ax3.set_yscale("log")
    ax4.set_yscale("log")
    ax1.vlines(min_threshold[0], 0, hist1_before[0].max(), colors="r")
    ax2.vlines(min_threshold[1], 0, hist2_before[0].max(), colors="r")
    plt.savefig(fileName + "_intensityHist.png")
    if not verbose:
        plt.close(fig)


def imageAdjust(image, fileName="test", lower_threshold=0.3, higher_threshold=0.9999, display=False):

    image1 = exposure.rescale_intensity(image, out_range=(0, 1))

    hist1_before = exposure.histogram(image1)

    image1 = exposure.rescale_intensity(image1, in_range=(lower_threshold, higher_threshold), out_range=(0, 1))

    hist1 = exposure.histogram(image1)

    if display:
        plt.figure(figsize=(30, 30))
        plt.imsave(fileName + "_adjusted.png", image1, cmap="hot")

    return image1, hist1_before, hist1


def showCCimage(image1, image2, verbose=False, fileName="test"):
    image_product = np.fft.fft2(image1) * np.fft.fft2(image2).conj()
    cc_image = _upsampled_dft(image_product, 150, 100, (shift * 100) + 75).conj()
    if verbose:
        plt.figure(figsize=(60, 30))
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
        RGB_falsecolor_image = np.dstack([image1, image2, np.zeros([2048, 2048])])
        ax1.imshow(RGB_falsecolor_image, origin="lower", interpolation="nearest")
        ax1.set_axis_off()
        ax1.set_title("Super-imposed images")

        ax2.imshow(cc_image.real)
        ax2.set_axis_off()
        ax2.set_title("Supersampled XC sub-area")
    else:
        plt.figure(figsize=(30, 30))
        plt.imsave(fileName + "_CC.png", cc_image)


def save2imagesRGB(image1, image2, fileName):
    RGB_falsecolor_image = np.dstack([image1, np.zeros([2048, 2048]), image2])
    plt.figure(figsize=(30, 30))
    plt.imsave(fileName, RGB_falsecolor_image)
    # cv2.write(fileName+'_RGB.png',RGB_falsecolor_image)


rootFolder = "/home/marcnol/Documents/Images/Embryo_debug_dataset/rawImages"
fileName1 = "scan_001_DAPI_017_ROI_converted_decon_ch01.tif"
# fileName1='scan_004_RT18_017_ROI_converted_decon_ch00.tif'
fileName2 = "scan_004_RT20_017_ROI_converted_decon_ch00.tif"

rootFolder = "/home/marcnol/data/Experiment_15/Embryo_006_ROI18"
fileName2 = "scan_001_DAPI_018_ROI_converted_decon_ch01.tif"
fileName1 = "scan_004_RT18_018_ROI_converted_decon_ch00.tif"

# rootFolder='/home/marcnol/Documents/Images/Experiment15_embryo001/rawImages'
# fileName1='scan_001_RT14_002_ROI_converted_decon_ch01.tif'
# fileName2='scan_001_RT1_002_ROI_converted_decon_ch01.tif'

verbose = False

outputFileName = os.path.basename(fileName2).split(".")[0]

parameterFiles = [
    "infoList_DAPI.json",
    "infoList_barcode.json",
    "infoList_fiducial.json",
]

label = parameterFiles[0]
param = Parameters()
param.initializeStandardParameters()
paramFile = rootFolder + os.sep + label
param.loadParametersFile(paramFile)
param.param["rootFolder"] = rootFolder

logFile = "alignImagesXcorrelation.log"
logFileName = rootFolder + os.sep + logFile
log1 = log(logFileName)
log1.eraseFile()
log1.report("Starting to log to: {}".format(logFileName))

# processes folders and files
dataFolder = folders(param.param["rootFolder"])
dataFolder.setsFolders()
dataFolder.createsFolders(rootFolder, param)
log1.report("folders read: {}".format(len(dataFolder.listFolders)))

# creates image object
Im1 = Image()
Im2 = Image()

# loads image
Im1.loadImage2D(fileName1, log1, dataFolder)
Im2.loadImage2D(fileName2, log1, dataFolder)

# adjusts image levels
image1_uncorrected = Im1.data_2D / Im1.data_2D.max()
image2_uncorrected = Im2.data_2D / Im2.data_2D.max()
lower_threshold1 = 2 * filters.threshold_otsu(image1_uncorrected)
lower_threshold2 = 2 * filters.threshold_otsu(image2_uncorrected)

image1, hist1_before, hist1_after = imageAdjust(
    image1_uncorrected, outputFileName + "_ref", lower_threshold1, higher_threshold=0.9999, display=verbose,
)
image2, hist2_before, hist2_after = imageAdjust(
    image2_uncorrected, outputFileName, lower_threshold2, higher_threshold=0.9999, display=verbose,
)

# shows eq histograms
displaysEquializationHistograms(
    hist1_before, hist1_after, hist2_before, hist2_after, (lower_threshold1, lower_threshold2), verbose, outputFileName,
)

# calculates shift
shift, error, diffphase = register_translation(image1, image2, 100)

# corrects image
# shift = (-22.4, 5.2)
# The shift corresponds to the pixel offset relative to the reference image
# offset_image = fourier_shift(np.fft.fftn(image2_uncorrected), shift)
# image2_corrected = np.fft.ifftn(offset_image)
image2_corrected = shiftImage(image2, shift)
image2_corrected = exposure.rescale_intensity(image2_corrected, out_range=(0, 1))

print(f"Detected subpixel offset (y, x): {shift}")


if verbose:
    showCCimage(image1, image2, verbose, outputFileName)

save2imagesRGB(image1_uncorrected, image2_uncorrected, outputFileName + "_overlay_uncorrected.png")
image1 = image1 > 0.1
image2_corrected = image2_corrected > 0.1
save2imagesRGB(image1, image2_corrected, outputFileName + "_overlay_corrected.png")
