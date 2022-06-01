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
from scipy.ndimage import shift as shift_image
from imageProcessing import Image
from fileManagement import Folders, Parameters
from fileManagement import Log
import cv2


def displaysEquializationHistograms(
    hist1_before, hist1_after, hist2_before, hist2_after, min_threshold, vebose=False, file_name="test",
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
    plt.savefig(file_name + "_intensityHist.png")
    if not verbose:
        plt.close(fig)


def image_adjust(image, file_name="test", lower_threshold=0.3, higher_threshold=0.9999, display=False):

    image1 = exposure.rescale_intensity(image, out_range=(0, 1))

    hist1_before = exposure.histogram(image1)

    image1 = exposure.rescale_intensity(image1, in_range=(lower_threshold, higher_threshold), out_range=(0, 1))

    hist1 = exposure.histogram(image1)

    if display:
        plt.figure(figsize=(30, 30))
        plt.imsave(file_name + "_adjusted.png", image1, cmap="hot")

    return image1, hist1_before, hist1


def show_cc_image(image1, image2, verbose=False, file_name="test"):
    image_product = np.fft.fft2(image1) * np.fft.fft2(image2).conj()
    cc_image = _upsampled_dft(image_product, 150, 100, (shift * 100) + 75).conj()
    if verbose:
        plt.figure(figsize=(60, 30))
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
        rgb_falsecolor_image = np.dstack([image1, image2, np.zeros([2048, 2048])])
        ax1.imshow(rgb_falsecolor_image, origin="lower", interpolation="nearest")
        ax1.set_axis_off()
        ax1.set_title("Super-imposed images")

        ax2.imshow(cc_image.real)
        ax2.set_axis_off()
        ax2.set_title("Supersampled XC sub-area")
    else:
        plt.figure(figsize=(30, 30))
        plt.imsave(file_name + "_CC.png", cc_image)


def save_2_images_rgb(image1, image2, file_name):
    rgb_falsecolor_image = np.dstack([image1, np.zeros([2048, 2048]), image2])
    plt.figure(figsize=(30, 30))
    plt.imsave(file_name, rgb_falsecolor_image)
    # cv2.write(file_name+'_RGB.png',rgb_falsecolor_image)


root_folder = "/home/marcnol/Documents/Images/Embryo_debug_dataset/raw_images"
filename_1 = "scan_001_DAPI_017_ROI_converted_decon_ch01.tif"
# filename_1='scan_004_RT18_017_ROI_converted_decon_ch00.tif'
filename_2 = "scan_004_RT20_017_ROI_converted_decon_ch00.tif"

root_folder = "/home/marcnol/data/Experiment_15/Embryo_006_ROI18"
filename_2 = "scan_001_DAPI_018_ROI_converted_decon_ch01.tif"
filename_1 = "scan_004_RT18_018_ROI_converted_decon_ch00.tif"

# root_folder='/home/marcnol/Documents/Images/Experiment15_embryo001/raw_images'
# filename_1='scan_001_RT14_002_ROI_converted_decon_ch01.tif'
# filename_2='scan_001_RT1_002_ROI_converted_decon_ch01.tif'

verbose = False

output_filename = os.path.basename(filename_2).split(".")[0]

parameterFiles = [
    "infoList_DAPI.json",
    "infoList_barcode.json",
    "infoList_fiducial.json",
]

label = parameterFiles[0]
current_param = Parameters()
current_param.initialize_standard_parameters()
param_file = root_folder + os.sep + label
current_param.load_parameters_file(param_file)
current_param.param_dict["rootFolder"] = root_folder

log_file = "alignImagesXcorrelation.log"
logFileName = root_folder + os.sep + log_file
current_log = Log(logFileName)
current_log.erase_file()
current_log.report("Starting to log to: {}".format(logFileName))

# processes folders and files
data_folder = Folders(current_param.param_dict["rootFolder"])
data_folder.set_folders()
data_folder.create_folders(root_folder, current_param)
current_log.report("folders read: {}".format(len(data_folder.list_folders)))

# creates image object
img_1 = Image()
img_2 = Image()

# loads image
img_1.load_image_2d(filename_1, current_log, data_folder)
img_2.load_image_2d(filename_2, current_log, data_folder)

# adjusts image levels
image1_uncorrected = img_1.data_2d / img_1.data_2d.max()
image2_uncorrected = img_2.data_2d / img_2.data_2d.max()
lower_threshold1 = 2 * filters.threshold_otsu(image1_uncorrected)
lower_threshold2 = 2 * filters.threshold_otsu(image2_uncorrected)

image1, hist1_before, hist1_after = image_adjust(
    image1_uncorrected, output_filename + "_ref", lower_threshold1, higher_threshold=0.9999, display=verbose,
)
image2, hist2_before, hist2_after = image_adjust(
    image2_uncorrected, output_filename, lower_threshold2, higher_threshold=0.9999, display=verbose,
)

# shows eq histograms
displaysEquializationHistograms(
    hist1_before, hist1_after, hist2_before, hist2_after, (lower_threshold1, lower_threshold2), verbose, output_filename,
)

# calculates shift
shift, error, diffphase = register_translation(image1, image2, 100)

# corrects image
# shift = (-22.4, 5.2)
# The shift corresponds to the pixel offset relative to the reference image
# offset_image = fourier_shift(np.fft.fftn(image2_uncorrected), shift)
# image2_corrected = np.fft.ifftn(offset_image)
image2_corrected = shift_image(image2, shift)
image2_corrected = exposure.rescale_intensity(image2_corrected, out_range=(0, 1))

print(f"Detected subpixel offset (y, x): {shift}")


if verbose:
    show_cc_image(image1, image2, verbose, output_filename)

save_2_images_rgb(image1_uncorrected, image2_uncorrected, output_filename + "_overlay_uncorrected.png")
image1 = image1 > 0.1
image2_corrected = image2_corrected > 0.1
save_2_images_rgb(image1, image2_corrected, output_filename + "_overlay_corrected.png")
