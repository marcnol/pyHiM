#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 16:00:31 2020

@author: marcnol

This series of scripts was written to test and prototype:
    
    - alignment of masks using CV2
    - [Not ported] testing translations and afine transformations
    - [Not ported] testing multi-scale translational alignment by zooming out of the center of the image
    - [ported] Displaying image difference to assess alignment
    - [ported] Test block translational alignment

"""
import numpy as np
import cv2
import os
import matplotlib.pyplot as plt
from datetime import datetime
from skimage.util.shape import view_as_blocks
from numpy import linalg as LA
from tqdm import trange
from skimage import data, segmentation, measure
from skimage import filters
from skimage.feature import register_translation
from skimage.registration import phase_cross_correlation
from scipy.ndimage import shift as shift_image
import cv2
from skimage.exposure import match_histograms

## Functions


def align_cv2(im1, im2, warp_mode):

    # Find size of image1
    sz = im1.shape

    # Define 2x3 or 3x3 matrices and initialize the matrix to identity
    if warp_mode == cv2.MOTION_HOMOGRAPHY:
        warp_matrix = np.eye(3, 3, dtype=np.float32)
    else:
        warp_matrix = np.eye(2, 3, dtype=np.float32)

    # Specify the number of iterations.
    number_of_iterations = 1000

    # Specify the threshold of the increment
    # in the correlation coefficient between two iterations
    termination_eps = 1e-10

    # Define termination criteria
    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, number_of_iterations, termination_eps)

    # # Run the ECC algorithm. The results are stored in warp_matrix.
    # try:
    #     cc, warp_matrix = cv2.findTransformECC(im1,im2,warp_matrix, warp_mode, criteria, inputMask=None, gaussFiltSize=1)
    # except TypeError:
    #     cc, warp_matrix = cv2.findTransformECC(im1,im2,warp_matrix, warp_mode, criteria)

    # Run the ECC algorithm. The results are stored in warp_matrix.
    try:
        cc, warp_matrix = cv2.findTransformECC(
            im1, im2, warp_matrix, warp_mode, criteria, inputMask=None, gaussFiltSize=1
        )
    except TypeError:
        cc, warp_matrix = cv2.findTransformECC(im1, im2, warp_matrix, warp_mode, criteria)
    except cv2.error:
        cc = 0
        print("Warning: find transform failed. Set warp as identity")

    return cc, warp_matrix


def apply_correction(im2, warp_matrix):

    sz = im2.shape

    # Use warpAffine for Translation, Euclidean and Affine
    im2_aligned = cv2.warpAffine(im2, warp_matrix, (sz[1], sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP)

    return im2_aligned


def showImages(image1, image2, im2_aligned, file_name):

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches((60, 30))
    cmap = "seismic"
    ax1.imshow(image1 - image2, cmap=cmap)
    ax1.axis("off")

    ax2.imshow(image1 - im2_aligned, cmap=cmap)
    ax2.axis("off")

    plt.savefig(file_name)


def align_images_by_blocks(img_1, img_2, block_size, upsample_factor=100):

    block_1 = view_as_blocks(image1, block_size)
    block_2 = view_as_blocks(image2, block_size)
    warp_matrix = np.eye(2, 3, dtype=np.float32)

    shift_image_norm = np.zeros((block_1.shape[0], block_1.shape[1]))
    shifted_image = np.zeros((block_1.shape[0], block_1.shape[1], 2))
    rms_image = np.zeros((block_1.shape[0], block_1.shape[1]))

    for i in trange(block_1.shape[0]):
        for j in range(block_1.shape[1]):

            # shift, error, diffphase = register_translation(block_1[i,j], block_2[i,j],upsample_factor=upsample_factor)
            shift, error, diffphase = phase_cross_correlation(
                block_1[i, j], block_2[i, j], upsample_factor=upsample_factor
            )

            shift_image_norm[i, j] = LA.norm(shift)
            shifted_image[i, j, 0], shifted_image[i, j, 1] = shift[0], shift[1]
            img_2_aligned = shift_image(img_2, shift)

            # cc, warp_matrix = align_cv2(block_1[i,j], block_2[i,j], warp_mode)
            # shift_image_norm[i,j] = LA.norm(warp_matrix[:,2])
            # shifted_image[i,j,0],shifted_image[i,j,1] = warp_matrix[:,2][0], warp_matrix[:,2][1]
            # img_2_aligned = apply_correction(img_2,warp_matrix)

            rms_image[i, j] = np.sum(np.sum(np.abs(img_1 - img_2_aligned), axis=1))

    # calculates optimal shifts by polling blocks showing the best RMS
    tolerance = 0.1

    # threshold = filters.threshold_otsu(rms_image)
    threshold = (1 + tolerance) * np.min(rms_image)
    mask = rms_image < threshold
    contours = measure.find_contours(rms_image, threshold)
    try:
        contour = sorted(contours, key=lambda x: len(x))[-1]
    except IndexError:
        contour = np.array([0, 0])

    mean_shifts = [np.mean(shifted_image[mask, 0]), np.mean(shifted_image[mask, 1])]
    warp_matrix[:, 2] = mean_shifts
    std_shifts = [np.std(shifted_image[mask, 0]), np.std(shifted_image[mask, 1])]
    mean_shift_norm = np.mean(shift_image_norm[mask])

    mean_error = np.mean(rms_image[mask])

    relative_shifts = np.abs(shift_image_norm - mean_shift_norm)

    # if it does not have enough pollsters to fall back to then it does a global cross correlation!
    if np.sum(mask) < 4:
        cc, warp_matrix_global = align_cv2(img_1, img_2, warp_mode)
        mean_shifts_global = warp_matrix_global[:, 2]
        img_2_aligned_global = apply_correction(img_2, warp_matrix_global)
        mean_error_global = np.sum(np.sum(np.abs(img_1 - img_2_aligned_global), axis=1))

        if mean_error_global < mean_error:
            mean_shifts = mean_shifts_global
            warp_matrix = warp_matrix_global
            mean_error = mean_error_global
            print("Falling back to global registration")

    print(
        "Mean polled XY shifts: {:.2f}({:.2f}) px | {:.2f}({:.2f}) px".format(
            mean_shifts[0], std_shifts[0], mean_shifts[1], std_shifts[1]
        )
    )

    # We apply the consensous correction and we calculate how well it does for each block to get a benchmark
    rmsBlock = np.zeros((block_1.shape[0], block_1.shape[1]))  # new
    for i in trange(block_1.shape[0]):
        for j in range(block_1.shape[1]):
            Block2aligned = shift_image(block_2[i, j], mean_shifts)  # new
            rmsBlock[i, j] = np.sum(np.sum(np.abs(block_1[i, j] - Block2aligned), axis=1)) / (
                np.sum(block_1[i, j]) + np.sum(Block2aligned)
            )  # new

    return mean_shifts, warp_matrix, mean_error, mask, relative_shifts, rms_image, rmsBlock, contour


def plotting_block_alignment_results(relative_shifts, rms_image, rmsBlock, contour, file_name="BlockALignmentResults.png"):

    # plotting
    fig, axes = plt.subplots(2, 2)
    ax = axes.ravel()
    fig.set_size_inches((10, 10))

    cbwindow = 3
    p_1 = ax[0].imshow(relative_shifts, cmap="terrain", vmin=0, vmax=cbwindow)
    ax[0].plot(contour.T[1], contour.T[0], linewidth=2, c="k")
    ax[0].set_title("abs(global-block) shifts, px")
    fig.colorbar(p_1, ax=ax[0], fraction=0.046, pad=0.04)

    p_2 = ax[1].imshow(rms_image, cmap="terrain", vmin=np.min(rms_image), vmax=np.max(rms_image))
    ax[1].plot(contour.T[1], contour.T[0], linewidth=2, c="k")
    ax[1].set_title("RMS-image")
    fig.colorbar(p_2, ax=ax[1], fraction=0.046, pad=0.04)

    # add
    p3 = ax[2].imshow(rmsBlock, cmap="terrain", vmin=np.min(rmsBlock), vmax=np.max(rmsBlock))
    ax[2].set_title("RMS-block")
    fig.colorbar(p3, ax=ax[2], fraction=0.046, pad=0.04)

    for x in range(len(ax)):
        ax[x].axis("off")


#%% loading and preprocessing of example images

## main
root_folder = "/home/marcnol/data/Embryo_debug_dataset/Experiment_18/zProject"
# imFileName1=root_folder+os.sep+"scan_001_RT27_001_ROI_converted_decon_ch00_2d.npy"
# imFileName2=root_folder+os.sep+"scan_001_RT37_001_ROI_converted_decon_ch00_2d.npy"

imFileName1 = root_folder + os.sep + "scan_001_RT27_003_ROI_converted_decon_ch00_2d.npy"
imFileName2 = root_folder + os.sep + "scan_001_RT40_003_ROI_converted_decon_ch00_2d.npy"

image1 = np.load(imFileName1).squeeze()
image2 = np.load(imFileName2).squeeze()

images = (image1, image2)

image1, image2 = list(map(lambda x: np.float32(x / x.max()), images))

# image1 = image1 / image1.max()
# image2 = image2 / image2.max()

# image1 = np.float32(image1)
# image2 = np.float32(image2)

#%% Testing translational and afine transformations with CV2

# Define the motion model
# model='afine'
model = "translate"
# model='euclidean'

if model == "translate":
    warp_mode = cv2.MOTION_TRANSLATION
elif model == "afine":
    warp_mode = cv2.MOTION_AFFINE
elif model == "euclidean":
    warp_mode = cv2.MOTION_EUCLIDEAN

output_filename = "/home/marcnol/Downloads/" + "test_" + model + ".png"

begin_time = datetime.now()

cc, warp_matrix = align_cv2(image1, image2, warp_mode)
im2_aligned = apply_correction(image2, warp_matrix)
print("Elapsed time: {}".format(datetime.now() - begin_time))

begin_time = datetime.now()

# shift, error, diffphase = register_translation(image1,image2,upsample_factor=100)
shift, error, diffphase = phase_cross_correlation(image1, image2, upsample_factor=100)
im2_aligned = shift_image(image2, shift)

print("Elapsed time: {}".format(datetime.now() - begin_time))


print("Done registering with CV: Shifts={}".format(warp_matrix[:, 2]))
print("Done registering with Skkimage: Shifts={}".format(shift))

# showing rgb image
sz = image1.shape

fig, (ax1, ax2) = plt.subplots(1, 2)
fig.set_size_inches((60, 30))

null_image = np.zeros(sz)

rgb = np.dstack([image1, image2, null_image])
ax1.imshow(rgb)
ax1.axis("off")

RGB_aligned = np.dstack([image1, im2_aligned, null_image])
ax2.imshow(RGB_aligned)
ax2.axis("off")

plt.savefig(output_filename)


#%% Testing difference image

showImages(image1, image2, im2_aligned, "/home/marcnol/Downloads/" + "difference_" + model + ".png")

#%% test multiscale alignments

image_size = image1.shape[0]
c = int(np.round(image_size / 2))

factors = [16, 8, 4, 2, 1]
windows = [int(image_size / x) for x in factors]
shifts = []
rms = []
CC = []

for window in windows:
    if window == image_size:
        img_1 = image1
        img_2 = image2
    else:
        img_1 = image1[c - window : c + window, c - window : c + window]
        img_2 = image2[c - window : c + window, c - window : c + window]

    cc, warp_matrix = align_cv2(img_1, img_2, warp_mode)

    img_2_aligned = apply_correction(img_2, warp_matrix)

    shifts.append(warp_matrix[:, 2])
    rms.append(np.sum(np.sum(np.abs(img_1 - img_2_aligned), axis=1)))
    CC.append(cc)

    print("Aligned images for window: {}".format(window))
    file_name = "/home/marcnol/Downloads/" + "difference_" + str(window) + ".png"

    showImages(img_1, img_2, img_2_aligned, file_name)


# display

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.set_size_inches((15, 5))

ax1.plot(windows, CC, "o-")
ax1.set_ylabel("cc")

ax2.plot(windows, rms, "o-")
ax2.set_ylabel("rms")

ax3.plot(windows, [x[0] for x in shifts] - shifts[-1][0], "o-")
ax3.plot(windows, [x[1] for x in shifts] - shifts[-1][1], "+-")
ax3.set_ylabel("relative shifts")


#%% test blocks alignments

# image2=match_histograms(image2,image1)
upsample_factor = 100
block_size = (256, 256)

mean_shifts, warp_matrix, mean_error, mask, relative_shifts, rms_image, rmsBlock, contour = align_images_by_blocks(
    image1, image2, block_size, upsample_factor
)

plotting_block_alignment_results(relative_shifts, rms_image, rmsBlock, contour)
