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
from scipy.ndimage import shift as shiftImage
import cv2
from skimage.exposure import match_histograms

## Functions


def alignCV2(im1, im2, warp_mode):

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


def applyCorrection(im2, warp_matrix):

    sz = im2.shape

    # Use warpAffine for Translation, Euclidean and Affine
    im2_aligned = cv2.warpAffine(im2, warp_matrix, (sz[1], sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP)

    return im2_aligned


def showImages(image1, image2, im2_aligned, fileName):

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches((60, 30))
    cmap = "seismic"
    ax1.imshow(image1 - image2, cmap=cmap)
    ax1.axis("off")

    ax2.imshow(image1 - im2_aligned, cmap=cmap)
    ax2.axis("off")

    plt.savefig(fileName)


def alignImagesByBlocks(I1, I2, blockSize, upsample_factor=100):

    Block1 = view_as_blocks(image1, blockSize)
    Block2 = view_as_blocks(image2, blockSize)
    warp_matrix = np.eye(2, 3, dtype=np.float32)

    shiftImageNorm = np.zeros((Block1.shape[0], Block1.shape[1]))
    shiftedImage = np.zeros((Block1.shape[0], Block1.shape[1], 2))
    rmsImage = np.zeros((Block1.shape[0], Block1.shape[1]))

    for i in trange(Block1.shape[0]):
        for j in range(Block1.shape[1]):

            # shift, error, diffphase = register_translation(Block1[i,j], Block2[i,j],upsample_factor=upsample_factor)
            shift, error, diffphase = phase_cross_correlation(
                Block1[i, j], Block2[i, j], upsample_factor=upsample_factor
            )

            shiftImageNorm[i, j] = LA.norm(shift)
            shiftedImage[i, j, 0], shiftedImage[i, j, 1] = shift[0], shift[1]
            I2_aligned = shiftImage(I2, shift)

            # cc, warp_matrix = alignCV2(Block1[i,j], Block2[i,j], warp_mode)
            # shiftImageNorm[i,j] = LA.norm(warp_matrix[:,2])
            # shiftedImage[i,j,0],shiftedImage[i,j,1] = warp_matrix[:,2][0], warp_matrix[:,2][1]
            # I2_aligned = applyCorrection(I2,warp_matrix)

            rmsImage[i, j] = np.sum(np.sum(np.abs(I1 - I2_aligned), axis=1))

    # calculates optimal shifts by polling blocks showing the best RMS
    tolerance = 0.1

    # threshold = filters.threshold_otsu(rmsImage)
    threshold = (1 + tolerance) * np.min(rmsImage)
    mask = rmsImage < threshold
    contours = measure.find_contours(rmsImage, threshold)
    try:
        contour = sorted(contours, key=lambda x: len(x))[-1]
    except IndexError:
        contour = np.array([0, 0])

    meanShifts = [np.mean(shiftedImage[mask, 0]), np.mean(shiftedImage[mask, 1])]
    warp_matrix[:, 2] = meanShifts
    stdShifts = [np.std(shiftedImage[mask, 0]), np.std(shiftedImage[mask, 1])]
    meanShiftNorm = np.mean(shiftImageNorm[mask])

    meanError = np.mean(rmsImage[mask])

    relativeShifts = np.abs(shiftImageNorm - meanShiftNorm)

    # if it does not have enough pollsters to fall back to then it does a global cross correlation!
    if np.sum(mask) < 4:
        cc, warp_matrix_global = alignCV2(I1, I2, warp_mode)
        meanShifts_global = warp_matrix_global[:, 2]
        I2_aligned_global = applyCorrection(I2, warp_matrix_global)
        meanError_global = np.sum(np.sum(np.abs(I1 - I2_aligned_global), axis=1))

        if meanError_global < meanError:
            meanShifts = meanShifts_global
            warp_matrix = warp_matrix_global
            meanError = meanError_global
            print("Falling back to global registration")

    print(
        "Mean polled XY shifts: {:.2f}({:.2f}) px | {:.2f}({:.2f}) px".format(
            meanShifts[0], stdShifts[0], meanShifts[1], stdShifts[1]
        )
    )

    # We apply the consensous correction and we calculate how well it does for each block to get a benchmark
    rmsBlock = np.zeros((Block1.shape[0], Block1.shape[1]))  # new
    for i in trange(Block1.shape[0]):
        for j in range(Block1.shape[1]):
            Block2aligned = shiftImage(Block2[i, j], meanShifts)  # new
            rmsBlock[i, j] = np.sum(np.sum(np.abs(Block1[i, j] - Block2aligned), axis=1)) / (
                np.sum(Block1[i, j]) + np.sum(Block2aligned)
            )  # new

    return meanShifts, warp_matrix, meanError, mask, relativeShifts, rmsImage, rmsBlock, contour


def plottingBlockALignmentResults(relativeShifts, rmsImage, rmsBlock, contour, fileName="BlockALignmentResults.png"):

    # plotting
    fig, axes = plt.subplots(2, 2)
    ax = axes.ravel()
    fig.set_size_inches((10, 10))

    cbwindow = 3
    p1 = ax[0].imshow(relativeShifts, cmap="terrain", vmin=0, vmax=cbwindow)
    ax[0].plot(contour.T[1], contour.T[0], linewidth=2, c="k")
    ax[0].set_title("abs(global-block) shifts, px")
    fig.colorbar(p1, ax=ax[0], fraction=0.046, pad=0.04)

    p2 = ax[1].imshow(rmsImage, cmap="terrain", vmin=np.min(rmsImage), vmax=np.max(rmsImage))
    ax[1].plot(contour.T[1], contour.T[0], linewidth=2, c="k")
    ax[1].set_title("RMS-image")
    fig.colorbar(p2, ax=ax[1], fraction=0.046, pad=0.04)

    # add
    p3 = ax[2].imshow(rmsBlock, cmap="terrain", vmin=np.min(rmsBlock), vmax=np.max(rmsBlock))
    ax[2].set_title("RMS-block")
    fig.colorbar(p3, ax=ax[2], fraction=0.046, pad=0.04)

    for x in range(len(ax)):
        ax[x].axis("off")


#%% loading and preprocessing of example images

## main
rootFolder = "/home/marcnol/data/Embryo_debug_dataset/Experiment_18/zProject"
# imFileName1=rootFolder+os.sep+"scan_001_RT27_001_ROI_converted_decon_ch00_2d.npy"
# imFileName2=rootFolder+os.sep+"scan_001_RT37_001_ROI_converted_decon_ch00_2d.npy"

imFileName1 = rootFolder + os.sep + "scan_001_RT27_003_ROI_converted_decon_ch00_2d.npy"
imFileName2 = rootFolder + os.sep + "scan_001_RT40_003_ROI_converted_decon_ch00_2d.npy"

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

outputFileName = "/home/marcnol/Downloads/" + "test_" + model + ".png"

begin_time = datetime.now()

cc, warp_matrix = alignCV2(image1, image2, warp_mode)
im2_aligned = applyCorrection(image2, warp_matrix)
print("Elapsed time: {}".format(datetime.now() - begin_time))

begin_time = datetime.now()

# shift, error, diffphase = register_translation(image1,image2,upsample_factor=100)
shift, error, diffphase = phase_cross_correlation(image1, image2, upsample_factor=100)
im2_aligned = shiftImage(image2, shift)

print("Elapsed time: {}".format(datetime.now() - begin_time))


print("Done registering with CV: Shifts={}".format(warp_matrix[:, 2]))
print("Done registering with Skkimage: Shifts={}".format(shift))

# showing RGB image
sz = image1.shape

fig, (ax1, ax2) = plt.subplots(1, 2)
fig.set_size_inches((60, 30))

nullImage = np.zeros(sz)

RGB = np.dstack([image1, image2, nullImage])
ax1.imshow(RGB)
ax1.axis("off")

RGB_aligned = np.dstack([image1, im2_aligned, nullImage])
ax2.imshow(RGB_aligned)
ax2.axis("off")

plt.savefig(outputFileName)


#%% Testing difference image

showImages(image1, image2, im2_aligned, "/home/marcnol/Downloads/" + "difference_" + model + ".png")

#%% test multiscale alignments

imageSize = image1.shape[0]
c = int(np.round(imageSize / 2))

factors = [16, 8, 4, 2, 1]
windows = [int(imageSize / x) for x in factors]
shifts = list()
rms = list()
CC = list()

for window in windows:
    if window == imageSize:
        I1 = image1
        I2 = image2
    else:
        I1 = image1[c - window : c + window, c - window : c + window]
        I2 = image2[c - window : c + window, c - window : c + window]

    cc, warp_matrix = alignCV2(I1, I2, warp_mode)

    I2_aligned = applyCorrection(I2, warp_matrix)

    shifts.append(warp_matrix[:, 2])
    rms.append(np.sum(np.sum(np.abs(I1 - I2_aligned), axis=1)))
    CC.append(cc)

    print("Aligned images for window: {}".format(window))
    fileName = "/home/marcnol/Downloads/" + "difference_" + str(window) + ".png"

    showImages(I1, I2, I2_aligned, fileName)


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
blockSize = (256, 256)

meanShifts, warp_matrix, meanError, mask, relativeShifts, rmsImage, rmsBlock, contour = alignImagesByBlocks(
    image1, image2, blockSize, upsample_factor
)

plottingBlockALignmentResults(relativeShifts, rmsImage, rmsBlock, contour)
