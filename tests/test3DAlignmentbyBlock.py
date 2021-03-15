#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 16:21:32 2021

@author: marcnol

this test will:
    - load test 3D fiducial files
    - reinterpolate second file in XY using dictionnary to get rough alignment
    - break in blocks
    - loop thru blocks and do:
        - cross correlate in 3D to find 3D shift
        - verify xy cross-correlation is close to zero
    - clean up final matrix considering outliers.
    - Consider smoothing by gaussian filtering?
    - Validate results plotting image with aligned masks
"""
#%% imports

import numpy as np
import matplotlib.pylab as plt

from skimage import io
import os, sys
import numpy as np
import matplotlib.pylab as plt
from imageProcessing.imageProcessing import (
    _reinterpolatesFocalPlane,
    imageShowWithValues,
    imageShowWithValuesSingle,
    imageAdjust,
    _removesInhomogeneousBackground,
    appliesXYshift3Dimages,
    imageBlockAlignment3D,
    plots3DshiftMatrices,
    combinesBlocksImageByReprojection,
    plots4images,
    makesShiftMatrixHiRes,
)

from scipy.stats import sigmaclip
from skimage.util.shape import view_as_blocks

from skimage.registration import phase_cross_correlation
from scipy.ndimage import shift as shiftImage
from tqdm import trange, tqdm
from skimage import exposure


#%%    - load test 3D fiducial files

if "atlantis" in os.uname()[1]:
    rootFolder = "/home/marcnol/data/Embryo_debug_dataset/Experiment_18"
else:
    rootFolder = "/home/marcnol/grey/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0"
files=["scan_001_RT27_001_ROI_converted_decon_ch00.tif","scan_001_RT29_001_ROI_converted_decon_ch00.tif"]
# files = ["scan_001_RT27_001_ROI_converted_decon_ch00.tif", "scan_001_RT41_001_ROI_converted_decon_ch00.tif"]

filenames = [rootFolder + os.sep + x for x in files]

print("\nReading files: \n{}".format(" \n\n ".join(filenames)))
images0 = [io.imread(x).squeeze() for x in filenames]

#%% settings

blockSizeXY = 128
upsample_factor=100
lower_threshold = 0.9
higher_threshold=0.9999
axes2Plot = range(3)

#%% makes

# images0= [x/x.max() for x in images0]
images0 = [exposure.rescale_intensity(x, out_range=(0, 1)) for x in images0]

print("Removing inhomogeneous background...")
images = [_removesInhomogeneousBackground(x) for x in images0]

print("Rescaling grey levels...")
images = [imageAdjust(x, lower_threshold=lower_threshold, higher_threshold=higher_threshold)[0] for x in images]

#%% calculates XY shift  and applies it
images_2D = [np.sum(x, axis=0) for x in images]

print("Calculating shifts...")
upsample_factor = 100
shift, error, diffphase = phase_cross_correlation(images_2D[0], images_2D[1], upsample_factor=upsample_factor)
print("shifts XY = {}".format(shift))

images_2D.append(shiftImage(images_2D[1], shift))

#%% shows original images and background substracted

allimages = images0 + images
plots4images(allimages, titles=['reference','cycle <i>','processed reference','processed cycle <i>'])

#%% reinterpolate second file in XY using dictionnary to get rough alignment

# this will only be done in the real pyHiM script. For now this has been simulated in the block above

images.append(appliesXYshift3Dimages(images[1], shift))

#%% shows images so far

# XY
fig, axes = plt.subplots(1, len(images))
fig.set_size_inches((len(images) * 5, 5))
ax = axes.ravel()

titles = ["reference", "uncorrected", "corrected"]
fig.suptitle("XZ {} images".format("-".join(titles)))

subVolume = ((0, 60), (0, 2048), (0, 2048))
subVolume = ((0, 60), (1000, 1250), (800, 1000))

cmap = "Greys"

for axis, img, title in zip(ax, images, titles):
    axis.imshow(
        np.sum(
            img[
                subVolume[0][0] : subVolume[0][1], subVolume[1][0] : subVolume[1][1], subVolume[2][0] : subVolume[2][1]
            ],
            axis=0,
        ),
        cmap=cmap,
        vmax=1,
    )

    axis.set_title(title)
    axis.set_xlabel("x")
    axis.set_ylabel("y")

# XZ
fig2, axes = plt.subplots(1, len(images))
fig2.set_size_inches((len(images) * 2, 5))
ax2 = axes.ravel()
fig2.suptitle("XZ {} images".format("-".join(titles)))

subVolume = ((0, 60), (1000, 1250), (1200, 1000))

for axis, img, title in zip(ax2, images, titles):
    axis.imshow(
        img[subVolume[0][0] : subVolume[0][1], subVolume[1][0] : subVolume[1][1], subVolume[2][0]].transpose(),
        cmap=cmap,
    )

    axis.set_title(title)
    axis.set_xlabel("z")
    axis.set_ylabel("x")

#%% 3D image alignment by block

shiftMatrices, block_ref, block_target = imageBlockAlignment3D(images, blockSizeXY=blockSizeXY, upsample_factor=upsample_factor)

#%% checks by plotting results

plots3DshiftMatrices(shiftMatrices, fontsize=8)


#%% - verify xy cross-correlation is close to zero


#%% - clean up final matrix considering outliers.


#%% - Consider smoothing by gaussian filtering?


#%% combines blocks into a single matrix for display instead of plotting a matrix of subplots each with a block

outputs = []
for axis in axes2Plot:
    outputs.append(combinesBlocksImageByReprojection(block_ref, block_target, shiftMatrices=shiftMatrices, axis1=axis))

SSIM_matrices = [x[1] for x in outputs]

fig = plt.figure(constrained_layout=False)
fig.set_size_inches((20 * 2, 20))
gs = fig.add_gridspec(2, 2)
ax = [fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])]

titles = ["Z-projection", "X-projection", "Y-projection"]

for axis, output, i in zip(ax, outputs, axes2Plot):
    axis.imshow(output[0])
    axis.set_title(titles[i])

fig.tight_layout()

newFig = plots3DshiftMatrices(SSIM_matrices, fontsize=6, log=False,valfmt="{x:.2f}")
newFig.suptitle("SSIM block matrices")

#%% makes HR shift matrix to save
# numberBlocks = block_ref.shape[0]
# blockSizeXY = block_ref.shape[3]

# shiftMatrix=np.zeros((3,blockSizeXY*shiftMatrices[0].shape[0],blockSizeXY*shiftMatrices[0].shape[1]))
# for _ax,m in enumerate(shiftMatrices):
#     print("size={}".format(m.shape))
#     for i in range(numberBlocks):
#         for j in range(numberBlocks):
#             shiftMatrix[_ax,i * blockSizeXY: (i + 1) * blockSizeXY,j * blockSizeXY: (j + 1) * blockSizeXY] = m[i,j]

# aaa=[shiftMatrix[0,:,:],shiftMatrix[1,:,:],shiftMatrix[2,:,:]]
shiftMatrix = makesShiftMatrixHiRes(shiftMatrices, block_ref.shape)

plots3DshiftMatrices(shiftMatrix, fontsize=8)

#%% - Validate results plotting image with aligned masks

def plotBlocks(block_ref, block_target, shiftMatrices, cmap="RdBu", axis1=0):
    a = 1
    Nx, Ny = block_ref.shape[0:2]
    fig, axes = plt.subplots(Nx, Ny)
    fig.set_size_inches((Nx * a, Ny * a))
    ax = axes.ravel()

    iplot = 0
    for i in trange(block_ref.shape[0]):
        for j in range(block_ref.shape[1]):
            axis = ax[iplot]

            # aligns block

            imgs = list()
            imgs.append(block_ref[i, j])

            shift3D = np.array([x[i, j] for x in shiftMatrices])
            imgs.append(shiftImage(block_target[i, j], shift3D))

            imgs = [x / x.max() for x in imgs]

            imgs = [np.sum(x, axis=axis1) for x in imgs]
            # imgs=[_removesInhomogeneousBackground(x) for x in imgs]
            imgs = [imageAdjust(x, lower_threshold=0.5, higher_threshold=0.9999)[0] for x in imgs]

            if axis1 > 0:
                imgs = [x.transpose() for x in imgs]

            imgs.append(np.zeros(imgs[0].shape))
            RGB = np.dstack(imgs)
            axis.imshow(RGB)

            axis.axes.xaxis.set_visible(False)
            axis.axes.yaxis.set_visible(False)

            iplot += 1

    # plt.subplots_adjust(wspace=0.05,hspace=0.05)
    fig.tight_layout()


for axis in range(2):
    plotBlocks(block_ref, block_target, shiftMatrices,cmap="RdBu", axis1=axis)

# fig.tight_layout()

# write this first

# - make list of axis
# - make figure with NxN subplots,
# - plot each block in a different axis
# - remove axis labels
# - repeat for xy and for xz



# plt.imshow(shiftMatrix[0,:,:])
