#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 16:58:03 2020

@author: marcnol
"""


# project a number of images and makes an RGB image


import os
import glob
import matplotlib.pylab as plt
import numpy as np
from skimage import io
import cv2
import matplotlib.gridspec as gridspec

from imageProcessing.imageProcessing import Image, imageAdjust
import imageProcessing.makeProjections
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground
from astropy.visualization import SqrtStretch, simple_norm

#%% loads all images
rootFolder = "/home/marcnol/grey/Olivier/2020_08_13_Experiment_tetraspeck/Deconvolved/Tetraspeck_Scan_1"
rootFolder = (
    "/mnt/grey/DATA/Olivier/2020_08_14_Experiment_Chromatic_aberration_Embryo_Locus2L/Deconvolved/Embryo_Scan_11"
)
files = glob.glob(rootFolder + os.sep + "*.tif")

print("Files to load {}\n".format("\n".join([os.path.basename(x) for x in files])))

images = list()
boxSize = (64, 64)
# boxSize=(16, 16)
for file in files:

    # loads image
    print("Loading {}...".format(os.path.basename(file)))
    im = io.imread(file).squeeze()

    imageSize = im.shape

    # makes actual 2d projection
    (zmin, zmax) = (0, imageSize[0])
    zRange = (round((zmin + zmax) / 2), range(zmin, zmax))

    print("Summing {}...".format(os.path.basename(file)))

    im2D = np.max(im[zRange[1][0] : zRange[1][-1]], axis=0)

    # im2D = im2D - im2D[0,0]
    # im2D = im2D / im2D.max()
    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()

    bkg = Background2D(im2D, boxSize, filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,)

    im2D_substracted = im2D - bkg.background
    im2D_substracted = im2D_substracted - np.min(im2D_substracted)
    images.append(im2D)


#%% displays
outputFileNameRoot = "/home/marcnol/Downloads" + os.sep + "RGB"
Nplots = len(files)
fig = plt.figure(constrained_layout=False, figsize=(10 * Nplots, 10), dpi=300, facecolor="w", edgecolor="k")
nCols, nRows = Nplots, 1
spec2 = gridspec.GridSpec(ncols=nCols, nrows=nRows, figure=fig)

FigList, Yticks, Xticks = [], [], []
for iRow in range(nRows):
    for iCol in range(nCols):
        FigList.append(fig.add_subplot(spec2[iRow, iCol]))

outputFileName = outputFileNameRoot + "single_ch" + ".png"

percents = [97.0, 99.9, 99.9]
# percents = [99.99, 99., 99.9, 99.]

for ax, im, file, percent in zip(FigList, images, files, percents):
    title = os.path.basename(file).split(".")[0].split("_")[-1]
    norm = simple_norm(im, "sqrt", percent=percent)
    ax.imshow(im, cmap="Greys", origin="lower", norm=simple_norm(im, "sqrt", percent=percent))
    ax.set_title(title)

fig.savefig(outputFileName)

# plt.close(fig)
#%% plots RGB images

thresholds = [[0.9, 0.99], [0.9, 0.999], [0.85, 0.999], [0.85, 0.90]]

cmaps = ["Blues", "Greens", "Reds"]
for i in range(len(images) - 2):
    fig, ax = plt.subplots()
    outputFileName = outputFileNameRoot + str(i) + ":" + str(i + 2) + ".png"

    fig.set_size_inches((30, 30))

    images2 = []
    for j, x in enumerate(images[i : i + 3]):
        x = x - x.min()
        x = x / x.max()
        x, _, _, _, _ = imageAdjust(x, lower_threshold=thresholds[j][0], higher_threshold=thresholds[j][1])
        images2.append(x)

    RGB = np.dstack(images[i : i + 3])
    # for I,cmap in zip(images2,cmaps):
    #     ax.imshow(I, cmap=cmap, origin="lower", alpha = 0.3, norm=simple_norm(I, "sqrt", percent=percent))

    RGB = np.dstack(images2)
    ax.imshow(RGB)
    ax.axis("off")
    fig.savefig(outputFileName)

    # plt.close(fig)


#%%

# Olivier Messina <Olivier469110@hotmail.fr>
#
# Aug 14, 2020, 5:19 PM
#
# to me
# Experiment TetraSpeck 0,1 µm

# 1. TetraSpeck diluted 1/500 15uL
# → Channel used : 405 / 488 / 561 / 647
# → Acquisition parameters : z range 15um // step 250nm.
# → Laser power : 405(60%) / 488(60%) / 561(60%) / 647(60%)

# /mnt/grey/DATA/Olivier/2020_08_13_Experiment_tetraspeck/002_Scan_Tetraspeck_scan_2_15um

# /mnt/grey/DATA/Olivier/2020_08_13_Experiment_tetraspeck/003_Scan_Tetraspeck_scan_3_15um

# 2. TetraSpeck diluted 1/500 40uL
# → Channel used : 405 / 488 / 561 / 647
# → Acquisition parameters : z range 15um // step 250nm.
# → Laser power : 405(60%) / 488(60%) / 561(60%) / 647(60%)

# /mnt/grey/DATA/Olivier/2020_08_13_Experiment_tetraspeck/004_Scan_Tetraspeck_scan_4_15um

# /mnt/grey/DATA/Olivier/2020_08_13_Experiment_tetraspeck/005_Scan_Tetraspeck_scan_5_15um

# /mnt/grey/DATA/Olivier/2020_08_13_Experiment_tetraspeck/006_Scan_Tetraspeck_scan_6_15um

# /mnt/grey/DATA/Olivier/2020_08_13_Experiment_tetraspeck/008_Scan_Tetraspeck_scan_8_15um


# Experiment Embryos Labelled new library locus 2L

# → Channel used : 488 (Primer) / 561(Bc3-4-5-6-7) / 647(Bc3-4-5-6-7) [Barcode with 1 binding site for the imaging oligo]
# → Acquisition parameters : z range 20um // step 250nm.
# → Laser power : 488(60%) / 561(60%) / 647(60%)

# /mnt/grey/DATA/Olivier/2020_08_14_Experiment_Chromatic_aberration_Embryo_Locus2L/010_Scan_Embryo_488_561_647_scan9_combined_RTs

# /mnt/grey/DATA/Olivier/2020_08_14_Experiment_Chromatic_aberration_Embryo_Locus2L/011_Scan_Embryo_488_561_647_scan10_combined_RTs

# /mnt/grey/DATA/Olivier/2020_08_14_Experiment_Chromatic_aberration_Embryo_Locus2L/012_Scan_Embryo_488_561_647_scan11_combined_RTs

# /mnt/grey/DATA/Olivier/2020_08_14_Experiment_Chromatic_aberration_Embryo_Locus2L/015_Scan_Embryo_488_561_647_scan14_combined_RTs

# /mnt/grey/DATA/Olivier/2020_08_14_Experiment_Chromatic_aberration_Embryo_Locus2L/016_Scan_Embryo_488_561_647_scan15_combined_RTs
