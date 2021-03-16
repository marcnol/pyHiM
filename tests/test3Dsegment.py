#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:05:02 2021

@author: marcnol

- tests several algorithms for 3D segmentation:
    - based on ASTROPY (plane by plane object-based segmentation and deblending)
    - blob_log from skimage
    - 3D centroids from 3D masks using:
        - regionprops (weighted centroids)
        - bigfish to get 3D gaussian fits

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


#%% This first test will load a pre-processed 3D image and will run _segments3DvolumesByThresholding to
# find the objects in 3D.

rootFolder="/mnt/grey/DATA/users/marcnol/models/StarDist3D/training3Dbarcodes/dataset1/"
file = rootFolder+'scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0.tif'

threshold_over_std,sigma ,boxSize, filter_size=1, 3, (32, 32),(3, 3)
nlevels=64
contrast=0.001
image3D = io.imread(file).squeeze()

binary, segmentedImage3D = _segments3DvolumesByThresholding(image3D,
                                                       threshold_over_std=threshold_over_std,
                                                       sigma = 3,
                                                       boxSize=(32, 32),
                                                       filter_size=(3, 3),
                                                       nlevels=nlevels,
                                                       contrast=contrast,
                                                       deblend3D=True)


display3D(image3D=image3D,labels=segmentedImage3D,z=40, rangeXY=1000)

# decompose output labeled image in blocks

fileName = rootFolder + "blockTest"
savesImageAsBlocks(segmentedImage3D,fileName,blockSizeXY=256)

# Deblend image in 3D by watersheding

binary0=binary>0

# Now we want to separate the two objects in image
# Generate the markers as local maxima of the distance to the background
print("Constructing distance matrix from 3D binary mask...")
distance = ndi.distance_transform_edt(binary0)
coords = peak_local_max(distance, footprint=np.ones((10, 10, 25)), labels=binary0)
mask = np.zeros(distance.shape, dtype=bool)
mask[tuple(coords.T)] = True
markers, _ = ndi.label(mask)
print("Deblending sources in 3D by watersheding...")
labels = watershed(-distance, markers, mask=binary0)

display3D(image3D=image3D,labels=labels,z=40, rangeXY=1000)

#%% loads and display output of a previous segmentation run

folder = "/home/marcnol/grey/users/marcnol/models/StarDist3D/training3Dbarcodes/dataset1/"
file_image3D=folder+os.sep+"scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0.tif"

file_segmented=folder+os.sep+"scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0_segmented.tif"
image3D = io.imread(file_image3D).squeeze()
labeled = io.imread(file_segmented).squeeze()

display3D(image3D=image3D,labels=labeled,z=40, rangeXY=1000)


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

#%% big-FISH

from bigfish.detection.spot_modeling import fit_subpixel
from skimage.measure import regionprops

preProcesses=True

z = 40
threshold_over_std,sigma ,boxSize, filter_size=1, 3, (32, 32),(3, 3)
nlevels=64
contrast=0.001
lower_threshold = 0.9
higher_threshold = 0.9999

# loads pre-processed image
if preProcesses:
    # pre-processes image
    if "atlantis" in os.uname()[1]:
        rootFolder="/home/marcnol/data/Embryo_debug_dataset/test_dataset/"
        file = rootFolder+'scan_001_RT31_001_ROI_converted_decon_ch01.tif'
    else:
        rootFolder="/mnt/grey/DATA/users/marcnol/models/StarDist3D/training3Dbarcodes/dataset1/"
        file = rootFolder+'scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0.tif'

    imageRaw = io.imread(file).squeeze()
    # autoscales exposures
    image3D = exposure.rescale_intensity(imageRaw, out_range=(0, 1))

    # removes inhomogeneous background
    print("\nRemoving inhomogeneous background...")
    image3D= _removesInhomogeneousBackground(image3D)

    # rescales grey levels
    print("\nRescaling grey levels...")
    image3D,_,_,_,_ = imageAdjust(image3D, lower_threshold=lower_threshold,higher_threshold=higher_threshold)
    shift = [-6.185384615384616, -2.9223076923076925]
    shift = np.array(shift)
    image3D = appliesXYshift3Dimages(image3D, shift)

else:

    if "atlantis" in os.uname()[1]:
        rootFolder="/home/marcnol/data/Embryo_debug_dataset/"
        file = rootFolder+'scan_001_RT27_001_ROI_converted_decon_ch01_preProcessed_index:0.tif'
    else:
        rootFolder="/mnt/grey/DATA/users/marcnol/models/StarDist3D/training3Dbarcodes/dataset1/"
        file = rootFolder+'scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0.tif'

    print("loading image {}".format(os.path.basename(file)))

    image3D = io.imread(file).squeeze()

image3D_test=image3D

# segments in 3D using ASTROPY
binary, segmentedImage3D = _segments3DvolumesByThresholding(image3D_test,
                                                       threshold_over_std=threshold_over_std,
                                                       sigma = 3,
                                                       boxSize=(32, 32),
                                                       filter_size=(3, 3),
                                                       nlevels=nlevels,
                                                       contrast=contrast,
                                                       deblend3D=True)

# finds first estimate of centroids from 3D masks
image3D_test = exposure.rescale_intensity(image3D_test, out_range=(0, 1))

properties = regionprops(segmentedImage3D, intensity_image=image3D_test, )
centroids=[x.weighted_centroid for x in properties]

z=[x[0] for x in centroids]
y=[x[1] for x in centroids]
x=[x[2] for x in centroids]

localizationTable = np.zeros((len(z),3))
localizationTable[:,0]=x
localizationTable[:,1]=y
localizationTable[:,2]=z

# fits 3D gaussians using BIGFISH
spots = np.zeros((len(z),3))
spots[:,0]=z
spots[:,1]=y
spots[:,2]=x
spots=spots.astype('int64')
print("Running bigfish...")
spots_subpixel = fit_subpixel(image3D_test, spots, voxel_size_z=250, voxel_size_yx=100,psf_z=500, psf_yx=200)

# Plots results
display3D(image3D=image3D_test,localizationsList = [spots,spots_subpixel],labels=segmentedImage3D,z=40, rangeXY=1000, norm=True,cmap='Greys')

# represents image in 3D with localizations
img = segmentedImage3D
img = image3D
center = int(img.shape[1]/2)
window = 10

images = list()
images.append(np.sum(img,axis=0))
images.append(np.sum(img[:,:,center-window:center+window],axis=2))
images.append(np.sum(img[:,center-window:center+window,:],axis=1))

# produces and saves output figure
figures=list()
figures.append([display3D_assembled(images, localizations = [spots_subpixel,spots], plottingRange = [center,window]),'_3DimageNlocalizations.png'])

# saves figures
outputFileNames = [rootFolder+os.path.basename(file)+x[1] for x in figures]

for fig, file in zip(figures,outputFileNames):
    fig[0].savefig(file)

# gets object properties
def getMaskProperties(segmentedImage3D, image3D_aligned, threshold=10):
    """
    get object properties from labeled image and formats
    centroids in NPY array

    Parameters
    ----------
    segmentedImage3D : NPY 3D array
        labeled 3D image.
    image3D_aligned : NPY 3D array
        pre-processed 3D image.

    Returns
    -------
    spots : NPY int64 array
        list of spots with the format: zyx

    """
    properties = regionprops(segmentedImage3D, intensity_image=image3D_aligned)
    centroids=[x.weighted_centroid for x in properties]
    sharpness=[float(x.filled_area/x.bbox_area) for x in properties]
    roundness1=[x.equivalent_diameter for x in properties]
    roundness2=[x.extent for x in properties]
    npix=[x.area for x in properties]
    sky=[0.0 for x in properties]
    peak=[x.max_intensity for x in properties]
    flux=[x.max_intensity/threshold for x in properties] # peak intensity over the detection threshold
    mag=[-2.5*np.log10(x) for x in flux] # -2.5 log10(flux)

    z=[x[0] for x in centroids]
    y=[x[1] for x in centroids]
    x=[x[2] for x in centroids]
    spots = np.zeros((len(z),3))
    spots[:,0]=z
    spots[:,1]=y
    spots[:,2]=x
    spots=spots.astype('int64')

    return (
            spots,
            sharpness,
            roundness1,
            roundness2,
            npix,
            sky,
            peak,
            flux,
            mag,
            )

(
    spots,
    sharpness,
    roundness1,
    roundness2,
    npix,
    sky,
    peak,
    flux,
    mag,
    ) = getMaskProperties(segmentedImage3D, image3D,threshold = threshold_over_std)

z,window=40,5
plt.imshow(np.sum(image3D[z-window:z+window,:,:],axis=0),cmap='Greys',vmax=.8)
selection=np.abs(spots[:,0]-z)<window
color=np.array(flux)
plt.scatter(spots[selection,2],spots[selection,1],c=color[selection],marker='+',alpha=.7,cmap='jet')

#%% LEFTOVERS

#%% breaks in blocks and represents in 3D

blockSizeXY=128
numPlanes = image3D_test.shape[0]
blockSize = (numPlanes, blockSizeXY, blockSizeXY)
images = [image3D_test,segmentedImage3D]
blocks = [view_as_blocks(x, block_shape=blockSize).squeeze() for x in images]

fig_output = []
for axis in range(3):
    fig_output.append(combinesBlocksImageByReprojection(blocks[0], blocks[1],axis1=axis))

images = [x[0][:,:,0] for x in fig_output]

fig = display3D_assembled(images, localizations = [spots_subpixel,spots])

#%%  first test ti make a 3D plot. Problem is it depends on image size... so will ignore.

# definitions for the axes
left, width = 0.02, 0.8
bottom, height = 0.15, 0.8
bottom_h, left_h = 0.7, 0.65
left_2, width_2 = 0.12, 0.6

rect_XY = [left, bottom, width, height]
rect_ZX = [left_h, bottom, 0.25, height]
rect_YZ = [left_2, left, width_2, 0.15]

fig = plt.figure()
# fig.set_size_inches((20, 20))
ax=list()
ax.append(plt.axes(rect_XY))
ax.append(plt.axes(rect_ZX))
ax.append(plt.axes(rect_YZ))

orientations=[0,1,0]
for axis, output, i, orientation in zip(ax, fig_output, range(3),orientations):
    img = output[0]
    if orientation==0:
        axis.imshow(img)
    else:
        axis.imshow(img.transpose((1,0,2)))
    axis.axes.xaxis.set_visible(False)
    axis.axes.yaxis.set_visible(False)

