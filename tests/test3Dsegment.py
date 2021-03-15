#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:05:02 2021

@author: marcnol

- tests several algorithms for 3D segmentation:
    - based on ASTROPY (plane by plane object-based segmentation and deblending)
    - blob_log from skimage

"""

import os, argparse, sys, glob
from datetime import datetime
import matplotlib.pylab as plt
from skimage import io
import numpy as np
from tifffile import imsave
from tqdm import tqdm, trange
from skimage import exposure
from imageProcessing  import (
    _removesInhomogeneousBackground2D,
    imageAdjust,
    _segments3DvolumesByThresholding,
    savesImageAsBlocks,
    display3D,
    combinesBlocksImageByReprojection,
    )
from photutils import detect_sources
from astropy.convolution import Gaussian2DKernel
from skimage import measure
from astropy.visualization import SqrtStretch, simple_norm
from skimage.util.shape import view_as_blocks

          
    
#%% loads and segments a file

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


#%% decompose output labeled image in blocks

fileName = rootFolder + "blockTest"
savesImageAsBlocks(segmentedImage3D,fileName,blockSizeXY=256)
    
#%% Deblend image in 3D by watersheding
    
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi

from skimage.segmentation import watershed
from skimage.feature import peak_local_max

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


#%% loads and display output of a segmentation run

folder = "/home/marcnol/grey/users/marcnol/models/StarDist3D/training3Dbarcodes/dataset1/"
file_image3D=folder+os.sep+"scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0.tif"

file_segmented=folder+os.sep+"scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0_segmented.tif"
image3D = io.imread(file_image3D).squeeze()
labeled = io.imread(file_segmented).squeeze()

display3D(image3D=image3D,labels=labeled,z=40, rangeXY=1000)

#%%

display3D(image3D=image3D,z=40, rangeXY=1000,norm=False)

#%% skimage

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

z = 40

rootFolder="/mnt/grey/DATA/users/marcnol/models/StarDist3D/training3Dbarcodes/dataset1/"
file = rootFolder+'scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0.tif'
threshold_over_std,sigma ,boxSize, filter_size=1, 3, (32, 32),(3, 3)
nlevels=64
contrast=0.001


print("loading image {}".format(os.path.basename(file)))

image3D = io.imread(file).squeeze()

image3D_test=image3D[z-5:z+5,:,:]

# segments in 3D using ASTROPY
binary, segmentedImage3D = _segments3DvolumesByThresholding(image3D_test,
                                                       threshold_over_std=threshold_over_std, 
                                                       sigma = 3, 
                                                       boxSize=(32, 32),
                                                       filter_size=(3, 3),
                                                       nlevels=nlevels,
                                                       contrast=contrast,
                                                       deblend3D=True)

#%%
# finds first estimate of centroids from 3D masks

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
spots_subpixel2 = spots_subpixel.copy()
# spots_subpixel2 = spots_subpixel2[:, [2, 1, 0]]

spots2 = spots.copy()
# spots2=spots2[:,[2,1,0]]

display3D(image3D=image3D_test,localizationsList = [spots2,spots_subpixel2],labels=segmentedImage3D,z=5, rangeXY=1000, norm=True,cmap='Greys')


#%%
blockSizeXY=128
numPlanes = image3D_test.shape[0]
blockSize = (numPlanes, blockSizeXY, blockSizeXY)
images = [image3D_test,segmentedImage3D]
blocks = [view_as_blocks(x, block_shape=blockSize).squeeze() for x in images]

fig_output = []
for axis in range(3):
    fig_output.append(combinesBlocksImageByReprojection(blocks[0], blocks[1],axis1=axis))

#%%
fig3 = plt.figure(constrained_layout=False)
ncols,nrows=5,5
widths = [1]*5
heights = [1]*5
gs = fig3.add_gridspec(ncols=ncols, nrows=nrows, width_ratios=widths,
                          height_ratios=heights)

fig3.set_size_inches((20, 20))
ax = [fig3.add_subplot(gs[0:-1,0:-1]), fig3.add_subplot(gs[4, 1:-1]), fig3.add_subplot(gs[1:-1,4])]

titles = ["Z-projection", "X-projection", "Y-projection"]
for axis, output, i in zip(ax, fig_output, range(3)):
    axis.imshow(output[0])
    axis.set_title(titles[i])