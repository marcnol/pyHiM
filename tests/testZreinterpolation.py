#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 10:28:24 2021

@author: marcnol

test several methods for reducing the number of Z-planes in an image

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
    _remove_inhomogeneous_background_2d,
    _remove_inhomogeneous_background,
    image_adjust,
    _segment_3d_volumes_by_thresholding,
    save_image_as_blocks,
    display_3d,
    combine_blocks_image_by_reprojection,
    display_3d_assembled,
    apply_xy_shift_3d_images,
    )
from photutils import detect_sources
from astropy.convolution import Gaussian2DKernel
from skimage import measure
from astropy.visualization import SqrtStretch, simple_norm
from skimage.util.shape import view_as_blocks

from scipy import ndimage as ndi

from skimage.segmentation import watershed
from skimage.feature import peak_local_max

root_folder="/home/marcnol/data/Embryo_debug_dataset/test_dataset/"
file = root_folder+'scan_001_RT27_001_ROI_converted_decon_ch01.tif'

image_3d = io.imread(file).squeeze()

#%% takes one every n planes
from datetime import datetime

def _remove_z_planes(image_3d, z_range):

    output = np.zeros(image_3d.shape)
    for i,index in enumerate(z_range):
        output[i,:,:] = image_3d[index,:,:]

    return output

def _interpolate_z_planes(image_3d, z_range):

    from scipy.interpolate import interpn

    output = np.zeros(image_3d.shape)

    # need to code using interpn
    output = image_3d

    return output

def reinterpolate_z(image_3d, z_range,mode='remove'):

    if 'interpolate' in mode:
        output = _interpolate_z_planes(image_3d, z_range)
    elif 'remove' in mode:
        output = _remove_z_planes(image_3d, z_range)

    return output


numberZplanes=image_3d.shape[0]
binning = 2
z_range = range(0,numberZplanes,binning)

begin_time = datetime.now()
image3Dr = reinterpolate_z(image_3d, z_range,mode='remove')
print("Elapsed time: {}".format(datetime.now() - begin_time))


images = [image_3d, image3Dr,np.zeros(image_3d.shape)]
axis = 0

images2d = [np.sum(img,axis=axis) for img in images]
images2d = [img/img.max() for img in images2d]
rgb = np.zeros((images2d[0].shape[0],images2d[0].shape[1],3))

for i, img in enumerate(images2d):
    rgb[:,:,i] = img

plt.imshow(rgb)
