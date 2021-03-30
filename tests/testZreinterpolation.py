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

rootFolder="/home/marcnol/data/Embryo_debug_dataset/test_dataset/"
file = rootFolder+'scan_001_RT27_001_ROI_converted_decon_ch01.tif'

image3D = io.imread(file).squeeze()

#%% takes one every n planes
from datetime import datetime

def _removesZplanes(image3D, Zrange):

    output = np.zeros(image3D.shape)
    for i,index in enumerate(Zrange):
        output[i,:,:] = image3D[index,:,:]

    return output

def _interpolatesZplanes(image3D, Zrange):

    from scipy.interpolate import interpn

    output = np.zeros(image3D.shape)

    # need to code using interpn
    output = image3D

    return output

def reinterpolateZ(image3D, Zrange,mode='remove'):

    if 'interpolate' in mode:
        output = _interpolatesZplanes(image3D, Zrange)
    elif 'remove' in mode:
        output = _removesZplanes(image3D, Zrange)

    return output


numberZplanes=image3D.shape[0]
binning = 2
Zrange = range(0,numberZplanes,binning)

begin_time = datetime.now()
image3Dr = reinterpolateZ(image3D, Zrange,mode='remove')
print("Elapsed time: {}".format(datetime.now() - begin_time))


images = [image3D, image3Dr,np.zeros(image3D.shape)]
axis = 0

images2d = [np.sum(img,axis=axis) for img in images]
images2d = [img/img.max() for img in images2d]
RGB = np.zeros((images2d[0].shape[0],images2d[0].shape[1],3))

for i, img in enumerate(images2d):
    RGB[:,:,i] = img

plt.imshow(RGB)
