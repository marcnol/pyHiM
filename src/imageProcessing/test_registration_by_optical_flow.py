#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 09:59:56 2023

@author: marcnol
"""

import numpy as np
from matplotlib import pyplot as plt
from skimage.color import rgb2gray
from skimage.data import stereo_motorcycle, vortex
from skimage.transform import warp
from skimage.registration import optical_flow_tvl1, optical_flow_ilk
import os
from skimage import io
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground

def remove_inhomogeneous_background(im, background_sigma):
    sigma_clip = SigmaClip(sigma=background_sigma)
    bkg_estimator = MedianBackground()
    bkg = Background2D(
        im,
        (64, 64),
        filter_size=(3, 3),
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
    )

    return im - bkg.background


#%%
# --- Load the sequence
 
folder= "/home/marcnol/grey/users/marcnol/test_HiM/testDataset"
im_list = ["scan_001_RT29_001_ROI_converted_decon_ch00.tif","scan_001_RT37_001_ROI_converted_decon_ch00.tif"]
im_list = [folder+os.sep+x for x in im_list]

img = [io.imread(x).squeeze() for x in im_list]

img = [np.max(x, axis=0) for x in img]


#%% 
background_sigma=3
img = [remove_inhomogeneous_background(x,background_sigma) for x in img]

image0, image1= img[0], img[1]

# --- Compute the optical flow
v, u = optical_flow_tvl1(image0, image1)

#%% Use the estimated optical flow for registration

nr, nc = image0.shape

row_coords, col_coords = np.meshgrid(np.arange(nr), np.arange(nc),
                                     indexing='ij')

image1_warp = warp(image1, np.array([row_coords + v, col_coords + u]),
                   mode='edge')


# build an RGB image with the unregistered sequence
seq_im = np.zeros((nr, nc, 3))
seq_im[..., 0] = image1
seq_im[..., 1] = image0
seq_im[..., 2] = image0

# build an RGB image with the registered sequence
reg_im = np.zeros((nr, nc, 3))
reg_im[..., 0] = image1_warp
reg_im[..., 1] = image0
reg_im[..., 2] = image0

# build an RGB image with the registered sequence
target_im = np.zeros((nr, nc, 3))
target_im[..., 0] = image0
target_im[..., 1] = image0
target_im[..., 2] = image0

# --- Show the result

fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(5, 10))

ax0.imshow(seq_im)
ax0.set_title("Unregistered sequence")
ax0.set_axis_off()

ax1.imshow(reg_im)
ax1.set_title("Registered sequence")
ax1.set_axis_off()

ax2.imshow(target_im)
ax2.set_title("Target")
ax2.set_axis_off()

fig.tight_layout()