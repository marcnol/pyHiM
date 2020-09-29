#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 16:00:31 2020

@author: marcnol
"""
import numpy as np
import cv2
import os
import matplotlib.pyplot as plt
from datetime import datetime


## Functions

def alignCV2(im1,im2,warp_mode):
    
    # Find size of image1
    sz = im1.shape
    
    # Define 2x3 or 3x3 matrices and initialize the matrix to identity
    if warp_mode == cv2.MOTION_HOMOGRAPHY :
        warp_matrix = np.eye(3, 3, dtype=np.float32)
    else :
        warp_matrix = np.eye(2, 3, dtype=np.float32)
    
    # Specify the number of iterations.
    number_of_iterations = 5000;
    
    # Specify the threshold of the increment
    # in the correlation coefficient between two iterations
    termination_eps = 1e-10;
    
    # Define termination criteria
    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, number_of_iterations,  termination_eps)
    
    # Run the ECC algorithm. The results are stored in warp_matrix.
    try:
        cc, warp_matrix = cv2.findTransformECC(im1,im2,warp_matrix, warp_mode, criteria, inputMask=None, gaussFiltSize=1)
    except TypeError:
        cc, warp_matrix = cv2.findTransformECC(im1,im2,warp_matrix, warp_mode, criteria)

    return cc, warp_matrix

def applyCorrection(im2,warp_matrix):
    
    sz = im2.shape

    # Use warpAffine for Translation, Euclidean and Affine
    im2_aligned = cv2.warpAffine(im2, warp_matrix, (sz[1],sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP);

    return im2_aligned 

#%%
    
## main
rootFolder = "/home/marcnol/data/Embryo_debug_dataset/Experiment_18/zProject"
imFileName1=rootFolder+os.sep+"scan_001_RT27_001_ROI_converted_decon_ch00_2d.npy"
imFileName2=rootFolder+os.sep+"scan_001_RT37_001_ROI_converted_decon_ch00_2d.npy"

image1 = np.load(imFileName1).squeeze()
image2 = np.load(imFileName2).squeeze()

image1 = image1 / image1.max()
image2 = image2 / image2.max()

image1 = np.float32(image1)
image2 = np.float32(image2)

#%%

# Define the motion model
# model='afine'
model='translate'
# model='euclidean'

if model=='translate':
    warp_mode = cv2.MOTION_TRANSLATION
elif model=='afine':
    warp_mode = cv2.MOTION_AFFINE
elif model=='euclidean':
    warp_mode = cv2.MOTION_EUCLIDEAN

outputFileName = '/home/marcnol/Downloads/'+'test_'+model+'.png'

begin_time = datetime.now()

cc, warp_matrix = alignCV2(image1,image2,warp_mode)
im2_aligned = applyCorrection(image2,warp_matrix)
    
print('Done registering!')

print("Elapsed time: {}".format(datetime.now() - begin_time))

#%%

sz = image1.shape

fig, (ax1,ax2) = plt.subplots(1,2)
fig.set_size_inches((60, 30))

nullImage = np.zeros(sz)

RGB = np.dstack([image1, image2, nullImage ])
ax1.imshow(RGB)
ax1.axis("off")

RGB_aligned= np.dstack([image1, im2_aligned, nullImage ])
ax2.imshow(RGB_aligned)
ax2.axis("off")

plt.savefig(outputFileName)



