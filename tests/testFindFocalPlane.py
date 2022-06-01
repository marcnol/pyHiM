#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:46:43 2021

@author: marcnol
"""
import os

import numpy as np
from numpy.testing import assert_array_almost_equal
import matplotlib.pyplot as plt

from skimage.filters import threshold_local, gaussian
from skimage.util.apply_parallel import apply_parallel
from datetime import datetime
from skimage import io
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
    calculate_focus_per_block,
    _reinterpolate_focal_plane,
    gaussian,
    )
import cv2


if "atlantis" in os.uname()[1]:
    root_folder="/home/marcnol/data/Embryo_debug_dataset/test_dataset/"
    file = root_folder+'scan_001_RT31_001_ROI_converted_decon_ch01.tif'
else:
    root_folder = "/home/marcnol/grey/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/"
    file = root_folder+'scan_001_RT27_001_ROI_converted_decon_ch01.tif'
    # file = root_folder+'scan_001_RT41_001_ROI_converted_decon_ch01.tif'
imageRaw = io.imread(file).squeeze()


#%% testing several measures that can be useful to determine focal plane
raw_images=[imageRaw[i,:,:] for i in range(imageRaw.shape[0])]

std= [np.std(img) for img in raw_images]
signal= [np.max(img) for img in raw_images]
laplacian_variance = [cv2.Laplacian(img, cv2.CV_64F).var() for img in raw_images]


std=(std-min(std))/max(std)
signal=signal/max(signal)
laplacian_variance  = laplacian_variance/max(laplacian_variance)

# the laplacian variance seems to be the best
plt.plot(signal,'o-b')
plt.plot(std,'o-g')
plt.plot(laplacian_variance,'o-r')


#%% implement function to fit the laplacian variance in 1D to get best estimate of focal plane and range

def fit1DGaussian(x,y,title='',verbose=True):
    from sherpa.data import Data1D
    from sherpa.plot import DataPlot
    from sherpa.models.basic import Gauss1D
    from sherpa.stats import LeastSq
    from sherpa.optmethods import LevMar
    from sherpa.fit import Fit


    d=Data1D('laplacianProfile',x,y)
    # print(d)

    # dplot=DataPlot()
    # dplot.prepare(d)
    # dplot.plot()

    Gauss1Dmodel = Gauss1D()
    opt = LevMar()
    stat = LeastSq()

    gFit = Fit(d,Gauss1Dmodel,stat=stat, method=opt)
    fit_result = gFit.fit()

    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)

    if fit_result.succeeded:
        if verbose:
            print("<<Fitting successful>>")
            print(fit_result.format())

            ax.plot(d.x,d.y,'ko',label='data')
            ax.plot(d.x,Gauss1Dmodel(d.x),linewidth=2, label='gaussian fit')
            ax.legend(loc=2)
            ax.set_title(title)

        return dict(zip(fit_result.parnames,fit_result.parvals)), fig
    else:
        return {}


import scipy.optimize as spo

def fit_1d_gaussian_scipy(x,y,title='',verbose=True):

    fit_result={}

    try:
        fitgauss = spo.curve_fit(gaussian, x, y)
        # print("Estimation of focal plane (px): ", int(fitgauss[0][1]))
        fit_result["gauss1d.pos"] = fitgauss[0][1]
        fit_result["gauss1d.ampl"] = fitgauss[0][0]
        fit_result["gauss1d.fwhm"] = 2.355*fitgauss[0][2]
    except RuntimeError:
        print("Warning, too many iterations")
        return {}
        
    if verbose:
        fig=plt.figure()
        ax = fig.add_subplot(1,1,1)

        print("<<Fitting successful>>")

        ax.plot(x,y,'ko',label='data')
        ax.plot(x,gaussian(x,fitgauss[0][0],fitgauss[0][1],fitgauss[0][2]),linewidth=2, label='gaussian fit')
        ax.legend(loc=2)
        ax.set_title(title)

    return fit_result, fig
    
    
y = laplacian_variance
x = range(len(laplacian_variance))
fit_result,fig2 = fit1DGaussian(x,y,title='z-profile sherpa',verbose=True)
focal_plane = fit_result['gauss1d.pos']

fitResult2,fig3 = fit_1d_gaussian_scipy(x,y,title='z-profile scipy',verbose=True)

#%% Partial recoding using block decomposition using existing functions from imageProcessing.py
from scipy.stats import sigmaclip

focal_plane_matrix, fwhm, block = calculate_focus_per_block(imageRaw, block_size_xy=256)
focal_planes_to_process = focal_plane_matrix[~np.isnan(focal_plane_matrix)]

focal_plane,_,_ = sigmaclip(focal_planes_to_process,high = 3, low=3)

fig=plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax1.imshow(focal_plane_matrix ,cmap='coolwarm',vmax=60,vmin=0)
ax2 = fig.add_subplot(2,1,2)
ax2.hist(focal_plane)

# NEED TO DO A GAUSSIAN FIT TO RECOVER BEST FIT
# ALSO GET ESTIMATE FOR FWHM !
print("Focal plane: {}".format(np.mean(focal_plane)))

#%% Final coding into a reused function from imageProcessing.py

focal_plane_matrix, z_range, block = _reinterpolate_focal_plane(imageRaw,block_size_xy = 256, window=10)
fig=plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.imshow(focal_plane_matrix ,cmap='coolwarm',vmax=60,vmin=0)
