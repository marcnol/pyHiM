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
    _removesInhomogeneousBackground2D,
    _removesInhomogeneousBackground,
    imageAdjust,
    _segments3DvolumesByThresholding,
    savesImageAsBlocks,
    display3D,
    combinesBlocksImageByReprojection,
    display3D_assembled,
    appliesXYshift3Dimages,
    calculatesFocusPerBlock,
    _reinterpolatesFocalPlane,
    gaussian,
    )
import cv2


if "atlantis" in os.uname()[1]:
    rootFolder="/home/marcnol/data/Embryo_debug_dataset/test_dataset/"
    file = rootFolder+'scan_001_RT31_001_ROI_converted_decon_ch01.tif'
else:
    rootFolder = "/home/marcnol/grey/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/"
    file = rootFolder+'scan_001_RT27_001_ROI_converted_decon_ch01.tif'
    # file = rootFolder+'scan_001_RT41_001_ROI_converted_decon_ch01.tif'
imageRaw = io.imread(file).squeeze()


#%% testing several measures that can be useful to determine focal plane
rawImages=[imageRaw[i,:,:] for i in range(imageRaw.shape[0])]

std= [np.std(img) for img in rawImages]
signal= [np.max(img) for img in rawImages]
LaplacianVariance = [cv2.Laplacian(img, cv2.CV_64F).var() for img in rawImages]


std=(std-min(std))/max(std)
signal=signal/max(signal)
LaplacianVariance  = LaplacianVariance/max(LaplacianVariance)

# the laplacian variance seems to be the best
plt.plot(signal,'o-b')
plt.plot(std,'o-g')
plt.plot(LaplacianVariance,'o-r')


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
    fitResult = gFit.fit()

    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)

    if fitResult.succeeded:
        if verbose:
            print("<<Fitting successful>>")
            print(fitResult.format())

            ax.plot(d.x,d.y,'ko',label='data')
            ax.plot(d.x,Gauss1Dmodel(d.x),linewidth=2, label='gaussian fit')
            ax.legend(loc=2)
            ax.set_title(title)

        return dict(zip(fitResult.parnames,fitResult.parvals)), fig
    else:
        return dict()


import scipy.optimize as spo

def fit1DGaussian_scipy(x,y,title='',verbose=True):

    fitResult=dict()

    try:
        fitgauss = spo.curve_fit(gaussian, x, y)
        # print("Estimation of focal plane (px): ", int(fitgauss[0][1]))
        fitResult["gauss1d.pos"] = fitgauss[0][1]
        fitResult["gauss1d.ampl"] = fitgauss[0][0]
        fitResult["gauss1d.fwhm"] = 2.355*fitgauss[0][2]
    except RuntimeError:
        print("Warning, too many iterations")
        return dict()
        
    if verbose:
        fig=plt.figure()
        ax = fig.add_subplot(1,1,1)

        print("<<Fitting successful>>")

        ax.plot(x,y,'ko',label='data')
        ax.plot(x,gaussian(x,fitgauss[0][0],fitgauss[0][1],fitgauss[0][2]),linewidth=2, label='gaussian fit')
        ax.legend(loc=2)
        ax.set_title(title)

    return fitResult, fig
    
    
y = LaplacianVariance
x = range(len(LaplacianVariance))
fitResult,fig2 = fit1DGaussian(x,y,title='z-profile sherpa',verbose=True)
focalPlane = fitResult['gauss1d.pos']

fitResult2,fig3 = fit1DGaussian_scipy(x,y,title='z-profile scipy',verbose=True)

#%% Partial recoding using block decomposition using existing functions from imageProcessing.py
from scipy.stats import sigmaclip

focalPlaneMatrix, fwhm, block = calculatesFocusPerBlock(imageRaw, blockSizeXY=256)
focalPlanes2Process = focalPlaneMatrix[~np.isnan(focalPlaneMatrix)]

focalPlane,_,_ = sigmaclip(focalPlanes2Process,high = 3, low=3)

fig=plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax1.imshow(focalPlaneMatrix ,cmap='coolwarm',vmax=60,vmin=0)
ax2 = fig.add_subplot(2,1,2)
ax2.hist(focalPlane)

# NEED TO DO A GAUSSIAN FIT TO RECOVER BEST FIT
# ALSO GET ESTIMATE FOR FWHM !
print("Focal plane: {}".format(np.mean(focalPlane)))

#%% Final coding into a reused function from imageProcessing.py

focalPlaneMatrix, zRange, block = _reinterpolatesFocalPlane(imageRaw,blockSizeXY = 256, window=10)
fig=plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.imshow(focalPlaneMatrix ,cmap='coolwarm',vmax=60,vmin=0)
