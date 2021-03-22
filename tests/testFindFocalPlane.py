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
    )
import cv2


if "atlantis" in os.uname()[1]:
    rootFolder="/home/marcnol/data/Embryo_debug_dataset/test_dataset/"
    file = rootFolder+'scan_001_RT31_001_ROI_converted_decon_ch01.tif'
else:
    rootFolder = "/home/marcnol/grey/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/"
    file = rootFolder+'scan_001_RT27_001_ROI_converted_decon_ch01.tif'
    file = rootFolder+'scan_001_RT41_001_ROI_converted_decon_ch01.tif'
imageRaw = io.imread(file).squeeze()

# autoscales exposures
# image3D = exposure.rescale_intensity(imageRaw, out_range=(0, 1))

#%%
rawImages=[imageRaw[i,:,:] for i in range(imageRaw.shape[0])]

std= [np.std(img) for img in rawImages]
signal= [np.max(img) for img in rawImages]
LaplacianVariance = [cv2.Laplacian(img, cv2.CV_64F).var() for img in rawImages]


std=(std-min(std))/max(std)
signal=signal/max(signal)
LaplacianVariance  = LaplacianVariance/max(LaplacianVariance)
    
plt.plot(signal,'o-b')
plt.plot(std,'o-g')
plt.plot(LaplacianVariance,'o-r')

#%%

def fit1DGaussian(x,y,title='',verbose=True):
    from sherpa.data import Data1D
    from sherpa.plot import DataPlot
    from sherpa.models.basic import Gauss1D
    from sherpa.stats import LeastSq
    from sherpa.optmethods import LevMar
    from sherpa.fit import Fit

   
    d=Data1D('laplacianProfile',x,y)
    # print(d)
    
    dplot=DataPlot()
    dplot.prepare(d)
    dplot.plot()
    
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
    
y = LaplacianVariance
x = range(len(LaplacianVariance))
fitResult,fig2 = fit1DGaussian(x,y,title='z-profile',verbose=True)
focalPlane = fitResult['gauss1d.pos']