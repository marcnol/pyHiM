#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:36:20 2020

@author: marcnol
"""


import numpy as np
import os
from alignBarcodesMasks import plotDistanceHistograms,plotMatrix



    
fileNamMatrix='/home/marcnol/Documents/Experiment_19/buildsPWDmatrix'+os.sep+'buildsPWDmatrix_HiMscMatrix.npy'
fileNameBarcodes='/home/marcnol/Documents/Experiment_19/buildsPWDmatrix'+os.sep+'buildsPWDmatrix_uniqueBarcodes.ecsv'

SCmatrixCollated=np.load(fileNamMatrix)
uniqueBarcodes = np.loadtxt(fileNameBarcodes).astype(int)

pixelSize=0.1



#%% plots distance matrix
plotMatrix(SCmatrixCollated,uniqueBarcodes, pixelSize,cm='terrain',clim=1.4) # twilight_shifted_r

#%% plots histograms
plotDistanceHistograms(SCmatrixCollated,pixelSize)

#%% inverse matrix
# plots inverse distance matrix
plotMatrix(np.reciprocal(SCmatrixCollated),uniqueBarcodes, pixelSize,cm='terrain',clim=.1, figtitle='Inverse PWD',cmtitle='inverse distance, 1/nm') # twilight_shifted_r

#%% contact prpbability

threshold=0.25
SCthresholdMatrix=SCmatrixCollated<threshold


nX = nY =  SCmatrixCollated.shape[0]
nCells = SCmatrixCollated.shape[2]
SCmatrix=np.zeros((nX,nY))

for i in range(nX):
    for j in range(nY):
        if i!=j:
            #print('Printing [{}:{}]'.format(i,j))
            #axs[i,j].hist(pixelSize*SCmatrixCollated[i,j,:],bins=bins)
            distanceDistribution = pixelSize*SCmatrixCollated[i,j,:]
            probability = len(np.nonzero(distanceDistribution<threshold)[0])/nCells
            SCmatrix[i,j] = probability
            
cScale=SCmatrix.max()/10            
plotMatrix(SCmatrix,uniqueBarcodes, pixelSize,cm='terrain',clim=cScale, figtitle='HiM counts',cmtitle='probability',nCells=nCells) # twilight_shifted_r
            