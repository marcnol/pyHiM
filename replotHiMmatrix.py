#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:36:20 2020

@author: marcnol
"""


import numpy as np
import os
from alignBarcodesMasks import plotDistanceHistograms,plotMatrix
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import LeaveOneOut
from alignBarcodesMasks import distributionMaximumKernelDensityEstimation,calculateContactProbabilityMatrix


#%% Lists and loads datasets from different embryos

# HiRes doc TAD
ListRootFolders=['/mnt/disk2/marcnol/data/Experiment_19/026_Embryo/buildsPWDmatrix',\
                 '/mnt/disk2/marcnol/data/Experiment_19/009_Embryo/buildsPWDmatrix',\
                 '/mnt/disk2/marcnol/data/Experiment_19/006_Embryo/buildsPWDmatrix']
fileTag='wt_docTAD'    
    
'''    
# HiRes doc locus    
ListRootFolders=[\
                 #'/mnt/disk2/marcnol/data/Experiment_3/019_Embryo/buildsPWDmatrix',\
                 '/mnt/disk2/marcnol/data/Experiment_3/007_Embryo/buildsPWDmatrix',\
                 '/mnt/disk2/marcnol/data/Experiment_3/016_Embryo/buildsPWDmatrix',\
                 #'/mnt/disk2/marcnol/data/Experiment_3/000_Embryo/buildsPWDmatrix'\
                 '/mnt/disk2/marcnol/data/Experiment_3/001_Embryo/buildsPWDmatrix',\
                 '/mnt/disk2/marcnol/data/Experiment_3/002_Embryo/buildsPWDmatrix',\
                 '/mnt/disk2/marcnol/data/Experiment_3/003_Embryo/buildsPWDmatrix',\
                 '/mnt/disk2/marcnol/data/Experiment_3/004_Embryo/buildsPWDmatrix',\
                 '/mnt/disk2/marcnol/data/Experiment_3/005_Embryo/buildsPWDmatrix'\
                 ]
fileTag='wt_HresDocLocus'    
'''
# HiRes doc TAD -zld
ListRootFolders=[\
              '/mnt/disk2/marcnol/data/Experiment_5/0_Embryo/buildsPWDmatrix',\
             '/mnt/disk2/marcnol/data/Experiment_5/1_Embryo/buildsPWDmatrix',\
             # '/mnt/disk2/marcnol/data/Experiment_5/2_Embryo/buildsPWDmatrix',\
             '/mnt/disk2/marcnol/data/Experiment_5/3_Embryo/buildsPWDmatrix',\
             '/mnt/disk2/marcnol/data/Experiment_5/4_Embryo/buildsPWDmatrix',\
             '/mnt/disk2/marcnol/data/Experiment_5/5_Embryo/buildsPWDmatrix',\
             '/mnt/disk2/marcnol/data/Experiment_5/6_Embryo/buildsPWDmatrix',\
             '/mnt/disk2/marcnol/data/Experiment_5/7_Embryo/buildsPWDmatrix',\
             '/mnt/disk2/marcnol/data/Experiment_5/8_Embryo/buildsPWDmatrix',\
             '/mnt/disk2/marcnol/data/Experiment_5/9_Embryo/buildsPWDmatrix'\
             ]
fileTag='zld_docTAD'    

embryoStage = [\
                #'fucked-up',\
                'nc14',\
                'nc14',\
                #'mitosis'\
                ]
    
outputFolder='/mnt/disk2/marcnol/data/Experiment_3/'
os.chdir(outputFolder)

pixelSize=0.1
SCmatrixCollated, uniqueBarcodes = [], []

for rootFolder in ListRootFolders:
    
    fileNamMatrix=rootFolder+os.sep+'buildsPWDmatrix_HiMscMatrix.npy'
    fileNameBarcodes=rootFolder+os.sep+'buildsPWDmatrix_uniqueBarcodes.ecsv'
    
    SCmatrixCollated.append(np.load(fileNamMatrix))
    uniqueBarcodes.append(np.loadtxt(fileNameBarcodes).astype(int))

print('{} datasets loaded'.format(len(SCmatrixCollated)))


#%% plots distance matrix for each dataset
for iSCmatrixCollated, iuniqueBarcodes in zip(SCmatrixCollated,uniqueBarcodes):
    plotMatrix(iSCmatrixCollated,iuniqueBarcodes, pixelSize,cm='terrain',clim=1.4,mode='KDE') # twilight_shifted_r 1.4, mode: median KDE

#%% plots histograms for each dataset
for iSCmatrixCollated, iuniqueBarcodes in zip(SCmatrixCollated,uniqueBarcodes):
    plotDistanceHistograms(iSCmatrixCollated,pixelSize,mode='KDE',limitNplots=15)

#%% plots inverse distance matrix for each dataset
for iSCmatrixCollated, iuniqueBarcodes in zip(SCmatrixCollated,uniqueBarcodes):
    # plotMatrix(np.reciprocal(iSCmatrixCollated),iuniqueBarcodes, pixelSize,cm='terrain',clim=.1, figtitle='Inverse PWD',cmtitle='inverse distance, 1/nm') # twilight_shifted_r
    plotMatrix(iSCmatrixCollated,
               iuniqueBarcodes, 
               pixelSize,
               cm='terrain',
               clim=6, 
               figtitle='Inverse PWD',
               cmtitle='inverse distance, 1/nm',
               inverseMatrix=True, 
               mode='KDE') # twilight_shifted_r

#%% Plots contact prpbability matrices for each dataset

threshold=.25

for iSCmatrixCollated, iuniqueBarcodes in zip(SCmatrixCollated,uniqueBarcodes):
    SCmatrix, nCells=calculateContactProbabilityMatrix(iSCmatrixCollated,iuniqueBarcodes,pixelSize,threshold)
    cScale=SCmatrix.max()/15            
    plotMatrix(SCmatrix,iuniqueBarcodes, pixelSize,cm='terrain',clim=cScale, figtitle='HiM counts',cmtitle='probability',nCells=nCells) # twilight_shifted_r
    
#%% combines matrices from different embryos and calculates integrated contact probability matrix
threshold=0.25
nBarcodes=uniqueBarcodes[0].shape[0]
SCmatrixAllDatasets=[] #np.zeros((nBarcodes,nBarcodes))

for iSCmatrixCollated, iuniqueBarcodes in zip(SCmatrixCollated,uniqueBarcodes):
    if len(SCmatrixAllDatasets)>0:
        SCmatrixAllDatasets=np.concatenate((SCmatrixAllDatasets,iSCmatrixCollated),axis=2)
    else:
        SCmatrixAllDatasets=iSCmatrixCollated
    
    commonSetUniqueBarcodes = iuniqueBarcodes
    
SCmatrix, nCells=calculateContactProbabilityMatrix(SCmatrixAllDatasets,commonSetUniqueBarcodes,pixelSize,threshold)
cScale=SCmatrix.max()/10            

plotMatrix(SCmatrix,iuniqueBarcodes, pixelSize,cm='terrain',clim=cScale, figtitle='HiM counts',cmtitle='probability',nCells=nCells) # twilight_shifted_r
    
np.savetxt(outputFolder+os.sep+'CombinedMatrix'+fileTag+'.dat', SCmatrix, fmt='%.4f', \
           delimiter=' ', newline='\n', header='Combined contact probability matrix', \
           footer='', comments='# ', encoding=None)
    
np.savetxt(outputFolder+os.sep+'UniqueBarcodes'+fileTag+'.dat', iuniqueBarcodes, fmt='%.4f', \
           delimiter=' ', newline='\n', header='unique barcodes', \
           footer='', comments='# ', encoding=None)
    
    
    
    
    
    
    
    
        