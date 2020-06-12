#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:03:05 2020

@author: marcnol
"""


from scipy.io import loadmat
import os
import numpy as np


rootFolder='/home/marcnol/data/Experiment_Julian'

# loads SC matrix
data=loadmat(rootFolder+os.sep+'test_python.mat')
fileNameBarcodes= rootFolder+os.sep+'barcodes.ecsv'

SCmatrixCollated=data['distanceMatrixCumulative']

# loads barcodes
uniqueBarcodes=[]
uniqueBarcodes.append(np.loadtxt(fileNameBarcodes).astype(int))

# loads cell attributes
cellAttributesMatrix = data['cellAttributesMatrix']
SClabeledCollated = cellAttributesMatrix[0,:]