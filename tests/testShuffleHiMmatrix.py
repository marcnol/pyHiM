#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 20:48:16 2021

@author: marcnol
"""

import numpy as np
import matplotlib.pylab as plt

def shuffleMatrix(matrix, index):

    newSize = len(index)
    newMatrix = np.zeros((newSize, newSize))

    if newSize <= matrix.shape[0]:
        for i in range(newSize):
            for j in range(newSize):
                if index[i] < matrix.shape[0] and index[j] < matrix.shape[0]:
                    newMatrix[i, j] = matrix[index[i], index[j]]
    else:
        print("Error: shuffle size {} is larger than matrix dimensions {}".format(newSize, matrix.shape[0]))
        print("Shuffle: {} ".format(index))

    return newMatrix

def calculateContactProbabilityMatrix(iSCmatrixCollated, iuniqueBarcodes, pixelSize, threshold=0.25, norm="nCells"):

    nX = nY = iSCmatrixCollated.shape[0]
    nCells = iSCmatrixCollated.shape[2]
    SCmatrix = np.zeros((nX, nY))

    for i in range(nX):
        for j in range(nY):
            if i != j:
                distanceDistribution = pixelSize * iSCmatrixCollated[i, j, :]

                # normalizes # of contacts by the # of cells
                if norm == "nCells":
                    probability = len(np.nonzero(distanceDistribution < threshold)[0]) / nCells
                    # print('Using nCells normalisation')

                # normalizes # of contacts by the # of PWD detected in each bin
                elif norm == "nonNANs":
                    numberNANs = len(np.nonzero(np.isnan(distanceDistribution))[0])
                    if nCells == numberNANs:
                        probability = 1
                    else:
                        probability = len(np.nonzero(distanceDistribution < threshold)[0]) / (nCells - numberNANs)
                    # print('Using NonNANs normalisation {}'.format(nCells-numberNANs))
                SCmatrix[i, j] = probability

    return SCmatrix, nCells

fileName="/home/marcnol/data/Embryo_debug_dataset/run_markus/buildsPWDmatrix/buildsPWDmatrix_3D_HiMscMatrix.npy"
# fileName="/home/marcnol/data/Embryo_debug_dataset/run_markus_raw_3D/buildsPWDmatrix/buildsPWDmatrix_3D_HiMscMatrix.npy"
# fileName='/home/marcnol/data/Embryo_debug_dataset/run_markus_raw_3D/buildsPWDmatrix/buildsPWDmatrix_3D_HiMscMatrix.npy'

shuffle = [0, 7, 1, 8, 2, 9, 10, 17, 11, 18, 12, 19, 13, 3, 14, 4, 15, 5, 16, 6]

matrixSC =np.load(fileName)
pixelSize = 0.1
matrix,_ = calculateContactProbabilityMatrix(matrixSC, shuffle, pixelSize, threshold=0.25, norm="nCells")

# matrix = np.nanmean(matrixSC,axis=2)
matrix_shuffled = shuffleMatrix(matrix, shuffle)
plt.imshow(matrix_shuffled,cmap='coolwarm',vmax=.01)