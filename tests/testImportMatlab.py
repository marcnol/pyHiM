#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:03:05 2020

@author: marcnol
"""


from scipy.io import loadmat
import os
import numpy as np
from astropy.table import Table, vstack

rootFolder = "/home/marcnol/data/Experiment_Julian"


def loadsSCdataMATLAB(ListData, datasetName, p):

    print("Dataset to load: {}\n\n".format(list(ListData.keys())[0]))

    SCmatrixCollated, uniqueBarcodes = [], []
    buildsPWDmatrixCollated, runName, SClabeledCollated = [], [], []

    for rootFolder in ListData[datasetName]["Folders"]:
        # [finds sample name]
        runName.append(os.path.basename(os.path.dirname(rootFolder)))
        # [loads and accumulates barcodes and scHiM matrix]
        fileNameMatrix = rootFolder + os.sep + "HiMscMatrix.mat"
        fileNameBarcodes = rootFolder + os.sep + "buildsPWDmatrix_uniqueBarcodes.ecsv"

        # loads SC matrix
        if os.path.exists(fileNameMatrix):
            data = loadmat(fileNameMatrix)
            SCmatrix1 = data["distanceMatrixCumulative"]
            SCmatrixCollated.append(SCmatrix1)
        else:
            print("*** Error: could not find {}".format(fileNameMatrix))

        # loads barcodes
        if os.path.exists(fileNameBarcodes):
            uniqueBarcodes.append(np.loadtxt(fileNameBarcodes).astype(int))
        else:
            print("*** Error: could not find {}".format(fileNameBarcodes))

        # loads cell attributes
        cellAttributesMatrix = data["cellAttributesMatrix"]
        ResultsTable = cellAttributesMatrix[0, :]

        SClabeled = np.zeros(len(ResultsTable))
        indexCellsWithLabel = [iRow for iRow, Row in enumerate(ResultsTable) if Row > 0]
        SClabeled[indexCellsWithLabel] = 1
        SClabeledCollated.append(SClabeled)

        print("\n>>>Merging rootFolder: {}".format(rootFolder))
        print("Cells added after merge: {}\n".format(SCmatrix1.shape[2]))

    print("{} datasets loaded\n".format(len(SCmatrixCollated)))

    return (
        SCmatrixCollated,
        uniqueBarcodes,
        runName,
        SClabeledCollated,
    )
