#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 08:40:30 2022

@author: marcnol

This script:
    - iterates over chromatin traces
        - calculates the pair-wise distances for each single-cell mask
        - outputs are:
            - Table with #cell #PWD #coordinates (e.g. buildsPWDmatrix_3D_order:0_ROI:1.ecsv)
            - NPY array with single cell PWD single cell matrices (e.g. buildsPWDmatrix_3D_HiMscMatrix.npy)
            - NPY array with barcode identities (e.g. buildsPWDmatrix_3D_uniqueBarcodes.ecsv)
            - the files with no "3D" tag contain data analyzed using 2D localizations.

    - Single-cell results are combined together to calculate:
        - Distribution of pairwise distance for each barcode combination
        - Ensemble mean pairwise distance matrix using mean of distribution
        - Ensemble mean pairwise distance matrix using Kernel density estimation
        - Ensemble Hi-M matrix using a predefined threshold
        - For each of these files, there is an image in PNG format saved. Images containing "3D" are for 3D other are for 2D.


"""

# =============================================================================
# IMPORTS
# =============================================================================

import glob, os, sys
import uuid
import re
import numpy as np
from tqdm.contrib import tzip
from tqdm import trange
import matplotlib.pyplot as plt

from sklearn.metrics import pairwise_distances

from astropy.table import Table, unique

from photutils.segmentation import SegmentationImage

from fileProcessing.fileManagement import (
    folders,
    writeString2File,
    printLog,
    getDictionaryValue,
)

from matrixOperations.HIMmatrixOperations import plotMatrix, plotDistanceHistograms, calculateContactProbabilityMatrix
from matrixOperations.build_traces import initialize_module
from matrixOperations.chromatin_trace_table import chromatin_trace_table


# to remove in a future version
import warnings

warnings.filterwarnings("ignore")

# =============================================================================
# CLASSES
# =============================================================================


class build_matrix:
    def __init__(self, param):

        self.param = param

        self.initialize_parameters()

        # initialize with default values
        self.pixelSize = [0.1,0.1,0.25] # default pixel size
        self.currentFolder = []

    def initialize_parameters(self):
        # initializes parameters from param

        self.tracing_method = getDictionaryValue(self.param.param["buildsPWDmatrix"], "tracing_method",  default="masking")
        self.zBinning = getDictionaryValue(self.param.param["acquisition"], "zBinning", default=1)
        self.pixelSizeXY = getDictionaryValue(self.param.param["acquisition"], "pixelSizeXY", default=0.1)
        self.pixelSizeZ_0 = getDictionaryValue(self.param.param["acquisition"], "pixelSizeZ", default=0.25)
        self.pixelSizeZ = self.zBinning * self.pixelSizeZ_0
        self.availableMasks = getDictionaryValue(self.param.param["buildsPWDmatrix"], "masks2process",  default={"nuclei":"DAPI"})
        self.logNameMD = self.param.param["fileNameMD"]
        self.mask_expansion = getDictionaryValue(self.param.param["buildsPWDmatrix"], "mask_expansion", default=8)
        self.availableMasks = self.param.param["buildsPWDmatrix"]["masks2process"]

    def calculatesPWDsingleMask(self, ROI, CellID, groupKeys, x, y, z):
        """
        Calculates PWD between barcodes detected in a given mask. For this:
            - converts xyz pixel coordinates into nm using self.pixelSize dictionary
            - calculates pair-wise distance matrix in nm
            - converts it into pixel units using self.pixelSize['x'] as an isotropic pixelsize.

        Parameters
        ----------
        ROI : string
            ROI used
        CellID: string
            ID of the cell
        x_uncorrected: float
            x coordinates uncorrected
        y_uncorrected: float
            y coordinates uncorrected
        z_uncorrected: float
            z coordinates uncorrected

        Returns
        -------
        Returns pairwise distance matrix between corrected barcodes in isotropic pixel units

        """
        R_nm = self.buildsVector(groupKeys, x, y, z)

        P = pairwise_distances(R_nm)

        P = P/self.pixelSize["x"]

        return P

    def buildsdistanceMatrix(self, mode="mean"):
        """
        Builds pairwise distance matrix from a coordinates table

        Parameters
        ----------
        mode : string, optional
            The default is "mean": calculates the mean distance if there are several combinations possible.
            "min": calculates the minimum distance if there are several combinations possible.
            "last": keeps the last distance calculated

        Returns
        -------
        self.SCmatrix the single-cell PWD matrix
        self.meanSCmatrix the ensamble PWD matrix (mean of SCmatrix without nans)
        self.uniqueBarcodes list of unique barcodes

        """
        # detects number of unique traces from trace table
        numberMatrices = len(unique(self.trace_table.data,keys='Trace_ID'))

        # finds unique barcodes from trace table
        uniqueBarcodes = unique(self.trace_table.data,keys='Barcode #')['Barcode #'].data
        numberUniqueBarcodes = uniqueBarcodes.shape[0]

        printLog(f"$ Found {numberUniqueBarcodes} barcodes and {numberMatrices} traces.","INFO")

        # Initializes SCmatrix
        SCmatrix = np.zeros((numberUniqueBarcodes, numberUniqueBarcodes, numberMatrices))
        SCmatrix[:] = np.NaN

        # loops over traces

        data_traces = self.trace_table.data.group_by("Trace_ID")
        for trace, trace_id in zip(data_traces.groups, data_traces.groups.keys):

#        for iCell, scPWDitem in zip(range(numberMatrices), self.SCdistanceTable):
            barcodes2Process = trace["Barcode #"].data

-----

            # loops over barcodes detected in cell mask: barcode1
            for barcode1, ibarcode1 in zip(barcodes2Process, range(len(barcodes2Process))):
                indexBarcode1 = np.nonzero(uniqueBarcodes == barcode1)[0][0]

                # loops over barcodes detected in cell mask: barcode2
                for barcode2, ibarcode2 in zip(barcodes2Process, range(len(barcodes2Process))):
                    indexBarcode2 = np.nonzero(uniqueBarcodes == barcode2)[0][0]

                    if barcode1 != barcode2:

                        # attributes distance from the PWDmatrix field in the scPWDitem table
                        newdistance = scPWDitem["PWDmatrix"][ibarcode1][ibarcode2]

                        # inserts value into SCmatrix
                        if mode == "last":
                            SCmatrix[indexBarcode1][indexBarcode2][iCell] = newdistance
                        elif mode == "mean":
                            SCmatrix[indexBarcode1][indexBarcode2][iCell] = np.nanmean(
                                [newdistance, SCmatrix[indexBarcode1][indexBarcode2][iCell],]
                            )
                        elif mode == "min":
                            SCmatrix[indexBarcode1][indexBarcode2][iCell] = np.nanmin(
                                [newdistance, SCmatrix[indexBarcode1][indexBarcode2][iCell],]
                            )

        self.SCmatrix = SCmatrix
        self.meanSCmatrix = np.nanmean(SCmatrix, axis=2)
        # self.uniqueBarcodes = uniqueBarcodes


    def calculatesNmatrix(self,SCmatrix):

        numberCells = SCmatrix.shape[2]

        if numberCells > 0:
            Nmatrix = np.sum(~np.isnan(SCmatrix), axis=2)
        else:
            numberBarcodes = SCmatrix.shape[0]
            Nmatrix = np.zeros((numberBarcodes, numberBarcodes))

        return Nmatrix


    def plotsAllmatrices(self, file):
        """
        Plots all matrices after analysis

        Parameters
        ----------
        SCmatrixCollated : npy array
            PWD matrix for single cells.
        Nmatrix : npy array
            2d matrix with number of measurements per barcode combination.
        uniqueBarcodes : npy array
            barcode identities.
        pixelSize : npy array
            pixelsize in um.
        numberROIs : int
            self explanatory.
        outputFileName : str
            self explanatory.
        logNameMD : str
            Markdown filename.
        localizationDimension : int
            indicates dimension of barcode localization.

        Returns
        -------
        None.

        """
        # adapts clim depending on whether 2 or 3 dimensions are used for barcode localizations
        if localizationDimension == 2:
            clim = 1.6
        else:
            clim = 2.2

        # plots PWD matrix
        # uses KDE
        plotMatrix(
            SCmatrixCollated,
            uniqueBarcodes,
            pixelSize,
            numberROIs,
            outputFileName,
            logNameMD,
            figtitle="PWD matrix - KDE",
            mode="KDE",  # median or KDE
            clim=clim,
            cm="terrain",
            fileNameEnding="_PWDmatrixKDE.png",
        )  # need to validate use of KDE. For the moment it does not handle well null distributions

        # uses median
        plotMatrix(
            SCmatrixCollated,
            uniqueBarcodes,
            pixelSize,
            numberROIs,
            outputFileName,
            logNameMD,
            figtitle="PWD matrix - median",
            mode="median",  # median or KDE
            clim=clim,
            cm="coolwarm",
            fileNameEnding="_PWDmatrixMedian.png",
        )  # need to validate use of KDE. For the moment it does not handle well null distributions

        # calculates and plots contact probability matrix from merged samples/datasets
        HiMmatrix, nCells = calculateContactProbabilityMatrix(
            SCmatrixCollated, uniqueBarcodes, pixelSize, norm="nonNANs",
        )  # norm: nCells (default), nonNANs

        cScale = HiMmatrix.max()
        plotMatrix(
            HiMmatrix,
            uniqueBarcodes,
            pixelSize,
            numberROIs,
            outputFileName,
            logNameMD,
            figtitle="Hi-M matrix",
            mode="counts",
            clim=cScale,
            cm="coolwarm",
            fileNameEnding="_HiMmatrix.png",
        )

        # plots Nmatrix
        plotMatrix(
            Nmatrix,
            uniqueBarcodes,
            pixelSize,
            numberROIs,
            outputFileName,
            logNameMD,
            figtitle="N-matrix",
            mode="counts",
            clim=np.max(Nmatrix),
            cm="Blues",
            fileNameEnding="_Nmatrix.png",
        )

        plotDistanceHistograms(
            SCmatrixCollated, pixelSize, outputFileName, logNameMD, mode="KDE", kernelWidth=0.25, optimizeKernelWidth=False
        )

    def save_matrices(self):
        # saves output
        np.save(outputFileName + "_" + maskIdentifier + "_HiMscMatrix.npy", SCmatrixCollated)
        np.savetxt(outputFileName + "_" + maskIdentifier + "_uniqueBarcodes.ecsv", uniqueBarcodes, delimiter=" ", fmt="%d")
        np.save(outputFileName + "_" + maskIdentifier + "_Nmatrix.npy", Nmatrix)

    def launch_analysis(self, file):
        """
        run analysis for a chromatin trace table.

        Returns
        -------
        None.

        """
        # outputFileName + "_mask:" + str(self.maskIdentifier) + "_ROI:" + str(self.nROI) + ".ecsv"

        # creates and loads trace table
        self.trace_table = chromatin_trace_table()
        self.trace_table.load(file)

        # decodes ROI

        # runs calculation of PWD matrix
        self.buildsdistanceMatrix("min")  # mean min last

        # calculates N-matrix: number of PWD distances for each barcode combination
        self.calculatesNmatrix()

        # runs plotting operations
        self.plotsAllmatrices(file)

        # saves matrix
        self.save_matrices(file)

    def run(self):

        # initializes sessionName, dataFolder, currentFolder
        label = "barcode"
        self.dataFolder, self.currentFolder  = initialize_module(self.param, module_name="build_matrix",label = label)

        # reads chromatin traces
        files = [x for x in glob.glob(self.dataFolder.outputFiles["buildsPWDmatrix"] + "_mask*" + label + ".ecsv")]

        if len(files) < 1:
            printLog("$ No chromatin trace table found to process!","WARN")
            return

        for file in files:
            self.launch_analysis(file)

        printLog(f"$ {len(files)} chromatin trace tables processed in {self.currentFolder}")
