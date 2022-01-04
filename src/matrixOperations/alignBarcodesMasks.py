#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:23:36 2020

@author: marcnol

This script:
    - iterates over ROIs
        - assigns barcode localizations to masks
        - applies local drift correction, if available
        - removes localizations using flux and driftTolerance
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

from astropy.table import Table

from photutils.segmentation import SegmentationImage

from fileProcessing.fileManagement import (
    folders,
    writeString2File,
    printLog,
)

from matrixOperations.HIMmatrixOperations import plotMatrix, plotDistanceHistograms, calculateContactProbabilityMatrix

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")

# =============================================================================
# CLASSES
# =============================================================================


class cellID:
    def __init__(self, param, dataFolder, barcodeMapROI, Masks, ROI, ndims=2):
        self.param = param
        self.dataFolder = dataFolder
        self.barcodeMapROI = barcodeMapROI
        self.Masks = Masks
        self.NcellsAssigned = 0
        self.NcellsUnAssigned = 0
        self.NbarcodesinMask = 0
        self.ndims = ndims
        self.dictErrorBlockMasks = {}  # contains the results from blockAlignment, if existing

        self.SegmentationMask = SegmentationImage(self.Masks)
        self.numberMasks = self.SegmentationMask.nlabels
        self.ROI = ROI
        self.alignmentResultsTable = Table()
        self.alignmentResultsTableRead = False
        self.barcodesinMask = dict()
        self.logNameMD = ""
        self.foundMatch = []

        for mask in range(self.numberMasks + 1):
            self.barcodesinMask["maskID_" + str(mask)] = []

    def initializeLists(self):
        self.ROIs, self.cellID, self.nBarcodes, self.barcodeIDs, self.p, self.cuid, self.buid, self.barcodeCoordinates = (
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
        )

    def filterLocalizations_Quality(self, i, flux_min):
        """
        [filters barcode localizations either by brigthness or 3D localization accuracy]

        Parameters
        ----------
        i : int
            index in barcodeMap Table
        flux_min : float
            Minimum flux to keep barcode localization

        Returns
        -------
        keep : Boolean
            True if the test is passed.

        """
        if "3DfitKeep" in self.barcodeMapROI.groups[0].keys() and self.ndims == 3:
            # [reading the flag in barcodeMapROI assigned by the 3D localization routine]
            keep = self.barcodeMapROI.groups[0]["3DfitKeep"][i] and self.barcodeMapROI.groups[0]["flux"][i] > flux_min
        else:
            # [or by reading the flux from 2D localization]
            keep = self.barcodeMapROI.groups[0]["flux"][i] > flux_min

        return keep

    def filterLocalizations_BlockAlignment(self, i, toleranceDrift, blockSize):
        """
        [filters barcode per blockAlignmentMask, if existing]
        runs only if localAligment was not run!

        Parameters
        ----------
        i : int
            index in barcodeMap Table
        toleranceDrift : float
            tolerance to keep barcode localization, in pixel units
        blockSize : int
            size of blocks used for blockAlignment.

        Returns
        -------
        keepAlignment : Boolean
            True if the test is passed.

        """
        y_int = int(self.barcodeMapROI.groups[0]["xcentroid"][i])
        x_int = int(self.barcodeMapROI.groups[0]["ycentroid"][i])
        keepAlignment = True
        if not self.alignmentResultsTableRead: # only proceeds if localAlignment was not performed
            barcodeID = "barcode:" + str(self.barcodeMapROI.groups[0]["Barcode #"][i])
            barcodeROI = "ROI:" + str(self.barcodeMapROI.groups[0]["ROI #"][i])

            if len(self.dictErrorBlockMasks) > 0:
                if barcodeROI in self.dictErrorBlockMasks.keys():
                    if barcodeID in self.dictErrorBlockMasks[barcodeROI].keys():
                        errorMask = self.dictErrorBlockMasks[barcodeROI][barcodeID]
                        keepAlignment = (
                            errorMask[int(np.floor(x_int / blockSize)), int(np.floor(y_int / blockSize))]
                            < toleranceDrift
                        )

            # keeps it always if barcode is fiducial
            if (
                "RT" + str(self.barcodeMapROI.groups[0]["Barcode #"][i])
                in self.param.param["alignImages"]["referenceFiducial"]
            ):
                keepAlignment = True


        return keepAlignment

    def plots_distributionFluxes(self):
        fileName = (
            self.dataFolder.outputFolders["buildsPWDmatrix"]
            + os.sep
            + "BarcodeStats_ROI:"
            + str(self.nROI)
            + "_"
            + str(self.ndims)
            + "D.png"
        )

        fig, axes = plt.subplots(1, 2)
        ax = axes.ravel()
        fig.set_size_inches((10, 5))

        fluxes = self.barcodeMapROI.groups[0]["flux"]
        sharpness = self.barcodeMapROI.groups[0]["sharpness"]
        roundness = self.barcodeMapROI.groups[0]["roundness1"]
        peak = self.barcodeMapROI.groups[0]["peak"]
        mag = self.barcodeMapROI.groups[0]["mag"]

        # p1 = ax[0].hist(fluxes,bins=25)
        p1 = ax[0].scatter(fluxes, sharpness, c=peak, cmap="terrain", alpha=0.5)
        ax[0].set_title("color: peak intensity")
        ax[0].set_xlabel("flux")
        ax[0].set_ylabel("sharpness")

        p2 = ax[1].scatter(roundness, mag, c=peak, cmap="terrain", alpha=0.5)
        ax[1].set_title("color: peak intensity")
        ax[1].set_xlabel("roundness")
        ax[1].set_ylabel("magnitude")
        fig.colorbar(p2, ax=ax[1], fraction=0.046, pad=0.04)

        fig.savefig(fileName)

        plt.close(fig)

        writeString2File(
            self.logNameMD, "Barcode stats for ROI:{}, dims:{} \n![]({})\n".format(self.nROI, self.ndims, fileName), "a"
        )

    def plots_barcodesAlignment(self, blockSize):
        """
        plots barcode localizations together with the blockAlignment map

        Returns
        -------
        None.

        """
        fileName = (
            self.dataFolder.outputFolders["buildsPWDmatrix"]
            + os.sep
            + "BarcodeAlignmentAccuracy_ROI:"
            + str(self.nROI)
            + "_"
            + str(self.ndims)
            + "D.png"
        )

        fig, axes = plt.subplots()
        fig.set_size_inches((20, 20))

        accuracy, x, y = [], [], []
        printLog("> Plotting barcode alignments...")
        for i in trange(len(self.barcodeMapROI.groups[0])):
            barcodeID = "barcode:" + str(self.barcodeMapROI.groups[0]["Barcode #"][i])
            barcodeROI = "ROI:" + str(self.barcodeMapROI.groups[0]["ROI #"][i])
            y_int = int(self.barcodeMapROI.groups[0]["xcentroid"][i])
            x_int = int(self.barcodeMapROI.groups[0]["ycentroid"][i])

            if len(self.dictErrorBlockMasks) > 0:
                if barcodeROI in self.dictErrorBlockMasks.keys():
                    if barcodeID in self.dictErrorBlockMasks[barcodeROI].keys():
                        errorMask = self.dictErrorBlockMasks[barcodeROI][barcodeID]
                        accuracy.append(errorMask[int(np.floor(x_int / blockSize)), int(np.floor(y_int / blockSize))])
                        x.append(self.barcodeMapROI.groups[0]["xcentroid"][i])
                        y.append(self.barcodeMapROI.groups[0]["ycentroid"][i])

        p1 = axes.scatter(x, y, s=5, c=accuracy, cmap="terrain", alpha=0.5, vmin=0, vmax=5)
        fig.colorbar(p1, ax=axes, fraction=0.046, pad=0.04)
        axes.set_title("barcode drift correction accuracy, px")

        axes.axis("off")

        fig.savefig(fileName)

        plt.close(fig)

        writeString2File(
            self.logNameMD, "Barcode stats for ROI:{}, dims:{} \n![]({})\n".format(self.nROI, self.ndims, fileName), "a"
        )

    def alignByMasking(self):
        """
        Assigns barcodes to masks and creates <NbarcodesinMask>
        This routine will only select which barcodes go to each cell mask

        Returns
        -------
        self.barcodesinMask # dictionnary with the identities of barcodes contained in each mask.
            Keys: 'maskID_1', 'maskID_2', and so on

        self.NbarcodesinMask # vector containing the number of barcodes for each mask
        self.NcellsAssigned # number of cells assigned
        self.NcellsUnAssigned # number of cells unassigned
        """

        NbarcodesinMask = np.zeros(self.numberMasks + 2)
        flux_key = "flux_min_3D" if self.ndims==3 else "flux_min"
        if flux_key in self.param.param["buildsPWDmatrix"]:
            flux_min = self.param.param["buildsPWDmatrix"][flux_key]
        else:
            flux_min = 0
            printLog("# Flux min not found. Set to {}!".format(flux_min))

        if "toleranceDrift" in self.param.param["buildsPWDmatrix"]:
            toleranceDrift = self.param.param["buildsPWDmatrix"]["toleranceDrift"]
        else:
            toleranceDrift = 1
            printLog("# toleranceDrift not found. Set to {}!".format(toleranceDrift))

        if "blockSize" in self.param.param["alignImages"]:
            blockSize = self.param.param["alignImages"]["blockSize"]
        else:
            blockSize = 256
            printLog("# blockSize not found. Set to {}!".format(blockSize))

        printLog("\n$ ndims = {}\n$ Flux min = {} \n$ ToleranceDrift = {} px\n$ Reference barcode = {}".format(self.ndims,
                                                                                                    flux_min,
                                                                                                    toleranceDrift,
                                                                                                    self.param.param["alignImages"]["referenceFiducial"]))

        # Produces images of distribution of fluxes.
        self.plots_distributionFluxes()
        self.plots_barcodesAlignment(blockSize)

        keepQualityAll, keepAlignmentAll, NbarcodesROI = [], [], 0
        # loops over barcode Table rows in a given ROI
        printLog("> Aligning by masking...")
        for i in trange(len(self.barcodeMapROI.groups[0])): # i is the index of the barcode in barcodeMapROI
            barcode = self.barcodeMapROI.groups[0]["Barcode #"][i]
            ROI = self.barcodeMapROI.groups[0]["ROI #"][i]

            # [filters barcode localizations either by]
            keepQuality = self.filterLocalizations_Quality(i, flux_min)

            # [filters barcode per blockAlignmentMask, if existing]
            keepAlignment = self.filterLocalizations_BlockAlignment(i, toleranceDrift, blockSize)

            # applies all filters
            if keepQuality and keepAlignment:

                # keeps the particle if the test passed
                x_uncorrected = self.barcodeMapROI.groups[0]["ycentroid"][i] # control inversion between x-y
                y_uncorrected = self.barcodeMapROI.groups[0]["xcentroid"][i]

                if self.ndims==2:
                    z_uncorrected = self.barcodeMapROI.groups[0]["zcentroid"][i] = 0.0
                else:
                    z_uncorrected = self.barcodeMapROI.groups[0]["zcentroid"][i]

                y_int = int(y_uncorrected)
                x_int = int(x_uncorrected)

                # finds what mask label this barcode is sitting on
                maskID = self.Masks[x_int][y_int]

                # Corrects XYZ coordinate of barcode if localDriftCorrection is available
                zxy_uncorrected = [z_uncorrected, x_uncorrected , y_uncorrected]
                RTbarcode = "RT" + str(barcode)
                if  RTbarcode not in self.param.param["alignImages"]["referenceFiducial"]:
                    zxy_corrected = self.searchLocalShift(ROI, maskID, barcode, zxy_uncorrected,toleranceDrift)
                else:
                    # if it is the reference cycle, then it does not correct coordinates
                    zxy_corrected = zxy_uncorrected

                # rewrites corrected XYZ values to Table
                self.barcodeMapROI.groups[0]["ycentroid"][i] = zxy_corrected[1]
                self.barcodeMapROI.groups[0]["xcentroid"][i] = zxy_corrected[2]
                if self.ndims>2:
                    self.barcodeMapROI.groups[0]["zcentroid"][i] = zxy_corrected[0]

                # attributes CellID to a barcode
                self.barcodeMapROI["CellID #"][i] = maskID

                # if it is not background,
                if maskID > 0:
                    # increments counter of number of barcodes in the cell mask attributed
                    NbarcodesinMask[maskID] += 1

                    # stores the identify of the barcode to the mask
                    self.barcodesinMask["maskID_" + str(maskID)].append(i)

            # keeps statistics
            if int(self.barcodeMapROI.groups[0]["ROI #"][i]) == int(self.nROI):
                keepQualityAll.append(keepQuality)
                keepAlignmentAll.append(keepAlignment)
                NbarcodesROI += 1

        # Total number of masks assigned and not assigned
        self.NcellsAssigned = np.count_nonzero(NbarcodesinMask > 0)
        self.NcellsUnAssigned = self.numberMasks - self.NcellsAssigned

        # this list contains which barcodes are allocated to which masks
        self.NbarcodesinMask = NbarcodesinMask

        printLog("$ Number of localizations passing quality test: {} / {}".format(sum(keepQualityAll), NbarcodesROI))

        printLog("$ Number of localizations passing alignment test: {} / {}".format(sum(keepAlignmentAll), NbarcodesROI))

        printLog("$ Number of cells assigned: {} | discarded: {}".format(self.NcellsAssigned, self.NcellsUnAssigned))

    def searchLocalShift(self, ROI, CellID, barcode, zxy_uncorrected,toleranceDrift=1):

        if "mask2D" in self.param.param["alignImages"]["localAlignment"]:
            return self.searchLocalShift_mask2D(ROI, CellID, zxy_uncorrected)
        elif "block3D" in self.param.param["alignImages"]["localAlignment"] and self.alignmentResultsTableRead:
            return self.searchLocalShift_block3D(ROI, barcode, zxy_uncorrected,toleranceDrift)
        else: # no correction was applied because the localAlignmentTable was not found
            return zxy_uncorrected

    def searchLocalShift_block3D(self, ROI, barcode, zxy_uncorrected, toleranceDrift=1):
        """
        Searches for local drift for a specific barcode in a given ROI.
        If it exists then it adds to the uncorrected coordinates

        Parameters
        ----------
        ROI : string
            ROI used
        CellID: string
            ID of the cell
        x_uncorrected : float
            x coordinate.
        y_uncorrected : float
            y coordinate.

        Returns
        -------
        x_corrected : float
            corrected x coordinate.
        y_corrected : float
            corrected y coordinate.

        """
        _foundMatch = False

        # gets blockSize
        blockSizeXY = self.alignmentResultsTable[0]["blockXY"]

        # zxy coord in block reference coord system
        zxyBlock = [np.floor(a/blockSizeXY).astype(int) for a in zxy_uncorrected]

        for row in self.alignmentResultsTable:

            # I need to check that the XY coordinates from localization are the same as the ij indices from the block decomposition!

            if row["ROI #"] == ROI and row["label"] == "RT" + str(barcode) and row["block_i"] == zxyBlock[1] and row["block_j"] == zxyBlock[2]:
                _foundMatch = True
                shifts = [row["shift_z"],row["shift_x"],row["shift_y"]]

                # checks that drifts > toleranceDrift are not applied
                if max(shifts)<toleranceDrift:
                    zxy_corrected = [a+shift for a,shift in zip(zxy_uncorrected,shifts)]
                else:
                    zxy_corrected = zxy_uncorrected

                # check for quality of shift correction before applying it !!
                #!TODO

        # keeps uncorrected values if no match is found
        if not _foundMatch:
            printLog("# Did not find match for ROI #{} barcode #{}".format(ROI, barcode))
            zxy_corrected = zxy_uncorrected
            self.foundMatch.append(False)
        else:
            self.foundMatch.append(True)

        return zxy_corrected

    def searchLocalShift_mask2D(self, ROI, CellID, zxy_uncorrected):
        """
        Searches for local drift for current mask. If it exists then id adds it to the uncorrected coordinates

        Parameters
        ----------
        ROI : string
            ROI used
        CellID: string
            ID of the cell
        x_uncorrected : float
            x coordinate.
        y_uncorrected : float
            y coordinate.

        Returns
        -------
        x_corrected : float
            corrected x coordinate.
        y_corrected : float
            corrected y coordinate.

        """

        _foundMatch = False
        for row in self.alignmentResultsTable:
            if row["ROI #"] == ROI and row["CellID #"] == CellID:
                _foundMatch = True
                shifts=[0,row["shift_x"],row["shift_y"]]
                zxy_corrected = [a + shift for a, shift in zip(zxy_uncorrected,shifts)]

        # keeps uncorrected values if no match is found
        if not _foundMatch:
            printLog("# Did not find match for CellID #{} in ROI #{}".format(CellID, ROI))
            zxy_corrected = zxy_uncorrected
            self.foundMatch.append(False)
        else:
            self.foundMatch.append(True)

        return zxy_corrected

    def buildsVector(self, groupKeys, x, y, z):
        """
        Builds vector from coordinates

        Parameters
        ----------
        groupKeys : list
            list of headers in the barcodes table
        x : float
            x coordinates
        y : float
            y coordinates
        z : float
            z coordinates

        Returns
        -------
        R : np array
            vector with coordinates in nanometers.

        """

        R = np.column_stack((x*self.pixelSize["x"], y*self.pixelSize["y"], z*self.pixelSize["z"]))

        return R

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


    def buildsSCdistanceTable(self):
        """
        iterates over all masks, calculates PWD for each mask, assigns them to SCdistanceTable

        Returns
        -------
        SCdistanceTable

        """
        # sorts Table by cellID
        barcodeMapROI = self.barcodeMapROI
        barcodeMapROI_cellID = barcodeMapROI.group_by("CellID #")  # ROI data sorted by cellID

        self.initializeLists()

        # iterates over all cell masks in an ROI
        printLog("> Building SC distance Tables")
        for key, group in tzip(barcodeMapROI_cellID.groups.keys, barcodeMapROI_cellID.groups):
            if key["CellID #"] > 1:  # excludes cellID 0 as this is background

                groupKeys, CellID, ROI = group.keys(), key["CellID #"], group["ROI #"].data[0]

                # gets lists of x, y and z coordinates for barcodes assigned to a cell mask
                x, y, z = np.array(group["xcentroid"].data), np.array(group["ycentroid"].data), np.array(group["zcentroid"].data)

                # calculates the PWD between barcodes in CellID
                PWD = self.calculatesPWDsingleMask(ROI, CellID, groupKeys, x, y, z)
                R_nm = self.buildsVector(groupKeys, x, y, z)

                self.ROIs.append(group["ROI #"].data[0])
                self.cellID.append(key["CellID #"])
                self.nBarcodes.append(len(group))
                self.barcodeIDs.append(group["Barcode #"].data)
                self.buid.append(group["Buid"].data)
                self.p.append(PWD)
                self.barcodeCoordinates.append(R_nm)
                self.cuid.append(str(uuid.uuid4()))  # creates cell unique identifier

        printLog(
            "$ Local correction applied to {}/{} barcodes in ROI {}".format(
                np.nonzero(self.foundMatch)[0].shape[0], len(self.foundMatch), group["ROI #"].data[0]
            )
        )

        printLog("$ Coordinates dimensions: {}".format(self.ndims))

        SCdistanceTable = Table()
        SCdistanceTable["Cuid"] = self.cuid
        SCdistanceTable["ROI #"] = self.ROIs
        SCdistanceTable["CellID #"] = self.cellID
        SCdistanceTable["nBarcodes"] = self.nBarcodes
        SCdistanceTable["Barcode #"] = self.barcodeIDs
        SCdistanceTable["Buid"] = self.buid
        SCdistanceTable["PWDmatrix"] = self.p
        SCdistanceTable["barcode xyz, nm"] = self.barcodeCoordinates

        self.SCdistanceTable = SCdistanceTable

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
        # [ builds SCdistanceTable ]
        self.buildsSCdistanceTable()
        printLog("$ Cells with barcodes found: {}".format(len(self.SCdistanceTable)))

        # [ builds SCmatrix ]
        numberMatrices = len(self.SCdistanceTable)  # z dimensions of SCmatrix

        # uniqueBarcodes = np.unique(self.barcodeMapROI["Barcode #"].data)
        uniqueBarcodes = self.uniqueBarcodes

        # number of unique Barcodes for xy dimensions of SCmatrix
        numberUniqueBarcodes = uniqueBarcodes.shape[0]
        SCmatrix = np.zeros((numberUniqueBarcodes, numberUniqueBarcodes, numberMatrices))
        SCmatrix[:] = np.NaN

        # loops over cell masks
        for iCell, scPWDitem in zip(range(numberMatrices), self.SCdistanceTable):
            barcodes2Process = scPWDitem["Barcode #"]

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


# =============================================================================
# FUNCTIONS
# =============================================================================


def calculatesNmatrix(SCmatrix):

    numberCells = SCmatrix.shape[2]

    if numberCells > 0:
        Nmatrix = np.sum(~np.isnan(SCmatrix), axis=2)
    else:
        numberBarcodes = SCmatrix.shape[0]
        Nmatrix = np.zeros((numberBarcodes, numberBarcodes))

    return Nmatrix

def loadsLocalAlignment(param,dataFolder):

    if "None" in param.param["alignImages"]["localAlignment"]:
        printLog("\n\n$ localAlignment option set to {}".format(param.param["alignImages"]["localAlignment"]))
        return False, Table()
    else:
        return _loadsLocalAlignment(dataFolder,param.param["alignImages"]["localAlignment"])

def _loadsLocalAlignment(dataFolder,mode):

    localAlignmentFileName = dataFolder.outputFiles["alignImages"].split(".")[0] + "_" + mode + ".dat"
    if os.path.exists(localAlignmentFileName):
        alignmentResultsTable = Table.read(localAlignmentFileName, format="ascii.ecsv")
        alignmentResultsTableRead = True
        printLog("$ LocalAlignment file loaded: {}\n$ Will correct coordinates using {} alignment".format(localAlignmentFileName,mode))
        printLog("$ Number of records: {}".format(len(alignmentResultsTable)))
    else:
        printLog("\n\n# Warning: could not find localAlignment: {}\n Proceeding with only global alignments...".format(
                localAlignmentFileName
            )
        )
        alignmentResultsTableRead = False
        alignmentResultsTable = Table()

    return alignmentResultsTable, alignmentResultsTableRead

# def loadsLocalAlignment(dataFolder):
#     """
#     reads and returns localAlignmentTable, if it exists

#     Parameters
#     ----------
#     dataFolder : folder()
#         DESCRIPTION.

#     Returns
#     -------
#     alignmentResultsTable : Table()
#         DESCRIPTION.
#     alignmentResultsTableRead : Boolean
#         DESCRIPTION.

#     """
#     localAlignmentFileName = dataFolder.outputFiles["alignImages"].split(".")[0] + "_localAlignment.dat"
#     if os.path.exists(localAlignmentFileName):
#         alignmentResultsTable = Table.read(localAlignmentFileName, format="ascii.ecsv")
#         alignmentResultsTableRead = True
#         printLog("LocalAlignment file loaded !\nWill correct coordinates in XY")
#     else:
#         printLog(
#             "\n\n*** Warning: could not find localAlignment: {}\n Proceeding with only global alignments...".format(
#                 localAlignmentFileName
#             )
#         )
#         alignmentResultsTableRead = False
#         alignmentResultsTable = Table()

#     return alignmentResultsTable, alignmentResultsTableRead


def loadsBarcodeMap(fileNameBarcodeCoordinates, ndims):
    """
    Loads barcodeMap

    Parameters
    ----------
    fileNameBarcodeCoordinates : string
        filename with barcodeMap
    ndims : int
        either 2 or 3.

    Returns
    -------
    barcodeMap : Table()
    localizationDimension : int
        either 2 or 3.
    uniqueBarcodes: list
        lis of unique barcodes read from barcodeMap

    """
    if os.path.exists(fileNameBarcodeCoordinates):
        barcodeMap = Table.read(fileNameBarcodeCoordinates, format="ascii.ecsv")
        printLog("$ Successfully loaded barcode localizations file: {}".format(fileNameBarcodeCoordinates))

        uniqueBarcodes = np.unique(barcodeMap["Barcode #"].data)
        numberUniqueBarcodes = uniqueBarcodes.shape[0]

        printLog("Number Barcodes read from barcodeMap: {}".format(numberUniqueBarcodes))
        printLog("Unique Barcodes detected: {}".format(uniqueBarcodes))
    else:
        printLog("\n\n# ERROR: could not find coordinates file: {}".format(fileNameBarcodeCoordinates))
        sys.exit()

    return barcodeMap, ndims, uniqueBarcodes


def buildsDictionaryErrorAlignmentMasks(param, dataFolder):
    """
    Builds and returns dictionary with error alignment block masks produced during the alignment process if
    the 'blockAlignment' option was used

    Parameters
    ----------
    param : Parameters()
    dataFolder : folder()

    Returns
    -------
    dictErrorBlockMasks : dict

    """
    folder = dataFolder.outputFolders["alignImages"]
    fileList = glob.glob(folder + os.sep + "*_errorAlignmentBlockMap.npy")

    # decodes files and builds dictionnary
    fileNameRegExp = param.param["acquisition"]["fileNameRegExp"]
    fileNameRegExp = fileNameRegExp.split(".")[0]
    listRE = [re.search(fileNameRegExp, os.path.basename(x).split("_errorAlignmentBlockMap.npy")[0]) for x in fileList]

    dictErrorBlockMasks = dict()

    for file, regExp in zip(fileList, listRE):
        if "ROI:" + str(int(regExp["roi"])) not in dictErrorBlockMasks.keys():
            dictErrorBlockMasks["ROI:" + str(int(regExp["roi"]))] = {}
        if "barcode:" + regExp["cycle"].split("RT")[-1] not in dictErrorBlockMasks.keys():
            newMask = np.load(file)
            dictErrorBlockMasks["ROI:" + str(int(regExp["roi"]))][
                "barcode:" + regExp["cycle"].split("RT")[-1]
            ] = newMask

    return dictErrorBlockMasks


def plotsAllmatrices(
    SCmatrixCollated, Nmatrix, uniqueBarcodes, pixelSize, numberROIs, outputFileName, logNameMD, localizationDimension
):
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


def buildsPWDmatrix(
    param,
    currentFolder,
    fileNameBarcodeCoordinates,
    outputFileName,
    dataFolder,
    pixelSize={'x':0.1,'y':0.1,'z':0.0},
    logNameMD="log.md",
    ndims=2,
    maskIdentifier='DAPI'
):
    """
    Main function that:
        loads and processes barcode localization files, local alignment file, and masks
        initializes <cellROI> class and assigns barcode localizations to masks
        then constructs the single cell PWD matrix and outputs it toghether with the contact map and the N-map.

    Parameters
    ----------
    param : Parameters Class
    currentFolder : string
    fileNameBarcodeCoordinates : string
    outputFileName : string
    dataFolder : Folder Class
        information to find barcode localizations, local drift corrections and masks

    pixelSize : dict, optional
        pixelSize = {'x': pixelSizeXY,
                    'y': pixelSizeXY,
                    'z': pixelSizeZ}
        The default is 0.1 for x and y, 0.0 for z. Pixelsize in um

    logNameMD : str, optional
        Filename of Markdown output. The default is "log.md".
    ndims : int, optional
        indicates whether barcodes were localized in 2 or 3D. The default is 2.

    Returns
    -------
    None.

    """
    # Loads localAlignment if it exists
    alignmentResultsTable, alignmentResultsTableRead = loadsLocalAlignment(param,dataFolder)

    # Loads coordinate Tables
    barcodeMap, localizationDimension, uniqueBarcodes = loadsBarcodeMap(fileNameBarcodeCoordinates, ndims)

    # Builds dictionnary with filenames of errorAlignmentBlockMasks for each ROI and each barcode
    dictErrorBlockMasks = buildsDictionaryErrorAlignmentMasks(param, dataFolder)

    # processes tables
    barcodeMapROI = barcodeMap.group_by("ROI #")
    numberROIs = len(barcodeMapROI.groups.keys)
    printLog("\n$ ROIs detected: {}".format(numberROIs))

    # loops over ROIs
    filesinFolder = glob.glob(currentFolder + os.sep + "*.tif")
    # SCmatrixCollated, uniqueBarcodes, processingOrder = [], [], 0
    SCmatrixCollated, processingOrder = [], 0

    for ROI in range(numberROIs):
        nROI = barcodeMapROI.groups.keys[ROI][0]  # need to iterate over the first index

        printLog("----------------------------------------------------------------------")
        printLog("> Loading masks and pre-processing barcodes for Mask <{}> ROI# {}".format(maskIdentifier,nROI))
        printLog("----------------------------------------------------------------------")

        barcodeMapSingleROI = barcodeMap.group_by("ROI #").groups[ROI]

        # finds file with cell masks
        fileList2Process = [
            file
            for file in filesinFolder
            if file.split("_")[-1].split(".")[0] == param.param["acquisition"]["label_channel"] # typically "ch00" 
            and maskIdentifier in os.path.basename(file).split("_")
            and int(os.path.basename(file).split("_")[3]) == nROI
        ]

        if len(fileList2Process) > 0:

            # loads file with cell masks
            fileNameROImasks = os.path.basename(fileList2Process[0]).split(".")[0] + "_Masks.npy"
            fullFileNameROImasks = os.path.dirname(fileNameBarcodeCoordinates) + os.sep + fileNameROImasks
            if os.path.exists(fullFileNameROImasks):
                Masks = np.load(fullFileNameROImasks)

                # Assigns barcodes to Masks for a given ROI
                cellROI = cellID(param, dataFolder, barcodeMapSingleROI, Masks, ROI, ndims=localizationDimension)
                cellROI.ndims, cellROI.nROI, cellROI.logNameMD, cellROI.pixelSize = ndims, nROI, logNameMD, pixelSize
                cellROI.uniqueBarcodes = uniqueBarcodes

                if alignmentResultsTableRead:
                    cellROI.alignmentResultsTable = alignmentResultsTable

                cellROI.dictErrorBlockMasks = dictErrorBlockMasks
                cellROI.alignmentResultsTableRead = alignmentResultsTableRead

                # finds what barcodes are in each cell mask
                cellROI.alignByMasking()

                # builds the single cell distance Matrix
                cellROI.buildsdistanceMatrix("min")  # mean min last

                printLog(
                    "$ ROI: {}, N cells assigned: {} out of {}\n".format(
                        ROI, cellROI.NcellsAssigned - 1, cellROI.numberMasks
                    )
                )


                # saves Table with results per ROI
                cellROI.SCdistanceTable.write(
                    outputFileName + "_order:" + str(processingOrder) + "_ROI:" + str(nROI) + ".ecsv",
                    format="ascii.ecsv",
                    overwrite=True,
                )

                if len(SCmatrixCollated) > 0:
                    SCmatrixCollated = np.concatenate((SCmatrixCollated, cellROI.SCmatrix), axis=2)
                else:
                    SCmatrixCollated = cellROI.SCmatrix
                del cellROI

                processingOrder += 1

            # Could not find a file with masks to assign. Report and continue with next ROI
            ###############################################################################
            else:
                printLog(
                    "# Error, no mask file found for ROI: {}, segmentedMasks: {}\n".format(
                        nROI, fileNameBarcodeCoordinates
                    )
                )
                printLog("# File I was searching for: {}".format(fullFileNameROImasks))
                printLog("# Debug: ")
                for file in filesinFolder:
                    if (
                        file.split("_")[-1].split(".")[0] == param.param["acquisition"]["label_channel"] # typically "ch00" 
                        and maskIdentifier in file.split("_")
                        and int(os.path.basename(file).split("_")[3]) == nROI
                    ):
                        printLog("$ Hit found!")
                    printLog(
                        "fileSplit:{}, {} in filename: {}, ROI: {}".format(
                            file.split("_")[-1].split(".")[0],
                            maskIdentifier,
                            maskIdentifier in os.path.basename(file).split("_"),
                            int(os.path.basename(file).split("_")[3]),
                        )
                    )

    if processingOrder > 0:
        # calculates N-matrix: number of PWD distances for each barcode combination
        Nmatrix = calculatesNmatrix(SCmatrixCollated)

        # saves output
        np.save(outputFileName + "_" + maskIdentifier + "_HiMscMatrix.npy", SCmatrixCollated)
        np.savetxt(outputFileName + "_" + maskIdentifier + "_uniqueBarcodes.ecsv", uniqueBarcodes, delimiter=" ", fmt="%d")
        np.save(outputFileName + "_" + maskIdentifier + "_Nmatrix.npy", Nmatrix)
        pixelSizeXY = pixelSize['x']

        if SCmatrixCollated.shape[2]>0:
            #################################
            # makes and saves outputs plots #
            #################################
            plotsAllmatrices(
                SCmatrixCollated,
                Nmatrix,
                uniqueBarcodes,
                pixelSizeXY,
                numberROIs,
                outputFileName + "_" + maskIdentifier,
                logNameMD,
                localizationDimension,
            )
        else:
            printLog("# Nothing to plot. Single cell matrix is empty. Number of cells: {}".format(SCmatrixCollated.shape[2]))

def processesPWDmatrices(param, session1):
    """
    Function that assigns barcode localizations to masks and constructs single cell cummulative PWD matrix.

    Parameters
    ----------
    param : class
        Parameters
    log1 : class
        logging class.
    session1 : class
        session information

    Returns
    -------
    None.

    """
    sessionName = "buildsPWDmatrix"

    # processes folders and files
    dataFolder = folders(param.param["rootFolder"])
    printLog("\n===================={}====================\n".format(sessionName))
    printLog("$ folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(param.param["fileNameMD"], "## {}\n".format(sessionName), "a")
    label = 'barcode'

    for currentFolder in dataFolder.listFolders:
        # filesFolder=glob.glob(currentFolder+os.sep+'*.tif')
        dataFolder.createsFolders(currentFolder, param)
        printLog("> Processing Folder: {}".format(currentFolder))

        availableMasks = param.param["buildsPWDmatrix"]["masks2process"]
        printLog("> Masks labels: {}".format(availableMasks))

        for maskLabel in availableMasks.keys():
            
            maskIdentifier = availableMasks[maskLabel]
            
            fileNameBarcodeCoordinates = dataFolder.outputFiles["segmentedObjects"] + "_" + label + ".dat"
            if os.path.exists(fileNameBarcodeCoordinates):
                # 2D
                outputFileName = dataFolder.outputFiles["buildsPWDmatrix"]
                printLog("> 2D processing: {}".format(outputFileName))
    
                if "pixelSizeXY" in param.param["acquisition"].keys():
                    pixelSizeXY = param.param["acquisition"]["pixelSizeXY"]
                    pixelSize = {'x': pixelSizeXY,
                                 'y': pixelSizeXY,
                                 'z': 0.0}
                else:
                    pixelSize = {'x': 0.1,
                                 'y': 0.1,
                                 'z': 0.0}
    
    
                buildsPWDmatrix(
                    param,
                    currentFolder,
                    fileNameBarcodeCoordinates,
                    outputFileName,
                    dataFolder,
                    pixelSize,
                    param.param["fileNameMD"],
                    maskIdentifier =maskIdentifier, 
                )
    
            # 3D
            fileNameBarcodeCoordinates = dataFolder.outputFiles["segmentedObjects"] + "_3D_" + label + ".dat"
            if os.path.exists(fileNameBarcodeCoordinates):
                outputFileName = dataFolder.outputFiles["buildsPWDmatrix"] + "_3D"
                printLog("> 3D processing: {}".format(outputFileName))
    
                if ("pixelSizeZ" in param.param["acquisition"].keys()) and ("pixelSizeXY" in param.param["acquisition"].keys()):
                    pixelSizeXY = param.param["acquisition"]["pixelSizeXY"]
    
                    if 'zBinning' in param.param['acquisition']:
                        zBinning = param.param['acquisition']['zBinning']
                    else:
                        zBinning = 1
    
                    pixelSizeZ = zBinning*param.param["acquisition"]["pixelSizeZ"]
    
                    pixelSize = {'x': pixelSizeXY,
                                 'y': pixelSizeXY,
                                 'z': pixelSizeZ*zBinning}
                else:
                    pixelSize = {'x': 0.1,
                                 'y': 0.1,
                                 'z': 0.25}
    
                buildsPWDmatrix(
                    param,
                    currentFolder,
                    fileNameBarcodeCoordinates,
                    outputFileName,
                    dataFolder,
                    pixelSize,
                    param.param["fileNameMD"],
                    ndims=3,
                    maskIdentifier =maskIdentifier,                     
                )
    
            # tights loose ends
            session1.add(currentFolder, sessionName)
    
            printLog("HiM matrix in {} processed".format(currentFolder), "info")
