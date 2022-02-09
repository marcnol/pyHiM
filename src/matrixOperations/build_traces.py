#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 14:11:58 2022

@author: marcnol

This script will build chromatin traces using a segmentObjects_barcode table 

The methods that will be implemented are:
    1= assigment by mask (either DAPI mask or other)
    2= spatial clusterization using KDtree. This method is mask-free.
    
    

Method 1:
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


class build_traces:
    def __init__(self, param, dataFolder, barcodeMapROI, Masks, ROI, ndims=2):
        self.param = param
        self.dataFolder = dataFolder
        self.barcodeMapROI = barcodeMapROI
        self.Masks = Masks
        self.NcellsAssigned = 0
        self.NcellsUnAssigned = 0
        self.NbarcodesinMask = 0
        self.ndims = ndims
#        self.SegmentationMask = SegmentationImage(self.Masks)
        self.numberMasks = self.SegmentationMask.nlabels
        self.ROI = ROI
        self.barcodesinMask = dict()

        for mask in range(self.numberMasks + 1):
            self.barcodesinMask["maskID_" + str(mask)] = []

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
        if not self.alignmentResultsTableRead:  # only proceeds if localAlignment was not performed
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
        NbarcodesROI = 0
        
        # loops over barcode Table rows in a given ROI
        printLog("> Aligning by masking...")
        for i in trange(len(self.barcodeMapROI.groups[0])):  # i is the index of the barcode in barcodeMapROI
            barcode = self.barcodeMapROI.groups[0]["Barcode #"][i]

            # keeps the particle if the test passed
            x_corrected = self.barcodeMapROI.groups[0]["ycentroid"][i]
            y_corrected = self.barcodeMapROI.groups[0]["xcentroid"][i]
            y_int = int(y_corrected)
            x_int = int(x_corrected)
            
            if self.ndims == 2:
                z_corrected = self.barcodeMapROI.groups[0]["zcentroid"][i] = 0.0
            else:
                z_corrected = self.barcodeMapROI.groups[0]["zcentroid"][i]

            # finds what mask label this barcode is sitting on
            maskID = self.Masks[x_int][y_int]

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
                NbarcodesROI += 1

        # Total number of masks assigned and not assigned
        self.NcellsAssigned = np.count_nonzero(NbarcodesinMask > 0)
        self.NcellsUnAssigned = self.numberMasks - self.NcellsAssigned

        # this list contains which barcodes are allocated to which masks
        self.NbarcodesinMask = NbarcodesinMask

        printLog("$ Number of cells assigned: {} | discarded: {}".format(self.NcellsAssigned, self.NcellsUnAssigned))

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

        R = np.column_stack((x * self.pixelSize["x"], y * self.pixelSize["y"], z * self.pixelSize["z"]))

        return R

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
        printLog("> Building SC traces")
        for key, group in tzip(barcodeMapROI_cellID.groups.keys, barcodeMapROI_cellID.groups):
            if key["CellID #"] > 1:  # excludes cellID 0 as this is background

                groupKeys, CellID, ROI = group.keys(), key["CellID #"], group["ROI #"].data[0]

                # gets lists of x, y and z coordinates for barcodes assigned to a cell mask
                x, y, z = (
                    np.array(group["xcentroid"].data),
                    np.array(group["ycentroid"].data),
                    np.array(group["zcentroid"].data),
                )

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
    
# =============================================================================
# FUNCTIONS
# =============================================================================

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

def buildsPWDmatrix(
    param,
    currentFolder,
    fileNameBarcodeCoordinates,
    outputFileName,
    dataFolder,
    pixelSize={"x": 0.1, "y": 0.1, "z": 0.0},
    logNameMD="log.md",
    ndims=2,
    maskIdentifier="DAPI",
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
    # Loads coordinate Tables
    barcodeMap, localizationDimension, uniqueBarcodes = loadsBarcodeMap(fileNameBarcodeCoordinates, ndims)

    # processes tables
    barcodeMapROI = barcodeMap.group_by("ROI #")
    numberROIs = len(barcodeMapROI.groups.keys)
    printLog("\n$ ROIs detected: {}".format(numberROIs))

    # loops over ROIs
    filesinFolder = glob.glob(currentFolder + os.sep + "*.tif")
    SCmatrixCollated, processingOrder = [], 0

    for ROI in range(numberROIs):
        nROI = barcodeMapROI.groups.keys[ROI][0]  # need to iterate over the first index

        printLog("----------------------------------------------------------------------")
        printLog("> Loading masks and pre-processing barcodes for Mask <{}> ROI# {}".format(maskIdentifier, nROI))
        printLog("----------------------------------------------------------------------")

        barcodeMapSingleROI = barcodeMap.group_by("ROI #").groups[ROI]

        # finds file with cell masks
        fileList2Process = [
            file
            for file in filesinFolder
            if file.split("_")[-1].split(".")[0] == param.param["acquisition"]["label_channel"]  # typically "ch00"
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
                        file.split("_")[-1].split(".")[0]
                        == param.param["acquisition"]["label_channel"]  # typically "ch00"
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
    label = "barcode"

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
                    pixelSize = {"x": pixelSizeXY, "y": pixelSizeXY, "z": 0.0}
                else:
                    pixelSize = {"x": 0.1, "y": 0.1, "z": 0.0}

                buildsPWDmatrix(
                    param,
                    currentFolder,
                    fileNameBarcodeCoordinates,
                    outputFileName,
                    dataFolder,
                    pixelSize,
                    param.param["fileNameMD"],
                    maskIdentifier=maskIdentifier,
                )

            # 3D
            fileNameBarcodeCoordinates = dataFolder.outputFiles["segmentedObjects"] + "_3D_" + label + ".dat"
            if os.path.exists(fileNameBarcodeCoordinates):
                outputFileName = dataFolder.outputFiles["buildsPWDmatrix"] + "_3D"
                printLog("> 3D processing: {}".format(outputFileName))

                if ("pixelSizeZ" in param.param["acquisition"].keys()) and (
                    "pixelSizeXY" in param.param["acquisition"].keys()
                ):
                    pixelSizeXY = param.param["acquisition"]["pixelSizeXY"]

                    if "zBinning" in param.param["acquisition"]:
                        zBinning = param.param["acquisition"]["zBinning"]
                    else:
                        zBinning = 1

                    pixelSizeZ = zBinning * param.param["acquisition"]["pixelSizeZ"]

                    pixelSize = {"x": pixelSizeXY, "y": pixelSizeXY, "z": pixelSizeZ * zBinning}
                else:
                    pixelSize = {"x": 0.1, "y": 0.1, "z": 0.25}

                buildsPWDmatrix(
                    param,
                    currentFolder,
                    fileNameBarcodeCoordinates,
                    outputFileName,
                    dataFolder,
                    pixelSize,
                    param.param["fileNameMD"],
                    ndims=3,
                    maskIdentifier=maskIdentifier,
                )

            # tights loose ends
            session1.add(currentFolder, sessionName)

            printLog("HiM matrix in {} processed".format(currentFolder), "info")
