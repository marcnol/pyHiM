#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:23:36 2020

@author: marcnol

test fitting barcode spots to masks


TO SOLVE:
    - I need to find a simple way of specifying the genomic coordinates of 
    barcodes for the production of the Hi-M matrix. At the moment it is just 
    using barcodeID from the file name.
    

"""

# =============================================================================
# IMPORTS
# =============================================================================

import glob, os
import uuid
import numpy as np
import matplotlib.pyplot as plt
from tqdm.contrib import tzip
from numba import jit

from sklearn.metrics import pairwise_distances
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import LeaveOneOut
from sklearn.neighbors import KernelDensity

from astropy.table import Table

from photutils.segmentation import SegmentationImage

from fileProcessing.fileManagement import (
    folders,
    isnotebook,
    session,
    writeString2File,
    Parameters, 
    log)

# to remove in a future version
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# CLASSES
# =============================================================================


class cellID:
    def __init__(self, param, barcodeMapROI, Masks, ROI,ndims=2):
        self.param=param
        self.barcodeMapROI = barcodeMapROI
        self.Masks = Masks
        self.NcellsAssigned = 0
        self.NcellsUnAssigned = 0
        self.NbarcodesinMask = 0
        self.ndims=ndims
        
        self.SegmentationMask = SegmentationImage(self.Masks)
        self.numberMasks = self.SegmentationMask.nlabels
        self.ROI = ROI
        self.alignmentResultsTable=Table()
        self.barcodesinMask = dict()
        for mask in range(self.numberMasks + 1):
            self.barcodesinMask["maskID_" + str(mask)] = []

    def initializeLists(self):
        self.ROIs, self.cellID, self.nBarcodes, self.barcodeIDs, self.p, self.cuid, self.buid = [], [], [], [], [], [], []

    # def visualize(self):
    #     pass
    #     # imageBarcodes = np.zeros([2048, 2048])
    #     # MasksBarcodes = Masks
    #     # R = []

    #     # for i in range(len(self.barcodeMapROI.groups[0])):
    #     #     y_int = int(self.barcodeMapROI.groups[0]["xcentroid"][i])
    #     #     x_int = int(self.barcodeMapROI.groups[0]["ycentroid"][i])
    #     #     barcodeID = self.barcodeMapROI.groups[0]["Barcode #"][i]
    #     #     imageBarcodes[x_int][y_int] = barcodeID
    #     #     MasksBarcodes[x_int][y_int] += 20 * barcodeID
    #     #     R.append([y_int, x_int, barcodeID])

    #     # # Shows results
    #     # Ra = np.array(R)
    #     # # plt.imshow(Masks, origin="lower", cmap="jet")
    #     # plt.scatter(Ra[:, 0], Ra[:, 1], s=5, c=Ra[:, 2], alpha=0.5)

    def alignByMasking(self):
        '''
        Assigns barcodes to masks and creates <NbarcodesinMask>
        '''

        NbarcodesinMask = np.zeros(self.numberMasks + 2)
        if "flux_min" in self.param.param["segmentedObjects"]:
            flux_min = self.param.param["segmentedObjects"]["flux_min"]
        else:
            flux_min = 0
        # print("Flux min = {}".format(flux_min))
        
        # loops over barcode Table rows in a given ROI
        for i in range(len(self.barcodeMapROI.groups[0])):
            
            if "3DfitKeep" in self.barcodeMapROI.groups[0].keys() and self.ndims==3:
                keep = self.barcodeMapROI.groups[0]["3DfitKeep"][i]
            else:
                keep = self.barcodeMapROI.groups[0]["flux"][i] > flux_min 
            
            if keep:
                y_int = int(self.barcodeMapROI.groups[0]["xcentroid"][i])
                x_int = int(self.barcodeMapROI.groups[0]["ycentroid"][i])
    
                #finds what mask label this barcode is sitting on
                maskID = self.Masks[x_int][y_int]
    
                # attributes CellID to a barcode            
                self.barcodeMapROI["CellID #"][i] = maskID
                
                # if it is not background,
                if maskID > 0:
                    # increments counter of number of barcodes in the cell mask attributed
                    NbarcodesinMask[maskID] += 1
                    
                    # stores the identify of the barcode to the mask
                    self.barcodesinMask["maskID_" + str(maskID)].append(i)

        # Total number of masks assigned and not assigned
        self.NcellsAssigned = np.count_nonzero(NbarcodesinMask > 0)
        self.NcellsUnAssigned = self.numberMasks - self.NcellsAssigned

        # this list contains which barcodes are allocated to which masks
        self.NbarcodesinMask = NbarcodesinMask
        
    def searchLocalShift(self,ROI,CellID,x_uncorrected,y_uncorrected):
        '''
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

        '''
        _foundMatch=False
        x_corrected, y_corrected = [], []
        for row in self.alignmentResultsTable:
            if row["ROI #"]==ROI and row["CellID #"]==CellID:
                _foundMatch=True
                x_corrected, y_corrected = x_uncorrected + row["shift_x"], y_uncorrected + row["shift_y"]
                
        # keeps uncorrected values if no match is found                            
        if not _foundMatch:
            print("Did not find match for CellID #{} in ROI #{}".format(CellID,ROI))
            x_corrected,y_corrected = x_uncorrected,y_uncorrected
            self.foundMatch.append(False)                    
        else:
            self.foundMatch.append(True)
            
        return x_corrected, y_corrected

    def buildsVector(self,groupKeys,x,y,z):
        '''
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
            vector with coordinates.

        '''
        if self.ndims==3 and "zcentroidGauss" in groupKeys:
            R = np.column_stack((x,y,z,)) 
        else:
            R = np.column_stack((x,y,)) 

        return R

    def calculatesPWDsingleMask(self,ROI,CellID,groupKeys,x_uncorrected, y_uncorrected,z_uncorrected):
        '''
        Calculates PWD between barcodes detected in a given mask

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
        Returns pairwise distance matrix between corrected barcodes

        '''

        if len(self.alignmentResultsTable)>0:                  
 
            # searches for local alignment shift for this mask in this ROI
            x_corrected, y_corrected = self.searchLocalShift(ROI,CellID,x_uncorrected,y_uncorrected)
              
            # applies local drift correction
            R = self.buildsVector(groupKeys,x_corrected, y_corrected,z_uncorrected )
                    
        else:
            # does not apply local drift correction
            R = self.buildsVector(groupKeys, x_uncorrected, y_uncorrected,z_uncorrected )

        return pairwise_distances(R)


    def buildsSCdistanceTable(self):
        '''
        iterates over all masks, calculates PWD for each mask, assigns them to SCdistanceTable

        Returns
        -------
        SCdistanceTable

        '''
        # sorts Table by cellID
        barcodeMapROI = self.barcodeMapROI
        barcodeMapROI_cellID = barcodeMapROI.group_by("CellID #")  # ROI data sorted by cellID
        
        self.initializeLists()
        
        # iterates over all cell masks in an ROI
        self.foundMatch=[]
        for key, group in tzip(barcodeMapROI_cellID.groups.keys, barcodeMapROI_cellID.groups):
            if key["CellID #"] > 1:  # excludes cellID 0 as this is background

                x_uncorrected, y_uncorrected= np.array(group["xcentroid"].data), np.array(group["ycentroid"].data)
                groupKeys, CellID, ROI = group.keys(), key["CellID #"], group["ROI #"].data[0]
                
                if self.ndims==3 and "zcentroidGauss" in group.keys():
                    z_uncorrected = np.array(group["zcentroidGauss"].data)
                else:
                    z_uncorrected = []

                PWD=self.calculatesPWDsingleMask(ROI,CellID,groupKeys,x_uncorrected, y_uncorrected,z_uncorrected)

                self.ROIs.append(group["ROI #"].data[0])
                self.cellID.append(key["CellID #"])
                self.nBarcodes.append(len(group))
                self.barcodeIDs.append(group["Barcode #"].data)
                self.buid.append(group["Buid"].data)
                self.p.append(PWD)
                self.cuid.append(str(uuid.uuid4()))  # creates cell unique identifier
                # print("CellID #={}, nBarcodes={}".format(key['CellID #'],len(group)))
                
        print("Local correction applied to {}/{} barcodes in ROI {}".format(np.nonzero(self.foundMatch)[0].shape[0],len(self.foundMatch),group["ROI #"].data[0]))
        if self.ndims==3 and "zcentroidGauss" in group.keys():
            print("Coordinates dimensions: 3")
        else:
            print("Coordinates dimensions: 2")
            
        SCdistanceTable = Table()
        SCdistanceTable["Cuid"] = self.cuid
        SCdistanceTable["ROI #"] = self.ROIs
        SCdistanceTable["CellID #"] = self.cellID
        SCdistanceTable["nBarcodes"] = self.nBarcodes
        SCdistanceTable["Barcode #"] = self.barcodeIDs
        SCdistanceTable["Buid"] = self.buid
        SCdistanceTable["PWDmatrix"] = self.p

        self.SCdistanceTable=SCdistanceTable
    
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
        # print("building distance matrix")

        # [ builds SCdistanceTable ]
        self.buildsSCdistanceTable()
        print("Cells with barcodes found: {}".format(len(self.SCdistanceTable)))

        # [ builds SCmatrix ]

        numberMatrices = len(self.SCdistanceTable)  # z dimensions of SCmatrix
        uniqueBarcodes = np.unique(self.barcodeMapROI["Barcode #"].data)
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
        self.uniqueBarcodes = uniqueBarcodes

# =============================================================================
# FUNCTIONS
# =============================================================================

@jit(nopython=True)
def findsOptimalKernelWidth(distanceDistribution):
    bandwidths = 10 ** np.linspace(-1, 1, 100)
    grid = GridSearchCV(KernelDensity(kernel="gaussian"), {"bandwidth": bandwidths}, cv=LeaveOneOut())
    grid.fit(distanceDistribution[:, None])
    return grid.best_params_

@jit(nopython=True)
def retrieveKernelDensityEstimator(distanceDistribution0, x_d, optimizeKernelWidth=False):
    '''
    Gets the kernel density function and maximum from a distribution of PWD distances

    Parameters
    ----------
    distanceDistribution0 : nd array
        List of PWD distances.
    x_d : nd array
        x grid.
    optimizeKernelWidth : Boolean, optional
        whether to optimize bandwidth. The default is False.

    Returns
    -------
    np array
        KDE distribution.
    np array
        Original distribution without NaNs

    '''

    nan_array = np.isnan(distanceDistribution0)

    not_nan_array = ~nan_array

    distanceDistribution = distanceDistribution0[not_nan_array]

    # instantiate and fit the KDE model
    if optimizeKernelWidth:
        kernelWidth = findsOptimalKernelWidth(distanceDistribution)["bandwidth"]
    else:
        kernelWidth = 0.3

    kde = KernelDensity(bandwidth=kernelWidth, kernel="gaussian")
    
    # makes sure the list is not full of NaNs.
    if distanceDistribution.shape[0]>0:
        kde.fit(distanceDistribution[:, None])
    else:
        return np.array([0]), np.array([0])
    
    # score_samples returns the log of the probability density
    logprob = kde.score_samples(x_d[:, None])

    return logprob, distanceDistribution


@jit(nopython=True)
def distributionMaximumKernelDensityEstimation(SCmatrixCollated, bin1, bin2, pixelSize, optimizeKernelWidth=False):
    '''
    calculates the kernel distribution and its maximum from a set of PWD distances

    Parameters
    ----------
    SCmatrixCollated : np array 3 dims
        SC PWD matrix.
    bin1 : int
        first bin.
    bin2 : int
        first bin.
    pixelSize : float
        pixel size in um
    optimizeKernelWidth : Boolean, optional
        does kernel need optimization?. The default is False.

    Returns
    -------
    float
        maximum of kernel.
    np array
        list of PWD distances used.
    np array
        kernel distribution.
    x_d : np array
        x grid.

    '''
    distanceDistribution0 = pixelSize * SCmatrixCollated[bin1, bin2, :]
    x_d = np.linspace(0, 5, 2000)

    # checks that distribution is not empty
    if distanceDistribution0.shape[0] > 0:
        logprob, distanceDistribution = retrieveKernelDensityEstimator(distanceDistribution0, x_d, optimizeKernelWidth)
        if logprob.shape[0]>1:
            kernelDistribution = 10 * np.exp(logprob)
            maximumKernelDistribution = x_d[np.argmax(kernelDistribution)]
            return maximumKernelDistribution, distanceDistribution, kernelDistribution, x_d
        else:
            return np.nan, np.zeros(x_d.shape[0]), np.zeros(x_d.shape[0]), x_d
    else:
        return np.nan, np.zeros(x_d.shape[0]), np.zeros(x_d.shape[0]), x_d


def plotMatrix(
    SCmatrixCollated,
    uniqueBarcodes,
    pixelSize,
    numberROIs=1,
    outputFileName="test",
    logNameMD="log.md",
    clim=1.4,
    cm="seismic",
    figtitle="PWD matrix",
    cmtitle="distance, um",
    nCells=0,
    mode="median",
    inverseMatrix=False,
    cMin=0,
    cells2Plot=[],
):
    Nbarcodes = SCmatrixCollated.shape[0]
    # projects matrix by calculating median in the nCell direction

    # Calculates ensemble matrix from single cell matrices
    if len(SCmatrixCollated.shape) == 3:
        if len(cells2Plot) == 0:
            cells2Plot = range(SCmatrixCollated.shape[2])

        if mode == "median":
            # calculates the median of all values
            if max(cells2Plot) > SCmatrixCollated.shape[2]:
                print(
                    "Error with range in cells2plot {} as it is larger than the number of available cells {}".format(
                        max(cells2Plot), SCmatrixCollated.shape[2]
                    )
                )
                keepPlotting = False
            else:
                meanSCmatrix = pixelSize * np.nanmedian(SCmatrixCollated[:, :, cells2Plot], axis=2)
                nCells = SCmatrixCollated[:, :, cells2Plot].shape[2]
                # print("Dataset {} cells2plot: {}".format(figtitle, nCells))
                # print('nCells={}'.format(nCells))
                keepPlotting = True
        elif mode == "KDE":
            # performs a Kernel Estimation to calculate the max of the distribution
            keepPlotting = True
            meanSCmatrix = np.zeros((Nbarcodes, Nbarcodes))
            for bin1 in range(Nbarcodes):
                for bin2 in range(Nbarcodes):
                    if bin1 != bin2:
                        (maximumKernelDistribution, _, _, _,) = distributionMaximumKernelDensityEstimation(
                            SCmatrixCollated[:, :, cells2Plot], bin1, bin2, pixelSize, optimizeKernelWidth=False,
                        )
                        meanSCmatrix[bin1, bin2] = maximumKernelDistribution
    else:
        if mode == "counts":
            meanSCmatrix = SCmatrixCollated
            keepPlotting = True
        else:
            meanSCmatrix = pixelSize * SCmatrixCollated
            keepPlotting = True

    if keepPlotting:
        # Calculates the inverse distance matrix if requested in the argument.
        if inverseMatrix:
            meanSCmatrix = np.reciprocal(meanSCmatrix)

        # plots figure
        fig = plt.figure(figsize=(10, 10))
        pos = plt.imshow(meanSCmatrix, cmap=cm)  # colormaps RdBu seismic
        plt.xlabel("barcode #")
        plt.ylabel("barcode #")
        plt.title(
            figtitle
            + " | "
            + str(meanSCmatrix.shape[0])
            + " barcodes | n="
            + str(nCells)
            + " | ROIs="
            + str(numberROIs)
        )
        plt.xticks(np.arange(SCmatrixCollated.shape[0]), uniqueBarcodes)
        plt.yticks(np.arange(SCmatrixCollated.shape[0]), uniqueBarcodes)
        cbar = plt.colorbar(pos, fraction=0.046, pad=0.04)
        cbar.minorticks_on()
        cbar.set_label(cmtitle)
        plt.clim(cMin, clim)

        if len(outputFileName.split(".")) > 1:
            if outputFileName.split(".")[1] != "png":
                plt.savefig(outputFileName)
            else:
                plt.savefig(outputFileName.split(".")[0] + "_HiMmatrix.png")
        else:
            plt.savefig(outputFileName + "_HiMmatrix.png")

        if not isnotebook():
            plt.close()

        writeString2File(logNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")
    else:
        print("Error plotting figure. Not executing script to avoid crash.")


def plotDistanceHistograms(
    SCmatrixCollated, pixelSize, outputFileName="test", logNameMD="log.md", mode="hist", limitNplots=10,
):

    if not isnotebook():
        NplotsX = NplotsY = SCmatrixCollated.shape[0]
    else:
        if limitNplots == 0:
            NplotsX = NplotsY = SCmatrixCollated.shape[0]
        else:
            NplotsX = NplotsY = min(
                [limitNplots, SCmatrixCollated.shape[0]]
            )  # sets a max of subplots if you are outputing to screen!

    bins = np.arange(0, 4, 0.25)

    sizeX, sizeY = NplotsX * 4, NplotsY * 4

    fig, axs = plt.subplots(figsize=(sizeX, sizeY), ncols=NplotsX, nrows=NplotsY, sharex=True)

    for i in range(NplotsX):
        for j in range(NplotsY):
            if i != j:
                # print('Printing [{}:{}]'.format(i,j))
                if mode == "hist":
                    axs[i, j].hist(pixelSize * SCmatrixCollated[i, j, :], bins=bins)
                else:
                    (maxKDE, distanceDistribution, KDE, x_d,) = distributionMaximumKernelDensityEstimation(
                        SCmatrixCollated, i, j, pixelSize, optimizeKernelWidth=False
                    )
                    axs[i, j].fill_between(x_d, KDE, alpha=0.5)
                    axs[i, j].plot(
                        distanceDistribution, np.full_like(distanceDistribution, -0.01), "|k", markeredgewidth=1,
                    )
                    axs[i, j].vlines(maxKDE, 0, KDE.max(), colors="r")

    plt.xlabel("distance, um")
    plt.ylabel("counts")
    plt.savefig(outputFileName + "_PWDhistograms.png")

    if not isnotebook():
        plt.close()

    writeString2File(logNameMD, "![]({})\n".format(outputFileName + "_PWDhistograms.png"), "a")

@jit(nopython=True)
def calculateContactProbabilityMatrix(iSCmatrixCollated, iuniqueBarcodes, pixelSize, threshold=0.25, norm="nCells"):
    # SCthresholdMatrix=iSCmatrixCollated<threshold

    nX = nY = iSCmatrixCollated.shape[0]
    nCells = iSCmatrixCollated.shape[2]
    SCmatrix = np.zeros((nX, nY))

    for i in range(nX):
        for j in range(nY):
            if i != j:
                distanceDistribution = pixelSize * iSCmatrixCollated[i, j, :]
                if norm == "nCells":
                    probability = len(np.nonzero(distanceDistribution < threshold)[0]) / nCells
                    # print('Using nCells normalisation')
                elif norm == "nonNANs":
                    numberNANs = len(np.nonzero(np.isnan(distanceDistribution))[0])
                    if nCells == numberNANs:
                        probability = 1
                    else:
                        probability = len(np.nonzero(distanceDistribution < threshold)[0]) / (nCells - numberNANs)
                    # print('Using NonNANs normalisation {}'.format(nCells-numberNANs))
                SCmatrix[i, j] = probability

    return SCmatrix, nCells


def buildsPWDmatrix(param,
    currentFolder, fileNameBarcodeCoordinates, outputFileName, dataFolder, pixelSize=0.1, logNameMD="log.md", ndims=2
):

    # Loads localAlignment if it exists
    localAlignmentFileName=dataFolder.outputFiles["alignImages"].split(".")[0] + "_localAlignment.dat"
    if os.path.exists(localAlignmentFileName):
        alignmentResultsTable= Table.read(localAlignmentFileName, format="ascii.ecsv")
        alignmentResultsTableRead=True
    else:
        print("\n\n *** Warning: could not found localAlignment: {}\n Proceeding with only global alignments...".format(localAlignmentFileName))
        alignmentResultsTableRead=False

    # Loads coordinate Tables        
    if os.path.exists(fileNameBarcodeCoordinates):
        barcodeMap = Table.read(fileNameBarcodeCoordinates, format="ascii.ecsv")
        if ndims==3 and "zcentroidGauss" in barcodeMap.keys():
            localizationDimension = 3
        else:
            localizationDimension = 2
    else:
        print("\n\n *** ERROR: could not found coordinates file: {}".format(fileNameBarcodeCoordinates))
        return -1
    
    # processes tables 
    barcodeMapROI = barcodeMap.group_by("ROI #")

    SCmatrixCollated, uniqueBarcodes = [], []
    numberROIs = len(barcodeMapROI.groups.keys)
    filesinFolder = glob.glob(currentFolder + os.sep + "*.tif")

    print("\nROIs detected: {}".format(barcodeMapROI.groups.keys))
    processingOrder = 0
    
    # loops over ROIs
    for ROI in range(numberROIs):
        nROI = barcodeMapROI.groups.keys[ROI][0]  # need to iterate over the first index

        print("Loading masks and pre-processing barcodes for ROI# {}".format(nROI))

        barcodeMapSingleROI = barcodeMap.group_by("ROI #").groups[ROI]

        # finds file with cell masks
        fileList2Process = [
            file
            for file in filesinFolder
            if file.split("_")[-1].split(".")[0] == "ch00"
            and "DAPI" in os.path.basename(file).split("_")
            and int(os.path.basename(file).split("_")[3]) == nROI
        ]

        if len(fileList2Process) > 0:
            
            # loads file with cell masks
            fileNameROImasks = os.path.basename(fileList2Process[0]).split(".")[0] + "_Masks.npy"
            fullFileNameROImasks = os.path.dirname(fileNameBarcodeCoordinates) + os.sep + fileNameROImasks
            if os.path.exists(fullFileNameROImasks):
                Masks = np.load(fullFileNameROImasks)

                # Assigns barcodes to Masks for a given ROI
                cellROI = cellID(param,barcodeMapSingleROI, Masks, ROI,ndims=localizationDimension)
                cellROI.ndims=ndims
                
                if alignmentResultsTableRead:
                    cellROI.alignmentResultsTable = alignmentResultsTable
                    
                # finds what barcodes are in each cell mask    
                cellROI.alignByMasking()

                # builds the single cell distance Matrix
                cellROI.buildsdistanceMatrix("min")  # mean min last

                print(
                    "ROI: {}, N cells assigned: {} out of {}\n".format(ROI, cellROI.NcellsAssigned-1, cellROI.numberMasks)
                )

                uniqueBarcodes = cellROI.uniqueBarcodes

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

            else:
                print(
                    "Error, no DAPI mask file found for ROI: {}, segmentedMasks: {}\n".format(
                        nROI, fileNameBarcodeCoordinates
                    )
                )
                print("File I was searching for: {}".format(fullFileNameROImasks))
                print("Debug: ")
                for file in filesinFolder:
                    if (
                        file.split("_")[-1].split(".")[0] == "ch00"
                        and "DAPI" in file.split("_")
                        and int(os.path.basename(file).split("_")[3]) == nROI
                    ):
                        print("Hit found!")
                    print(
                        "fileSplit:{}, DAPI in filename: {}, ROI: {}".format(
                            file.split("_")[-1].split(".")[0],
                            "DAPI" in os.path.basename(file).split("_"),
                            int(os.path.basename(file).split("_")[3]),
                        )
                    )

    # saves output
    np.save(outputFileName + "_HiMscMatrix.npy", SCmatrixCollated)
    np.savetxt(outputFileName + "_uniqueBarcodes.ecsv", uniqueBarcodes, delimiter=" ", fmt="%d")

    # plots outputs
    
    # adapts clim depending on whether 2 or 3 dimensions are use for the barcode localizations
    if localizationDimension==2:
        clim=1.6
    else:
        clim=1.6
        
    plotMatrix(
        SCmatrixCollated, uniqueBarcodes, pixelSize, numberROIs, outputFileName, logNameMD, mode="median", clim=clim, cm='terrain'
    )  # need to validate use of KDE. For the moment it does not handle well null distributions

    plotDistanceHistograms(SCmatrixCollated, pixelSize, outputFileName, logNameMD)


def processesPWDmatrices(param, log1, session1):
    sessionName = "buildsPWDmatrix"

    # processes folders and files
    dataFolder = folders(param.param["rootFolder"])
    log1.addSimpleText("\n===================={}====================\n".format(sessionName))
    log1.report("folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(log1.fileNameMD, "## {}\n".format(sessionName), "a")

    for currentFolder in dataFolder.listFolders:
        # filesFolder=glob.glob(currentFolder+os.sep+'*.tif')
        dataFolder.createsFolders(currentFolder, param)
        log1.report("-------> Processing Folder: {}".format(currentFolder))


        fileNameBarcodeCoordinates = dataFolder.outputFiles["segmentedObjects"] + "_barcode.dat"
        
        # 2D
        outputFileName = dataFolder.outputFiles["buildsPWDmatrix"]
        log1.report("-------> 2D processing: {}".format(outputFileName ))        
        
        if "pixelSizeXY" in param.param['acquisition'].keys():
            pixelSize = param.param['acquisition']['pixelSizeXY']
        else:
            pixelSize = 0.1
            
        buildsPWDmatrix(
            param,currentFolder, fileNameBarcodeCoordinates, outputFileName, dataFolder, pixelSize, log1.fileNameMD,
        )

        # 3D
        outputFileName = dataFolder.outputFiles["buildsPWDmatrix"]+"_3D"
        log1.report("-------> 3D processing: {}".format(outputFileName ))        

        if ("pixelSizeZ" in param.param['acquisition'].keys()) and ("pixelSizeXY" in param.param['acquisition'].keys()):
            pixelSizeXY = param.param['acquisition']['pixelSizeXY']
            pixelSizeZ = param.param['acquisition']['pixelSizeZ']
            
            # need to solve the issue of voxelsize...
            # pixelSize = [pixelSizeXY,pixelSizeXY,pixelSizeZ]
            pixelSize = pixelSizeXY
        else:
            pixelSize = 0.1            
            
        buildsPWDmatrix(
            param,currentFolder, fileNameBarcodeCoordinates, outputFileName, dataFolder, pixelSize, log1.fileNameMD, ndims=3
        )

        # loose ends
        session1.add(currentFolder, sessionName)

        log1.report("HiM matrix in {} processed".format(currentFolder), "info")


