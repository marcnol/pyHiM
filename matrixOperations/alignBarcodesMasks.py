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
import argparse

from datetime import datetime
import uuid
import numpy as np
import matplotlib.pyplot as plt

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

# =============================================================================
# CLASSES
# =============================================================================


class cellID:
    def __init__(self, barcodeMapROI, Masks, ROI):
        self.barcodeMapROI = barcodeMapROI
        self.Masks = Masks
        self.NcellsAssigned = 0
        self.NcellsUnAssigned = 0
        self.NbarcodesinMask = 0

        self.SegmentationMask = SegmentationImage(self.Masks)
        self.numberMasks = self.SegmentationMask.nlabels
        self.ROI = ROI
        self.alignmentResultsTable=Table()
        self.barcodesinMask = dict()
        for mask in range(self.numberMasks + 1):
            self.barcodesinMask["maskID_" + str(mask)] = []

    def visualize(self):
        pass
        # imageBarcodes = np.zeros([2048, 2048])
        # MasksBarcodes = Masks
        # R = []

        # for i in range(len(self.barcodeMapROI.groups[0])):
        #     y_int = int(self.barcodeMapROI.groups[0]["xcentroid"][i])
        #     x_int = int(self.barcodeMapROI.groups[0]["ycentroid"][i])
        #     barcodeID = self.barcodeMapROI.groups[0]["Barcode #"][i]
        #     imageBarcodes[x_int][y_int] = barcodeID
        #     MasksBarcodes[x_int][y_int] += 20 * barcodeID
        #     R.append([y_int, x_int, barcodeID])

        # # Shows results
        # Ra = np.array(R)
        # # plt.imshow(Masks, origin="lower", cmap="jet")
        # plt.scatter(Ra[:, 0], Ra[:, 1], s=5, c=Ra[:, 2], alpha=0.5)

    def alignByMasking(self):
        # [ Assigns barcodes to masks and creates <NbarcodesinMask> ]
        NbarcodesinMask = np.zeros(self.numberMasks + 2)
        print("ROI:{}".format(self.ROI))
        for i in range(len(self.barcodeMapROI.groups[0])):
            y_int = int(self.barcodeMapROI.groups[0]["xcentroid"][i])
            x_int = int(self.barcodeMapROI.groups[0]["ycentroid"][i])
            # barcodeID = self.barcodeMapROI.groups[ROI]['Barcode #'][i]
            maskID = self.Masks[x_int][y_int]
            self.barcodeMapROI["CellID #"][i] = maskID
            if maskID > 0:
                NbarcodesinMask[maskID] += 1
                self.barcodesinMask["maskID_" + str(maskID)].append(i)

        # Total number of masks assigned and not assigned
        self.NcellsAssigned = np.count_nonzero(NbarcodesinMask > 0)
        self.NcellsUnAssigned = self.numberMasks - self.NcellsAssigned

        # this list contains which barcodes are allocated to which masks
        self.NbarcodesinMask = NbarcodesinMask

    def buildsdistanceMatrix(self, mode="mean"):
        """
        

        
        """
        print("building distance matrix")
        barcodeMapROI = self.barcodeMapROI

        # [ builds SCdistanceTable ]

        # sorts Table by cellID
        barcodeMapROI_cellID = barcodeMapROI.group_by("CellID #")  # ROI data sorted by cellID
        ROIs, cellID, nBarcodes, barcodeIDs, p, cuid, buid = [], [], [], [], [], [], []

        # iterates over all cell masks in an ROI
        foundMatch=[]
        for key, group in zip(barcodeMapROI_cellID.groups.keys, barcodeMapROI_cellID.groups):
            if key["CellID #"] > 1:  # excludes cellID 0 as this is background
                if len(self.alignmentResultsTable)>0:
                     # applies local drift correction
                    x_uncorrected,y_uncorrected = np.array(group["xcentroid"].data), np.array(group["ycentroid"].data)
                    
                    # searches for local alignment shift for this mask in this ROI
                    _foundMatch=False
                    for row in self.alignmentResultsTable:
                        if row["ROI #"]==group["ROI #"].data[0] and row["CellID #"]==key["CellID #"]:
                            _foundMatch=True
                            x_corrected, y_corrected = x_uncorrected + row["shift_x"], y_uncorrected + row["shift_y"]
                            # print(row)
                            
                    # keeps uncorrected values if no match is found                            
                    if not _foundMatch:
                        print("Did not find match for CellID #{} in ROI #{}".format(key["CellID #"],group["ROI #"].data[0]))
                        x_corrected,y_corrected = x_uncorrected,y_uncorrected
                        foundMatch.append(False)                    
                    else:
                        foundMatch.append(True)
                        
                    R = np.column_stack((x_corrected,y_corrected,)) 
                else:
                    # does not apply local drift correction
                    R = np.column_stack((np.array(group["xcentroid"].data), np.array(group["ycentroid"].data),))

                ROIs.append(group["ROI #"].data[0])
                cellID.append(key["CellID #"])
                nBarcodes.append(len(group))
                barcodeIDs.append(group["Barcode #"].data)
                buid.append(group["Buid"].data)
                p.append(pairwise_distances(R))
                cuid.append(str(uuid.uuid4()))  # creates cell unique identifier
                # print("CellID #={}, nBarcodes={}".format(key['CellID #'],len(group)))
                
        print("Local correction applied to {}/{} barcodes in ROI {}".format(np.nonzero(foundMatch)[0].shape[0],len(foundMatch),group["ROI #"].data[0]))
        SCdistanceTable = Table()  # [],names=('CellID', 'barcode1', 'barcode2', 'distances'))
        SCdistanceTable["Cuid"] = cuid
        SCdistanceTable["ROI #"] = ROIs
        SCdistanceTable["CellID #"] = cellID
        SCdistanceTable["nBarcodes"] = nBarcodes
        SCdistanceTable["Barcode #"] = barcodeIDs
        SCdistanceTable["Buid"] = buid
        SCdistanceTable["PWDmatrix"] = p

        # NEEED TO ADD THE IDENTITIES OF THE BARCODES

        self.SCdistanceTable = SCdistanceTable

        print("Cells with barcodes found: {}".format(len(SCdistanceTable)))

        # [ builds SCmatrix ]
        numberMatrices = len(SCdistanceTable)  # z dimensions of SCmatrix
        uniqueBarcodes = np.unique(barcodeMapROI["Barcode #"].data)
        # number of unique Barcodes for xy dimensions of SCmatrix
        numberUniqueBarcodes = uniqueBarcodes.shape[0]
        SCmatrix = np.zeros((numberUniqueBarcodes, numberUniqueBarcodes, numberMatrices))
        SCmatrix[:] = np.NaN

        for iCell, scPWDitem in zip(range(numberMatrices), SCdistanceTable):
            barcodes2Process = scPWDitem["Barcode #"]
            for barcode1, ibarcode1 in zip(barcodes2Process, range(len(barcodes2Process))):
                indexBarcode1 = np.nonzero(uniqueBarcodes == barcode1)[0][0]
                for barcode2, ibarcode2 in zip(barcodes2Process, range(len(barcodes2Process))):
                    indexBarcode2 = np.nonzero(uniqueBarcodes == barcode2)[0][0]
                    if barcode1 != barcode2:
                        newdistance = scPWDitem["PWDmatrix"][ibarcode1][ibarcode2]
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
def findsOptimalKernelWidth(distanceDistribution):
    bandwidths = 10 ** np.linspace(-1, 1, 100)
    grid = GridSearchCV(KernelDensity(kernel="gaussian"), {"bandwidth": bandwidths}, cv=LeaveOneOut())
    grid.fit(distanceDistribution[:, None])
    return grid.best_params_


def retrieveKernelDensityEstimator(distanceDistribution0, x_d, optimizeKernelWidth=False):

    nan_array = np.isnan(distanceDistribution0)

    not_nan_array = ~nan_array

    distanceDistribution = distanceDistribution0[not_nan_array]

    # instantiate and fit the KDE model
    if optimizeKernelWidth:
        kernelWidth = findsOptimalKernelWidth(distanceDistribution)["bandwidth"]
    else:
        kernelWidth = 0.3

    kde = KernelDensity(bandwidth=kernelWidth, kernel="gaussian")
    kde.fit(distanceDistribution[:, None])

    # score_samples returns the log of the probability density
    logprob = kde.score_samples(x_d[:, None])

    return logprob, distanceDistribution


def distributionMaximumKernelDensityEstimation(SCmatrixCollated, bin1, bin2, pixelSize, optimizeKernelWidth=False):
    distanceDistribution0 = pixelSize * SCmatrixCollated[bin1, bin2, :]
    x_d = np.linspace(0, 5, 2000)

    # checks that distribution is not empty
    if distanceDistribution0.shape[0] > 0:
        logprob, distanceDistribution = retrieveKernelDensityEstimator(distanceDistribution0, x_d, optimizeKernelWidth)
        kernelDistribution = 10 * np.exp(logprob)
        maximumKernelDistribution = x_d[np.argmax(kernelDistribution)]
        return maximumKernelDistribution, distanceDistribution, kernelDistribution, x_d
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
                print("Dataset {} cells2plot: {}".format(figtitle, nCells))
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


def buildsPWDmatrix(
    currentFolder, fileNameBarcodeCoordinates, outputFileName, dataFolder, pixelSize=0.1, logNameMD="log.md",
):

    # Loads localAlignment and coordinate Tables
    localAlignmentFileName=dataFolder.outputFiles["alignImages"].split(".")[0] + "_localAlignment.dat"
    if os.path.exists(localAlignmentFileName):
        alignmentResultsTable= Table.read(localAlignmentFileName, format="ascii.ecsv")
        alignmentResultsTableRead=True
    else:
        print("\n\n *** Warning: could not found localAlignment: {}\n Proceeding with only global alignments...".format(localAlignmentFileName))
        alignmentResultsTableRead=False
        
    if os.path.exists(fileNameBarcodeCoordinates):
        barcodeMap = Table.read(fileNameBarcodeCoordinates, format="ascii.ecsv")
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
    for ROI in range(numberROIs):
        nROI = barcodeMapROI.groups.keys[ROI][0]  # need to iterate over the first index

        print("Working on ROI# {}".format(nROI))

        barcodeMapSingleROI = barcodeMap.group_by("ROI #").groups[ROI]

        # finds file for masks
        fileList2Process = [
            file
            for file in filesinFolder
            if file.split("_")[-1].split(".")[0] == "ch00"
            and "DAPI" in os.path.basename(file).split("_")
            and int(os.path.basename(file).split("_")[3]) == nROI
        ]

        if len(fileList2Process) > 0:
            # loads Masks
            fileNameROImasks = os.path.basename(fileList2Process[0]).split(".")[0] + "_Masks.npy"
            fullFileNameROImasks = os.path.dirname(fileNameBarcodeCoordinates) + os.sep + fileNameROImasks
            if os.path.exists(fullFileNameROImasks):
                Masks = np.load(fullFileNameROImasks)

                # Assigns barcodes to Masks for a given ROI
                cellROI = cellID(barcodeMapSingleROI, Masks, ROI)
                
                if alignmentResultsTableRead:
                    cellROI.alignmentResultsTable = alignmentResultsTable
                    
                cellROI.alignByMasking()

                cellROI.buildsdistanceMatrix("min")  # mean min last

                print(
                    "ROI: {}, N cells assigned: {} out of {}".format(ROI, cellROI.NcellsAssigned, cellROI.numberMasks)
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
    plotMatrix(
        SCmatrixCollated, uniqueBarcodes, pixelSize, numberROIs, outputFileName, logNameMD, mode="median",
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
        # segmentedMasksFolder=dataFolder.outputFolders['segmentedObjects']
        outputFileName = dataFolder.outputFiles["buildsPWDmatrix"]
        pixelSize = 0.1

        buildsPWDmatrix(
            currentFolder, fileNameBarcodeCoordinates, outputFileName, dataFolder, pixelSize, log1.fileNameMD,
        )
        session1.add(currentFolder, sessionName)

        log1.report("HiM matrix in {} processed".format(currentFolder), "info")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = os.getcwd()
        rootFolder = "."
        rootFolder = "/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug"

    print("parameters> rootFolder: {}".format(rootFolder))
    now = datetime.now()

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    # session
    session1 = session(rootFolder, "processingPipeline")

    # setup logs
    log1 = log(rootFolder)
    log1.addSimpleText("\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format("processingPipeline"))
    log1.report("Hi-M analysis MD: {}".format(log1.fileNameMD))
    writeString2File(
        log1.fileNameMD, "# Hi-M analysis {}".format(now.strftime("%Y/%m/%d %H:%M:%S")), "w",
    )  # initialises MD file

    for ilabel in range(len(labels2Process)):
        label = labels2Process[ilabel]["label"]
        labelParameterFile = labels2Process[ilabel]["parameterFile"]
        log1.addSimpleText("**Analyzing label: {}**".format(label))

        # sets parameters
        param = Parameters(rootFolder, labelParameterFile)

        # [builds PWD matrix for all folders with images]
        if label == "DAPI":
            processesPWDmatrices(param, log1, session1)
            
    # parser = argparse.ArgumentParser()
    # parser.add_argument("-F", "--rootFolder", help="Folder with images")
    # args = parser.parse_args()

    # print("\n--------------------------------------------------------------------------")

    # if args.rootFolder:
    #     rootFolder = args.rootFolder
    # else:
    #     rootFolder = "."
    #     # rootFolder='/home/marcnol/data/Experiment_15/Embryo_006_ROI18'
    #     # rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset'
    # print("parameters> rootFolder: {}".format(rootFolder))

    # # sets name for coordinate file and for masks
    # segmentedMasksFolder = rootFolder + "/rawData/segmentedObjects/"
    # fileNameBarcodeCoordinates = segmentedMasksFolder + "segmentedObjects_barcode.dat"
    # # currentFolder=glob.glob(rootFolder+'/rawData/'+'*.tif')
    # currentFolder = rootFolder + "/rawData"

    # outputFileName = rootFolder.split("/")[-1]
    # pixelSize = 0.1
    # dataFolder = folders(param.param["rootFolder"])
    
    # buildsPWDmatrix(currentFolder, fileNameBarcodeCoordinates, outputFileName, dataFolder, pixelSize)