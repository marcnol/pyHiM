#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 21:21:02 2020

@author: marcnol
"""
import numpy as np
import os
import glob

from astropy.table import Table, vstack

from fileManagement import writeString2File
from alignBarcodesMasks import plotMatrix

def attributesLabels2cells(SNDtable, ResultsTable, label="doc"):

    sortedSNDTable = SNDtable.group_by("MaskID #")
    listKeys = list(sortedSNDTable.groups.keys["MaskID #"].data)
    indexKey = [index for i, index in zip(listKeys, range(len(listKeys))) if i == label]
    SNDTablewithLabel = sortedSNDTable.groups[indexKey[0]]
    print("\n>>> Matching labels")
    print("Found {} out of {} cells with {} in dataset".format(len(SNDTablewithLabel), len(sortedSNDTable), label))

    # sorts Results Table by ROI
    PWDTableSortedROI = ResultsTable.group_by("ROI #")

    CUIDs = []

    for ROI, group in zip(PWDTableSortedROI.groups.keys, PWDTableSortedROI.groups):
        # list of cellIDs in ROI
        cells2Process = group["CellID #"].data.compressed()
        cells2ProcessUID = group["Cuid"]

        # list of ROIs detected in SNDTablewithLabel
        ROIsinSNDTablewithLabel = list(SNDTablewithLabel.group_by("ROI #").groups.keys["ROI #"].data)

        # index of ROI within the keys of SNDTablewithLabel
        indexROIs = [
            index for i, index in zip(ROIsinSNDTablewithLabel, range(len(ROIsinSNDTablewithLabel))) if i == ROI["ROI #"]
        ]

        # subtable of cells with label and ROI that we are looking for
        SNDTablewithLabelROI = SNDTablewithLabel.group_by("ROI #").groups[indexROIs[0]]
        cellswithLabel = list(SNDTablewithLabelROI["CellID #"].data)

        # finds which cell indeces in Table have label
        listofSelectedCells = [
            index for iCell, index in zip(cells2Process, range(len(cells2Process))) if iCell in cellswithLabel
        ]

        if len(CUIDs) > 0:
            CUIDs = vstack([CUIDs, cells2ProcessUID[listofSelectedCells]])
        else:
            CUIDs = cells2ProcessUID[listofSelectedCells]

        print(
            "Processed ROI # {}, found {} out of {} cells with {}".format(
                ROI["ROI #"], len(listofSelectedCells), len(group), label
            )
        )

    # from list of CUIDs from cells that show label, I construct a binary vector of the same size as SCmatrix. Labeled cells have a 1.
    SClabeled = np.zeros(len(ResultsTable))
    CUIDsList = CUIDs["Cuid"].data.compressed()
    indexCellsWithLabel = [
        iRow for Row, iRow in zip(ResultsTable, range(len(ResultsTable))) if Row["Cuid"] in CUIDsList
    ]
    SClabeled[indexCellsWithLabel] = 1

    return SClabeled, CUIDsList


def loadsSCdata(ListData, datasetName, p):
    """
    loads SC datasets from a dict of folders (ListData)

    Parameters
    ----------
    ListData : dict
        dict of folders with data that can be loaded.
    dataset2Load : int, optional
        The item in the dictionary that will be loaded. The default is 3.

    Returns
    -------
    SCmatrixCollated : list of np arrays nBarcodes x nBarcodes x nCells
        Cummulative SC PWD matrix.
    uniqueBarcodes : list of np arrays 
        containing the barcode identities for each matrix.
    buildsPWDmatrixCollated : list of Tables
        Tables with all the data for cells and barcodes used to produce SCmatrixCollated.

    """
    # tags2process = list(ListData.keys())
    print("Dataset to load: {}\n\n".format(list(ListData.keys())[0]))

    SCmatrixCollated, uniqueBarcodes = [], []
    buildsPWDmatrixCollated, runName, SClabeledCollated = [], [], []

    for rootFolder in ListData[datasetName]["Folders"]:
        # [finds sample name]
        runName.append(os.path.basename(os.path.dirname(rootFolder)))

        # [makes list of files with Tables to load]
        # tries to load files from newer version of proceesingPipeline.py
        files2Process = glob.glob(rootFolder + "/buildsPWDmatrix_order*ROI*.ecsv")
        if len(files2Process) == 0:
            # it resorts to old format
            files2Process = glob.glob(rootFolder + "/buildsPWDmatrix_*ROI*.ecsv")
        else:
            print("Found {} ECSV files in {}".format(len(files2Process), rootFolder))

        # checks that something was found
        if len(files2Process) > 0:
            print(">>> Loading {} results tables".format(len(files2Process)))

            # [initializes variables]
            buildsPWDmatrix = Table()
            fileOrder, fileTimeStamp = (
                np.zeros(len(files2Process), dtype=int),
                np.zeros(len(files2Process)),
            )

            # [finds what order Table files should be loaded to agree with order in buildsPWDmatrix_HiMscMatrix.npy]
            for fileName, ifileName in zip(files2Process, range(len(files2Process))):
                if "_order" in fileName:
                    for isplit in fileName.split("_"):
                        if "order" in isplit:
                            fileOrder[ifileName] = int(isplit.split(":")[1])
                            print("order= {}--> {}".format(os.path.basename(fileName), fileOrder[ifileName]))
                    fileTimeStamp[ifileName] = os.path.getmtime(fileName)
                    choosingTimeStamp = False
                else:
                    fileTimeStamp[ifileName] = os.path.getmtime(fileName)
                    choosingTimeStamp = True

            if choosingTimeStamp:
                fileOrder = np.argsort(fileTimeStamp).astype(int)

            # [loads buildsPWDmatrix Tables]
            for ifileName in range(len(files2Process)):
                fileName = files2Process[fileOrder[ifileName]]
                newbuildsPWDmatrix = Table.read(fileName, format="ascii.ecsv")  # ascii.ecsv
                buildsPWDmatrix = vstack([buildsPWDmatrix, newbuildsPWDmatrix])
                print(
                    "[{}:{}] From {}, Read: {} cells, Cummulative: {} cells".format(
                        fileOrder[ifileName],
                        fileTimeStamp[fileOrder[ifileName]],
                        os.path.basename(fileName),
                        len(newbuildsPWDmatrix),
                        len(buildsPWDmatrix),
                    )
                )

            # [loads SNDassignedCells.ecsv files if available]
            fileNameSNDassignedCells = os.path.dirname(rootFolder) + os.sep + "segmentedObjects/SNDassignedCells.ecsv"
            if os.path.exists(fileNameSNDassignedCells):
                SNDassignedCells = Table.read(fileNameSNDassignedCells, format="ascii.ecsv")

                # attributes masks to single cells
                SClabeled, CUIDsList = attributesLabels2cells(SNDassignedCells, buildsPWDmatrix, label=p["label"])

                SClabeledCollated.append(SClabeled)
            else:
                # if not available it makes a mock SClabeled matrix so that pipeline always works
                SClabeled = np.ones(len(buildsPWDmatrix)).astype(int)
                SClabeledCollated.append(SClabeled)

            # [loads and accumulates barcodes and scHiM matrix]
            fileNamMatrix = rootFolder + os.sep + "buildsPWDmatrix_HiMscMatrix.npy"
            fileNameBarcodes = rootFolder + os.sep + "buildsPWDmatrix_uniqueBarcodes.ecsv"

            if os.path.exists(fileNamMatrix):
                SCmatrix1 = np.load(fileNamMatrix)
                SCmatrixCollated.append(SCmatrix1)
            else:
                print("*** Error: could not find {}".format(fileNamMatrix))

            if os.path.exists(fileNameBarcodes):
                uniqueBarcodes.append(np.loadtxt(fileNameBarcodes).astype(int))
            else:
                print("*** Error: could not find {}".format(fileNameBarcodes))

            buildsPWDmatrixCollated.append(buildsPWDmatrix)

            print("\n>>>Merging rootFolder: {}".format(rootFolder))
            print("Cells added after merge: {}\n".format(SCmatrix1.shape[2]))
        else:
            print("No file detected in the folder you provide: {}".format(rootFolder + "/buildsPWDmatrix_*ROI*.ecsv"))
    print("{} datasets loaded\n".format(len(SCmatrixCollated)))

    return (
        SCmatrixCollated,
        uniqueBarcodes,
        buildsPWDmatrixCollated,
        runName,
        SClabeledCollated,
    )


def listsSCtoKeep(p, mask):
    # print('{}'.format(p['action']))
    if p["action"] == "all":
        try:
            cells2Plot = range(len(mask))
        except TypeError:
            print(mask)
            cells2Plot = range(mask.shape[0])
    elif p["action"] == "labeled":
        a = [i for i in range(len(mask)) if mask[i] == 1]
        cells2Plot = a
    elif p["action"] == "unlabeled":
        a = [i for i in range(len(mask)) if mask[i] == 0]
        cells2Plot = a

    return cells2Plot


def normalizeMatrix(SCmatrix_wt):
    SCmatrix_wt_normalized = SCmatrix_wt
    nBins = SCmatrix_wt.shape[0]

    for iRow in range(nBins):
        rowSum = np.sum(SCmatrix_wt[iRow, :])
        for iCol in range(nBins):
            SCmatrix_wt_normalized[iRow, iCol] = SCmatrix_wt_normalized[iRow, iCol] / rowSum
            SCmatrix_wt_normalized[iCol, iRow] = SCmatrix_wt_normalized[iCol, iRow] / rowSum
    return SCmatrix_wt_normalized

def plotsEnsemble3wayContactMatrix(
    SCmatrixCollated, uniqueBarcodes, anchors, sOut, runName, iListData, p,fileNameMD="tmp.md", datasetName=""):

    # combines matrices from different embryos and calculates integrated contact probability matrix


    SCmatrixAllDatasets = []  # np.zeros((nBarcodes,nBarcodes))
    for iSCmatrixCollated, iuniqueBarcodes, mask, iTag in zip(SCmatrixCollated, uniqueBarcodes,p['SClabeledCollated'],runName):
        cells2Plot = listsSCtoKeep(p, mask)

        if max(cells2Plot)>iSCmatrixCollated.shape[2]:
            print('Error: max in cells2plot {} in dataset {} is larger than the number of available cells {}'.format(max(cells2Plot),iTag,iSCmatrixCollated.shape[2]))
        else:
            if len(SCmatrixAllDatasets) > 0:
                SCmatrixAllDatasets = np.concatenate((SCmatrixAllDatasets, iSCmatrixCollated[:,:,cells2Plot]), axis=2)
            else:
                SCmatrixAllDatasets = iSCmatrixCollated[:,:,cells2Plot]

            commonSetUniqueBarcodes = iuniqueBarcodes


    # print(commonSetUniqueBarcodes)
    for anchor in anchors:
        print('nCells processed: {}'.format(SCmatrixAllDatasets.shape[2]))
        SCmatrix = calculate3wayContactMatrix(
            SCmatrixAllDatasets,
            uniqueBarcodes,
            p['pixelSize'],
            anchor,
            sOut,
            threshold=iListData["ContactProbability_distanceThreshold"],
            norm="nonNANs",
        )  # norm: nonNANs (default)

        outputFileName = p['outputFolder'] + os.sep + datasetName + "_Cells:" + p['action'] + "_ensemble3wayContacts"
        outputFileName += "_anchor_" + str(anchor)
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")

        cScale = np.max(SCmatrix)
        # print(cScale)
        # print(SCmatrix.shape)

        plotMatrix(
            SCmatrix,
            uniqueBarcodes,
            p['pixelSize'],
            cm=iListData["ContactProbability_cm"],
            outputFileName=outputFileName,
            clim=cScale,
            cMin=iListData["ContactProbability_cmin"],
            figtitle="3way contacts",
            cmtitle=sOut,
            nCells=0,
            mode="counts"
        )  # twilight_shifted_r

def calculate3wayContactMatrix(iSCmatrixCollated, iuniqueBarcodes, pixelSize,
        anchor, sOut, threshold=0.25, norm="nonNANs"):
    
    nX = nY = iSCmatrixCollated.shape[0]
    SCmatrix = np.zeros((nX, nY))

    # transform distance matrix from pixel to µm
    mat = pixelSize * iSCmatrixCollated

    # print(nX, nY)
    for bait1 in range(nX):
        for bait2 in range(nY):
            if ( bait1 == bait2 ):
                continue

            # print("current bait1", bait1, "bait2", bait2)
            n_contacts, n_nonNaN = getMultiContact(mat,anchor,bait1,bait2,threshold)
            if ( sOut == "Counts"):
                SCmatrix[bait1,bait2] = n_contacts
            elif (sOut == "Probability"):
                SCmatrix[bait1,bait2] = n_contacts / n_nonNaN;
            else:
                print("Unexpected sOut.")
                return -1

            # print(n_contacts / n_nonNaN)
            # print(type(n_contacts), type(n_nonNaN))

    SCmatrix[np.isnan(SCmatrix)] = 0 # set NaN to zero
    return SCmatrix


def getMultiContact(mat,anchor,bait1,bait2,threshold):
    """
    Input:
    mat        : pwd matrix, including only the bins of used RTs
    anchor     : anchor bin
    bait1      : bait bin #1
    bait2      : bait bin #2
    threshold  : contact threshold
    Output:
    n_contacts : number of contacts between bins anchor, bait1, and bait2
    n_nonNaN   : number of cells where the distances anchor-bait1 and anchor-bait2 are present
    """

    A = mat[anchor,bait1,:]
    B = mat[anchor,bait2,:]

    # get fraction of points in quadrant
    n1 = np.sum( (A<threshold) & (B<threshold) );
    totN = np.sum( (~np.isnan(A)) & (~np.isnan(B)) );


    return n1, totN

