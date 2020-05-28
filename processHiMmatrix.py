#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:36:20 2020

@author: marcnol


This script takes JSON file with folders where datasets are 
stored and processes multiple PWD matrices together.

$ processHiMmatrix.py -F rootFolder


"""

# =============================================================================
# IMPORTS
# =============================================================================q

import numpy as np
import glob, time
import os
import json
from datetime import datetime
import argparse

from alignBarcodesMasks import plotDistanceHistograms, plotMatrix
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import LeaveOneOut
from alignBarcodesMasks import (
    distributionMaximumKernelDensityEstimation,
    calculateContactProbabilityMatrix,
)
from astropy.table import Table, Column, vstack
from fileManagement import writeString2File

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")

# =============================================================================
# FUNCTIONS
# =============================================================================q


def normalizeMatrix(SCmatrix_wt):
    SCmatrix_wt_normalized = SCmatrix_wt
    nBins = SCmatrix_wt.shape[0]

    for iRow in range(nBins):
        rowSum = np.sum(SCmatrix_wt[iRow, :])
        for iCol in range(nBins):
            SCmatrix_wt_normalized[iRow, iCol] = SCmatrix_wt_normalized[iRow, iCol] / rowSum
            SCmatrix_wt_normalized[iCol, iRow] = SCmatrix_wt_normalized[iCol, iRow] / rowSum
    return SCmatrix_wt_normalized


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
        indexROIs = [index for i, index in zip(ROIsinSNDTablewithLabel, range(len(ROIsinSNDTablewithLabel))) if i == ROI["ROI #"]]

        # subtable of cells with label and ROI that we are looking for
        SNDTablewithLabelROI = SNDTablewithLabel.group_by("ROI #").groups[indexROIs[0]]
        cellswithLabel = list(SNDTablewithLabelROI["CellID #"].data)

        # finds which cell indeces in Table have label
        listofSelectedCells = [index for iCell, index in zip(cells2Process, range(len(cells2Process))) if iCell in cellswithLabel]

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
    indexCellsWithLabel = [iRow for Row, iRow in zip(ResultsTable, range(len(ResultsTable))) if Row["Cuid"] in CUIDsList]
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

        # [finds loading order]
        files2Process = glob.glob(rootFolder + "/buildsPWDmatrix_*ROI*.ecsv")
        print(">>> Loading {} results tables".format(len(files2Process)))
        buildsPWDmatrix = Table()
        fileOrder, fileTimeStamp = (
            np.zeros(len(files2Process)),
            np.zeros(len(files2Process)),
        )
        for fileName, ifileName in zip(files2Process, range(len(files2Process))):
            if "_order" in fileName:
                for isplit in fileName.split("_"):
                    if "order" in isplit:
                        fileOrder[ifileName] = int(isplit.split(":")[1])
                        print("order= {}--> {}".format(os.path.basename(fileName), fileOrder[ifileName]))
                choosingTimeStamp = False
            else:
                fileTimeStamp[ifileName] = os.path.getmtime(fileName)
                choosingTimeStamp = True

        if choosingTimeStamp:
            fileOrder = np.argsort(fileTimeStamp)

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

        # [loads SNDassignedCells.ecsv files]
        fileNameSNDassignedCells = os.path.dirname(rootFolder) + os.sep + "segmentedObjects/SNDassignedCells.ecsv"
        if os.path.exists(fileNameSNDassignedCells):
            SNDassignedCells = Table.read(fileNameSNDassignedCells, format="ascii.ecsv")

            # attributes masks to single cells
            SClabeled, CUIDsList = attributesLabels2cells(SNDassignedCells, buildsPWDmatrix, label=p['label'])

            SClabeledCollated.append(SClabeled)
        else:
            SClabeled=np.ones(len(buildsPWDmatrix)).astype(int)
            SClabeledCollated.append(SClabeled)
            
        # [loads and accumulates barcodes and scHiM matrix]

        fileNamMatrix = rootFolder + os.sep + "buildsPWDmatrix_HiMscMatrix.npy"
        fileNameBarcodes = rootFolder + os.sep + "buildsPWDmatrix_uniqueBarcodes.ecsv"

        SCmatrix1 = np.load(fileNamMatrix)
        SCmatrixCollated.append(SCmatrix1)
        uniqueBarcodes.append(np.loadtxt(fileNameBarcodes).astype(int))
        buildsPWDmatrixCollated.append(buildsPWDmatrix)
        print("\n>>>Merging rootFolder: {}".format(rootFolder))
        print("Cells added after merge: {}\n".format(SCmatrix1.shape[2]))

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
    if p['action']=='all':
        cells2Plot = range(len(mask))
    elif p['action']=='labeled':
        a=[i for i in range(len(mask)) if mask[i]==1]
        cells2Plot = a
    elif p['action']=='unlabeled':
        a=[i for i in range(len(mask)) if mask[i]==0]
        cells2Plot = a
        
    
    return cells2Plot


def plotsSinglePWDmatrices(
    SCmatrixCollated, uniqueBarcodes, runName, iListData, p, fileNameMD="tmp.md", datasetName="",
):
    # plots distance matrix for each dataset
    for iSCmatrixCollated, iuniqueBarcodes, iTag, mask in zip(SCmatrixCollated, uniqueBarcodes, runName,p['SClabeledCollated']):
        outputFileName = outputFolder + os.sep + iTag + "_Cells:" + p['action'] +"_PWDmatrix"

        cells2Plot = listsSCtoKeep(p, mask)
        
        plotMatrix(
            iSCmatrixCollated,
            iuniqueBarcodes,
            pixelSize,
            outputFileName=outputFileName,
            figtitle="PWD:" + iTag,
            cm="terrain",
            clim=iListData["PWD_clim"],
            mode=iListData["PWD_mode"],
            nCells=iSCmatrixCollated.shape[2],
            cells2Plot=cells2Plot
        )  # twilight_shifted_r 1.4, mode: median KDE
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def plotsInversePWDmatrice(SCmatrixCollated, uniqueBarcodes, runName, iListData, p,fileNameMD, datasetName=""):
    # plots inverse distance matrix for each dataset
    for iSCmatrixCollated, iuniqueBarcodes, iTag, mask in zip(SCmatrixCollated, uniqueBarcodes, runName,p['SClabeledCollated']):
        outputFileName = outputFolder + os.sep + iTag + "_Cells:" + p['action'] +"_invPWDmatrix"

        cells2Plot = listsSCtoKeep(p, mask)

        plotMatrix(
            iSCmatrixCollated,
            iuniqueBarcodes,
            pixelSize,
            cm="terrain",
            outputFileName=outputFileName,
            clim=iListData["iPWD_clim"],
            mode=iListData["iPWD_mode"],
            figtitle="inverse PWD:" + iTag,
            cmtitle="inverse distance, 1/nm",
            inverseMatrix=True,
            nCells=iSCmatrixCollated.shape[2],
            cells2Plot=cells2Plot
        )  # twilight_shifted_r, mode: median KDE
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def plotsSingleContactProbabilityMatrix(
    SCmatrixCollated, uniqueBarcodes, runName, iListData, p,fileNameMD="tmp.md", datasetName="",
):
    # Plots contact probability matrices for each dataset

    for iSCmatrixCollated, iuniqueBarcodes, iTag, mask in zip(SCmatrixCollated, uniqueBarcodes, runName,p['SClabeledCollated']):
        cells2Plot = listsSCtoKeep(p, mask)

        SCmatrix, nCells = calculateContactProbabilityMatrix(
            iSCmatrixCollated[:,:,cells2Plot],
            iuniqueBarcodes,
            pixelSize,
            threshold=iListData["ContactProbability_distanceThreshold"],
            norm="nonNANs",
        )  # norm: nCells (default), nonNANs
        outputFileName = outputFolder + os.sep + iTag + "_Cells:" + p['action'] +"_contactProbability"

        print('{}'.format(len(cells2Plot)))
        cScale = SCmatrix.max() / iListData["ContactProbability_scale"]
        
        plotMatrix(
            SCmatrix,
            iuniqueBarcodes,
            pixelSize,
            cm="terrain",
            outputFileName=outputFileName,
            cMin=iListData["ContactProbability_cmin"],
            clim=cScale,
            figtitle="HiM:" + iTag,
            cmtitle="probability",
            nCells=nCells,
            cells2Plot=cells2Plot
        )  # twilight_shifted_r
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def plotsEnsembleContactProbabilityMatrix(
    SCmatrixCollated, uniqueBarcodes, runName, iListData, p,fileNameMD="tmp.md", datasetName="",
):

    # print('{}'.format(iListData))
    # combines matrices from different embryos and calculates integrated contact probability matrix

    # nBarcodes = uniqueBarcodes[0].shape[0]
    SCmatrixAllDatasets = []  # np.zeros((nBarcodes,nBarcodes))

    for iSCmatrixCollated, iuniqueBarcodes,mask in zip(SCmatrixCollated, uniqueBarcodes,p['SClabeledCollated']):
        cells2Plot = listsSCtoKeep(p, mask)

        if len(SCmatrixAllDatasets) > 0:
            SCmatrixAllDatasets = np.concatenate((SCmatrixAllDatasets, iSCmatrixCollated[:,:,cells2Plot]), axis=2)
        else:
            SCmatrixAllDatasets = iSCmatrixCollated[:,:,cells2Plot]

        commonSetUniqueBarcodes = iuniqueBarcodes

    print('nCells processed: {}'.format(SCmatrixAllDatasets.shape[2]))
    SCmatrix, nCells = calculateContactProbabilityMatrix(
        SCmatrixAllDatasets,
        commonSetUniqueBarcodes,
        pixelSize,
        threshold=iListData["ContactProbability_distanceThreshold"],
        norm="nonNANs",
    )  # norm: nCells (default), nonNANs
    cScale = SCmatrix.max() / iListData["ContactProbability_scale"]
    outputFileName = outputFolder + os.sep + datasetName + "_Cells:" + p['action'] + "_ensembleContactProbability"
    writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")

    plotMatrix(
        SCmatrix,
        iuniqueBarcodes,
        pixelSize,
        cm="terrain",
        outputFileName=outputFileName,
        clim=cScale,
        cMin=iListData["ContactProbability_cmin"],
        figtitle="HiM counts",
        cmtitle="probability",
        nCells=nCells,
    )  # twilight_shifted_r

    np.savetxt(
        outputFolder + os.sep + "CombinedMatrix" + ":"+ list(ListData.keys())[0] + "_Cells:" + p['action'] + ".dat",
        SCmatrix,
        fmt="%.4f",
        delimiter=" ",
        newline="\n",
        header="Combined contact probability matrix",
        footer="",
        comments="# ",
        encoding=None,
    )

    np.savetxt(
        outputFolder + os.sep + "UniqueBarcodes" + ":"+ list(ListData.keys())[0] + "_Cells:" + p['action'] + ".dat",
        iuniqueBarcodes,
        fmt="%.4f",
        delimiter=" ",
        newline="\n",
        header="unique barcodes",
        footer="",
        comments="# ",
        encoding=None,
    )


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    pixelSize = 0.1

    # [parsing arguments]
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument(
        "-P", "--parameters", help="Provide name of parameter files. folders2Load.json assumed as default",
    )
    parser.add_argument("-A", "--label", help="Add name of label (e.g. doc)")
    parser.add_argument("-W", "--action", help="Select: [all], [labeled] or [unlabeled] cells plotted ")
    
    args = parser.parse_args()
    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = "."
    
    p={}
    
    if args.parameters:
        p['parametersFileName'] = args.parameters
    else:
        p['parametersFileName']= "folders2Load.json"

    if args.label:
        p['label']= args.label
    else:
        p['label']= "doc"

    if args.action:
        p['action']= args.action
    else:
        p['action']= "all"

    # [ initialises MD file]
    now = datetime.now()
    dateTime = now.strftime("%d%m%Y_%H%M%S")
    fileNameRoot = "processHiMmatrixAnalysis_"

    # [ Lists and loads datasets from different embryos]
    fileNameListDataJSON = rootFolder + os.sep + p['parametersFileName']
    print("\n--------------------------------------------------------------------------")
    if os.path.exists(fileNameListDataJSON):
        with open(fileNameListDataJSON) as json_file:
            ListData = json.load(json_file)
        print("Loaded JSON file with {} datasets from {}\n".format(len(ListData), fileNameListDataJSON))

    # [ creates output folder]
    outputFolder = rootFolder + os.sep + "scHiMmatrices"
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
        print("Folder created: {}".format(outputFolder))

    # [loops over lists of datafolders]
    for datasetName in list(ListData.keys()):

        # [loads SC matrices]
        (SCmatrixCollated, uniqueBarcodes, buildsPWDmatrixCollated, runName, SClabeledCollated,) = loadsSCdata(
            ListData, datasetName, p
        )
        fileNameMD = rootFolder + os.sep + fileNameRoot + "_" + datasetName + "_Cells:" + p['action'] + "_" + dateTime + ".md"
        writeString2File(fileNameMD, "# Post-processing of Hi-M matrices", "w")
        writeString2File(fileNameMD, "**dataset: {}** - **Cells: {}**".format(datasetName, p['action'] ), "a")
        p['SClabeledCollated']=SClabeledCollated
        
        # [plots distance matrix for each dataset]
        writeString2File(fileNameMD, "## single cell PWD matrices", "a")
        print(">>> Producing {} PWD matrices for dataset {}\n".format(len(SCmatrixCollated), datasetName))
        plotsSinglePWDmatrices(
            SCmatrixCollated, uniqueBarcodes, runName, ListData[datasetName], p, fileNameMD, datasetName=datasetName,
        )

        # plots histograms for each dataset
        # for iSCmatrixCollated, iuniqueBarcodes in zip(SCmatrixCollated, uniqueBarcodes):
        #     plotDistanceHistograms(iSCmatrixCollated, pixelSize, mode="KDE", limitNplots=15)

        # [plots inverse distance matrix for each dataset]
        writeString2File(fileNameMD, "## single cell inverse PWD matrices", "a")
        print(">>> Producing {} inverse PWD matrices for dataset {}\n".format(len(SCmatrixCollated), datasetName))
        plotsInversePWDmatrice(
            SCmatrixCollated, uniqueBarcodes, runName, ListData[datasetName], p,fileNameMD, datasetName=datasetName,
        )

        # [Plots contact probability matrices for each dataset]
        writeString2File(fileNameMD, "## single cell Contact Probability matrices", "a")
        print(">>> Producing {} contact matrices for dataset {}\n".format(len(SCmatrixCollated), datasetName))
        plotsSingleContactProbabilityMatrix(
            SCmatrixCollated, uniqueBarcodes, runName, ListData[datasetName], p,fileNameMD=fileNameMD, datasetName=datasetName,
        )

        # [combines matrices from different embryos and calculates integrated contact probability matrix]
        writeString2File(fileNameMD, "## Ensemble contact probability", "a")
        print(">>> Producing ensemble contact matrix for dataset {}\n".format(datasetName))
        plotsEnsembleContactProbabilityMatrix(
            SCmatrixCollated, uniqueBarcodes, runName, ListData[datasetName], p,fileNameMD=fileNameMD, datasetName=datasetName,
        )

        # [deletes variables before starting new iteration]
        # del SCmatrixCollated, uniqueBarcodes, buildsPWDmatrixCollated, runName
        print("\nDone with dataset {}".format(datasetName))

print("Finished execution")

#%%

# SCmatrix_wt_normalized=normalizeMatrix(SCmatrix)

# cScale,cMin = SCmatrix_wt_normalized.max(), SCmatrix_wt_normalized.min()
# plotMatrix(SCmatrix_wt_normalized,uniqueBarcodes, pixelSize,cm='terrain',clim=cScale, figtitle='HiM counts',cmtitle='probability',nCells=nCells,cMin=cMin)
