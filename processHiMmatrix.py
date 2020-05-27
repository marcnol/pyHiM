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
warnings.filterwarnings('ignore')

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


def loadsSCdata(ListData, datasetName):
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

    SCmatrixCollated, uniqueBarcodes, buildsPWDmatrixCollated, runName = [], [], [], []

    for rootFolder in ListData[datasetName]["Folders"]:
        # [finds sample name]
        runName.append(os.path.basename(os.path.dirname(rootFolder)))

        # [finds loading order]
        files2Process = glob.glob(rootFolder + "/buildsPWDmatrix_*ROI*.ecsv")
        print("{}".format(files2Process))
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
        # fileNameSNDassignedCells='SNDassignedCells.ecsv'
        # SNDassignedCells = Table.read(fileNameSNDassignedCells, format="ascii.ecsv")

        # [loads and accumulates barcodes and scHiM matrix]

        fileNamMatrix = rootFolder + os.sep + "buildsPWDmatrix_HiMscMatrix.npy"
        fileNameBarcodes = rootFolder + os.sep + "buildsPWDmatrix_uniqueBarcodes.ecsv"

        SCmatrix1 = np.load(fileNamMatrix)
        SCmatrixCollated.append(SCmatrix1)
        uniqueBarcodes.append(np.loadtxt(fileNameBarcodes).astype(int))
        buildsPWDmatrixCollated.append(buildsPWDmatrix)
        print("Cells to merge to SCmatrixCollated: {}\n".format(SCmatrix1.shape[2]))
        print("** Merging rootFolder: {}\n".format(rootFolder))

    print("{} datasets loaded".format(len(SCmatrixCollated)))

    return SCmatrixCollated, uniqueBarcodes, buildsPWDmatrixCollated, runName


def plotsSinglePWDmatrices(SCmatrixCollated, uniqueBarcodes, runName, iListData, fileNameMD="tmp.md", datasetName=""):
    # plots distance matrix for each dataset
    for iSCmatrixCollated, iuniqueBarcodes, iTag in zip(SCmatrixCollated, uniqueBarcodes, runName):
        outputFileName = outputFolder + os.sep + iTag + "_PWDmatrix"
        plotMatrix(
            iSCmatrixCollated,
            iuniqueBarcodes,
            pixelSize,
            outputFileName=outputFileName,
            figtitle="PWD:" + iTag,
            cm="terrain",
            clim=iListData['PWD_clim'],
            mode=iListData['PWD_mode'],
            nCells=iSCmatrixCollated.shape[2],
        )  # twilight_shifted_r 1.4, mode: median KDE
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def plotsInversePWDmatrice(SCmatrixCollated, uniqueBarcodes, runName, iListData, fileNameMD, datasetName=""):
    # plots inverse distance matrix for each dataset
    for iSCmatrixCollated, iuniqueBarcodes, iTag in zip(SCmatrixCollated, uniqueBarcodes, runName):
        outputFileName = outputFolder + os.sep + iTag + "_invPWDmatrix"
        plotMatrix(
            iSCmatrixCollated,
            iuniqueBarcodes,
            pixelSize,
            cm="terrain",
            outputFileName=outputFileName,
            clim=iListData['iPWD_clim'],
            mode=iListData['iPWD_mode'],
            figtitle="inverse PWD:" + iTag,
            cmtitle="inverse distance, 1/nm",
            inverseMatrix=True,
            nCells=iSCmatrixCollated.shape[2],
        )  # twilight_shifted_r, mode: median KDE
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def plotsSingleContactProbabilityMatrix(
    SCmatrixCollated, uniqueBarcodes, runName, iListData, fileNameMD="tmp.md", datasetName="",
):
    # Plots contact probability matrices for each dataset

    for iSCmatrixCollated, iuniqueBarcodes, iTag in zip(SCmatrixCollated, uniqueBarcodes, runName):
        SCmatrix, nCells = calculateContactProbabilityMatrix(
            iSCmatrixCollated, iuniqueBarcodes, pixelSize, threshold=iListData['ContactProbability_distanceThreshold'], norm="nonNANs"
        )  # norm: nCells (default), nonNANs
        outputFileName = outputFolder + os.sep + iTag + "_contactProbability"
        cScale = SCmatrix.max() / iListData['ContactProbability_scale']
        plotMatrix(
            SCmatrix,
            iuniqueBarcodes,
            pixelSize,
            cm="terrain",
            outputFileName=outputFileName,
            cMin=iListData['ContactProbability_cmin'],
            clim=cScale,
            figtitle="HiM:" + iTag,
            cmtitle="probability",
            nCells=nCells,
        )  # twilight_shifted_r
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def plotsEnsembleContactProbabilityMatrix(
    SCmatrixCollated, uniqueBarcodes, runName, iListData, fileNameMD="tmp.md", datasetName="",
):

    # combines matrices from different embryos and calculates integrated contact probability matrix

    # nBarcodes = uniqueBarcodes[0].shape[0]
    SCmatrixAllDatasets = []  # np.zeros((nBarcodes,nBarcodes))

    for iSCmatrixCollated, iuniqueBarcodes in zip(SCmatrixCollated, uniqueBarcodes):
        if len(SCmatrixAllDatasets) > 0:
            SCmatrixAllDatasets = np.concatenate((SCmatrixAllDatasets, iSCmatrixCollated), axis=2)
        else:
            SCmatrixAllDatasets = iSCmatrixCollated

        commonSetUniqueBarcodes = iuniqueBarcodes

    SCmatrix, nCells = calculateContactProbabilityMatrix(
        SCmatrixAllDatasets, commonSetUniqueBarcodes, pixelSize, threshold=iListData['ContactProbability_distanceThreshold'], norm="nonNANs",
    )  # norm: nCells (default), nonNANs
    cScale = SCmatrix.max() / iListData['ContactProbability_scale']
    outputFileName = outputFolder + os.sep + datasetName + "_ensembleContactProbability"
    writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")

    plotMatrix(
        SCmatrix,
        iuniqueBarcodes,
        pixelSize,
        cm="terrain",
        outputFileName=outputFileName,
        clim=cScale,
        cMin=iListData['ContactProbability_cmin'],
        figtitle="HiM counts",
        cmtitle="probability",
        nCells=nCells,
    )  # twilight_shifted_r

    np.savetxt(
        outputFolder + os.sep + "CombinedMatrix" + list(ListData.keys())[0] + ".dat",
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
        outputFolder + os.sep + "UniqueBarcodes" + list(ListData.keys())[0] + ".dat",
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

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument(
        "-P", "--parameters", help="Provide name of parameter files. folders2Load.json assumed as default"
    )

    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")
    processingList = {}

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = "."

    if args.parameters:
        parametersFileName = args.parameters
    else:
        parametersFileName = "folders2Load.json"

    pixelSize = 0.1

    # initialises MD file
    now = datetime.now()
    dateTime = now.strftime("%d%m%Y_%H%M%S")
    fileNameRoot = "processHiMmatrixAnalysis_"

    # Lists and loads datasets from different embryos
    fileNameListDataJSON = rootFolder + os.sep + parametersFileName
    if os.path.exists(fileNameListDataJSON):
        with open(fileNameListDataJSON) as json_file:
            ListData = json.load(json_file)
        print("Loaded JSON file with {} datasets from {}\n".format(len(ListData), fileNameListDataJSON))

    # creates output folder
    outputFolder = rootFolder + os.sep + "scHiMmatrices"
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
        print("Folder created: {}".format(outputFolder))

    # [loops over lists of datafolders]
    for datasetName in list(ListData.keys()):

        SCmatrixCollated, uniqueBarcodes, buildsPWDmatrixCollated, runName = loadsSCdata(ListData, datasetName)
        fileNameMD = rootFolder + os.sep + fileNameRoot + "_" + datasetName + "_" + dateTime + ".md"
        writeString2File(fileNameMD, "# Post-processing of Hi-M matrices", "w")
        writeString2File(fileNameMD, "**dataset: {}**".format(datasetName), "a")

        # plots distance matrix for each dataset
        writeString2File(fileNameMD, "## single cell PWD matrices", "a")
        print('>>> Producing {} PWD matrices for dataset {}\n'.format(len(SCmatrixCollated),datasetName))
        plotsSinglePWDmatrices(
            SCmatrixCollated, uniqueBarcodes, runName, ListData[datasetName], fileNameMD, datasetName=datasetName,
        )

        # plots histograms for each dataset
        # for iSCmatrixCollated, iuniqueBarcodes in zip(SCmatrixCollated, uniqueBarcodes):
        #     plotDistanceHistograms(iSCmatrixCollated, pixelSize, mode="KDE", limitNplots=15)

        # plots inverse distance matrix for each dataset
        writeString2File(fileNameMD, "## single cell inverse PWD matrices", "a")
        print('>>> Producing {} inverse PWD matrices for dataset {}\n'.format(len(SCmatrixCollated),datasetName))
        plotsInversePWDmatrice(
            SCmatrixCollated, uniqueBarcodes, runName, ListData[datasetName], fileNameMD, datasetName=datasetName,
        )

        # Plots contact probability matrices for each dataset
        writeString2File(fileNameMD, "## single cell Contact Probability matrices", "a")
        print('>>> Producing {} contact matrices for dataset {}\n'.format(len(SCmatrixCollated),datasetName))
        plotsSingleContactProbabilityMatrix(
            SCmatrixCollated, uniqueBarcodes, runName, ListData[datasetName], fileNameMD=fileNameMD, datasetName=datasetName,
        )

        # combines matrices from different embryos and calculates integrated contact probability matrix
        writeString2File(fileNameMD, "## Ensemble contact probability", "a")
        print('>>> Producing ensemble contact matrix for dataset {}\n'.format(datasetName))
        plotsEnsembleContactProbabilityMatrix(
            SCmatrixCollated, uniqueBarcodes, runName, ListData[datasetName], fileNameMD=fileNameMD, datasetName=datasetName,
        )

        del SCmatrixCollated, uniqueBarcodes, buildsPWDmatrixCollated, runName
        print('Done with dataset {}'.format(datasetName))

print('Finished execution')

    #%%

    # SCmatrix_wt_normalized=normalizeMatrix(SCmatrix)

    # cScale,cMin = SCmatrix_wt_normalized.max(), SCmatrix_wt_normalized.min()
    # plotMatrix(SCmatrix_wt_normalized,uniqueBarcodes, pixelSize,cm='terrain',clim=cScale, figtitle='HiM counts',cmtitle='probability',nCells=nCells,cMin=cMin)



   #%% links matrices to SND masks

    # ListData["wt_docTAD"] = [
    #     "/mnt/disk2/marcnol/data/Experiment_19/026_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_19/009_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_19/006_Embryo/buildsPWDmatrix",
    # ]

    # ListData[
    #     "wt_HresDocLocus"
    # ] = [  #'/mnt/disk2/marcnol/data/Experiment_3/019_Embryo/buildsPWDmatrix',\
    #     "/mnt/disk2/marcnol/data/Experiment_3/007_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_3/016_Embryo/buildsPWDmatrix",  #'/mnt/disk2/marcnol/data/Experiment_3/000_Embryo/buildsPWDmatrix'\
    #     "/mnt/disk2/marcnol/data/Experiment_3/001_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_3/002_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_3/003_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_3/004_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_3/005_Embryo/buildsPWDmatrix",
    # ]
    # ListData["zld_docTAD"] = [
    #     "/mnt/disk2/marcnol/data/Experiment_5/0_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_5/1_Embryo/buildsPWDmatrix",  # '/mnt/disk2/marcnol/data/Experiment_5/2_Embryo/buildsPWDmatrix',\
    #     "/mnt/disk2/marcnol/data/Experiment_5/3_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_5/4_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_5/5_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_5/6_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_5/7_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_5/8_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_5/9_Embryo/buildsPWDmatrix",
    # ]

    # ListData["HiRes_snaTAD"] = [
    #     "/mnt/disk2/marcnol/data/Experiment_4/0_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_4/1_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_4/2_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_4/4_Embryo/buildsPWDmatrix",
    #     "/mnt/disk2/marcnol/data/Experiment_4/5_Embryo/buildsPWDmatrix",
    # ]
      