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
import csv

from alignBarcodesMasks import plotDistanceHistograms, plotMatrix, calculateContactProbabilityMatrix

from fileManagement import writeString2File

from HIMmatrixOperations import (
    loadsSCdata,
    plotsEnsemble3wayContactMatrix,
    listsSCtoKeep,
)

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")

# =============================================================================
# FUNCTIONS
# =============================================================================q

def joinsListArrays(ListArrays,axis=0):
    joinedArray=np.zeros(0)
    for iArray in ListArrays:
        if joinedArray.shape[0]==0:
            joinedArray = iArray
        else:
            joinedArray = np.concatenate((joinedArray,iArray) ,axis=axis)
    return joinedArray
    
   

def plotsSinglePWDmatrices(
    SCmatrixCollated, uniqueBarcodes, runName, iListData, p, fileNameMD="tmp.md", datasetName="",
):
    # plots distance matrix for each dataset
    for iSCmatrixCollated, iuniqueBarcodes, iTag, mask in zip(SCmatrixCollated, uniqueBarcodes, runName, p["SClabeledCollated"]):
        outputFileName = p["outputFolder"] + os.sep + iTag + "_Cells:" + p["action"] + "_PWDmatrix"

        cells2Plot = listsSCtoKeep(p, mask)
        # print('Dataset {} cells2plot: {}'.format(iTag,cells2Plot))
        plotMatrix(
            iSCmatrixCollated,
            iuniqueBarcodes,
            p["pixelSize"],
            outputFileName=outputFileName,
            figtitle="PWD:" + iTag,
            cm=iListData["PWD_cm"],
            clim=iListData["PWD_clim"],
            mode=iListData["PWD_mode"],
            nCells=iSCmatrixCollated.shape[2],
            cells2Plot=cells2Plot,
        )  # twilight_shifted_r 1.4, mode: median KDE coolwarm terrain
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def plotsInversePWDmatrice(SCmatrixCollated, uniqueBarcodes, runName, iListData, p, fileNameMD, datasetName=""):
    # plots inverse distance matrix for each dataset
    for iSCmatrixCollated, iuniqueBarcodes, iTag, mask in zip(SCmatrixCollated, uniqueBarcodes, runName, p["SClabeledCollated"]):
        outputFileName = p["outputFolder"] + os.sep + iTag + "_Cells:" + p["action"] + "_invPWDmatrix"

        cells2Plot = listsSCtoKeep(p, mask)
        # print('Dataset {} cells2plot: {}'.format(iTag,cells2Plot))

        plotMatrix(
            iSCmatrixCollated,
            iuniqueBarcodes,
            p["pixelSize"],
            cm=iListData["iPWD_cm"],
            outputFileName=outputFileName,
            clim=iListData["iPWD_clim"],
            mode=iListData["iPWD_mode"],
            figtitle="inverse PWD:" + iTag,
            cmtitle="inverse distance, 1/nm",
            inverseMatrix=True,
            nCells=iSCmatrixCollated.shape[2],
            cells2Plot=cells2Plot,
        )  # twilight_shifted_r, mode: median KDE
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def plotsSingleContactProbabilityMatrix(
    SCmatrixCollated, uniqueBarcodes, runName, iListData, p, fileNameMD="tmp.md", datasetName="",
):
    # Plots contact probability matrices for each dataset

    for iSCmatrixCollated, iuniqueBarcodes, iTag, mask in zip(SCmatrixCollated, uniqueBarcodes, runName, p["SClabeledCollated"]):
        cells2Plot = listsSCtoKeep(p, mask)
        if not cells2Plot:
            break
        
        if max(cells2Plot) > iSCmatrixCollated.shape[2]:
            print(
                "Error with range in cells2plot {} as it is larger than the number of available cells {}".format(
                    max(cells2Plot), iSCmatrixCollated.shape[2]
                )
            )
        else:
            SCmatrix, nCells = calculateContactProbabilityMatrix(
                iSCmatrixCollated[:, :, cells2Plot],
                iuniqueBarcodes,
                p["pixelSize"],
                threshold=iListData["ContactProbability_distanceThreshold"],
                norm="nonNANs",
            )  # norm: nCells (default), nonNANs
            outputFileName = p["outputFolder"] + os.sep + iTag + "_Cells:" + p["action"] + "_contactProbability"

            print("Dataset {} cells2plot: {}".format(iTag, cells2Plot))
            cScale = SCmatrix.max() / iListData["ContactProbability_scale"]

            plotMatrix(
                SCmatrix,
                iuniqueBarcodes,
                p["pixelSize"],
                cm=iListData["ContactProbability_cm"],
                outputFileName=outputFileName,
                cMin=iListData["ContactProbability_cmin"],
                clim=cScale,
                figtitle="HiM:" + iTag,
                cmtitle="probability",
                nCells=nCells,
                cells2Plot=cells2Plot,
            )  # twilight_shifted_r terrain coolwarm
            writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def plotsEnsembleContactProbabilityMatrix(
    SCmatrixCollated, uniqueBarcodes, runName, iListData, p, fileNameMD="tmp.md", datasetName="",
):

    # combines matrices from different embryos and calculates integrated contact probability matrix

    SCmatrixAllDatasets = []  # np.zeros((nBarcodes,nBarcodes))

    for iSCmatrixCollated, iuniqueBarcodes, mask, iTag in zip(SCmatrixCollated, uniqueBarcodes, p["SClabeledCollated"], runName):
        cells2Plot = listsSCtoKeep(p, mask)

        if not cells2Plot:
            break

        if max(cells2Plot) > iSCmatrixCollated.shape[2]:
            print(
                "Error: max in cells2plot {} in dataset {} is larger than the number of available cells {}".format(
                    max(cells2Plot), iTag, iSCmatrixCollated.shape[2]
                )
            )
        else:
            if len(SCmatrixAllDatasets) > 0:
                SCmatrixAllDatasets = np.concatenate((SCmatrixAllDatasets, iSCmatrixCollated[:, :, cells2Plot]), axis=2)
            else:
                SCmatrixAllDatasets = iSCmatrixCollated[:, :, cells2Plot]

            commonSetUniqueBarcodes = iuniqueBarcodes

    print("nCells processed: {}".format(SCmatrixAllDatasets.shape[2]))
    SCmatrix, nCells = calculateContactProbabilityMatrix(
        SCmatrixAllDatasets,
        commonSetUniqueBarcodes,
        p["pixelSize"],
        threshold=iListData["ContactProbability_distanceThreshold"],
        norm="nonNANs",
    )  # norm: nCells (default), nonNANs
    cScale = SCmatrix.max() / iListData["ContactProbability_scale"]
    outputFileName = p["outputFolder"] + os.sep + datasetName + "_Cells:" + p["action"] + "_ensembleContactProbability"
    writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")

    plotMatrix(
        SCmatrix,
        commonSetUniqueBarcodes,
        p["pixelSize"],
        cm=iListData["ContactProbability_cm"],
        outputFileName=outputFileName,
        clim=cScale,
        cMin=iListData["ContactProbability_cmin"],
        figtitle="HiM counts",
        cmtitle="probability",
        nCells=nCells,
    )  # twilight_shifted_r

    np.savetxt(
        p["outputFolder"] + os.sep + "CombinedMatrix" + ":" + list(ListData.keys())[0] + "_Cells:" + p["action"] + ".dat",
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
        p["outputFolder"] + os.sep + "UniqueBarcodes" + ":" + list(ListData.keys())[0] + "_Cells:" + p["action"] + ".dat",
        iuniqueBarcodes,
        fmt="%.4f",
        delimiter=" ",
        newline="\n",
        header="unique barcodes",
        footer="",
        comments="# ",
        encoding=None,
    )

    return SCmatrix, iuniqueBarcodes


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

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

    p = {}
    p["pixelSize"] = 0.1

    if args.parameters:
        p["parametersFileName"] = args.parameters
    else:
        p["parametersFileName"] = "folders2Load.json"

    if args.label:
        p["label"] = args.label
    else:
        p["label"] = "doc"

    if args.action:
        p["action"] = args.action
    else:
        p["action"] = "all"

    # [ initialises MD file]
    now = datetime.now()
    dateTime = now.strftime("%d%m%Y_%H%M%S")
    fileNameRoot = "processHiMmatrixAnalysis_"

    # [ Lists and loads datasets from different embryos]
    fileNameListDataJSON = rootFolder + os.sep + p["parametersFileName"]
    print("\n--------------------------------------------------------------------------")
    if os.path.exists(fileNameListDataJSON):
        with open(fileNameListDataJSON) as json_file:
            ListData = json.load(json_file)
        print("Loaded JSON file with {} datasets from {}\n".format(len(ListData), fileNameListDataJSON))

    # [ creates output folder]
    p["outputFolder"] = rootFolder + os.sep + "scHiMmatrices"
    if not os.path.exists(p["outputFolder"]):
        os.mkdir(p["outputFolder"])
        print("Folder created: {}".format(p["outputFolder"]))

    # [loops over lists of datafolders]
    for datasetName in list(ListData.keys()):

        # [loads SC matrices]
        (SCmatrixCollated, uniqueBarcodes, buildsPWDmatrixCollated, runName, SClabeledCollated,) = loadsSCdata(
            ListData, datasetName, p
        )

        fileNameMD = rootFolder + os.sep + fileNameRoot + "_" + datasetName + "_Cells:" + p["action"] + "_" + dateTime + ".md"
        writeString2File(fileNameMD, "# Post-processing of Hi-M matrices", "w")
        writeString2File(fileNameMD, "**dataset: {}** - **Cells: {}**".format(datasetName, p["action"]), "a")
        p["SClabeledCollated"] = SClabeledCollated

        if len(SCmatrixCollated) > 0:
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
                SCmatrixCollated, uniqueBarcodes, runName, ListData[datasetName], p, fileNameMD, datasetName=datasetName,
            )

            # [Plots contact probability matrices for each dataset]
            writeString2File(fileNameMD, "## single cell Contact Probability matrices", "a")
            print(">>> Producing {} contact matrices for dataset {}\n".format(len(SCmatrixCollated), datasetName))
            plotsSingleContactProbabilityMatrix(
                SCmatrixCollated,
                uniqueBarcodes,
                runName,
                ListData[datasetName],
                p,
                fileNameMD=fileNameMD,
                datasetName=datasetName,
            )

            # [combines matrices from different embryos and calculates integrated contact probability matrix]
            writeString2File(fileNameMD, "## Ensemble contact probability", "a")
            print(">>> Producing ensemble contact matrix for dataset {}\n".format(datasetName))
            SCmatrixCollatedEnsemble, uniqueBarcodes = plotsEnsembleContactProbabilityMatrix(
                SCmatrixCollated,
                uniqueBarcodes,
                runName,
                ListData[datasetName],
                p,
                fileNameMD=fileNameMD,
                datasetName=datasetName,
            )

            anchors=ListData[datasetName]['3wayContacts_anchors']
            anchors = [a - 1 for a in anchors]  # convert to zero-based
            sOut = "Probability"  # Probability or Counts
            writeString2File(fileNameMD, "## Ensemble 3way contacts", "a")
            print(">>> Producing ensemble 3way contact matrix for dataset {}\n".format(datasetName))
            plotsEnsemble3wayContactMatrix(
                SCmatrixCollated,
                uniqueBarcodes,
                anchors,
                sOut,
                runName,
                ListData[datasetName],
                p,
                fileNameMD=fileNameMD,
                datasetName=datasetName,
            )

            # [deletes variables before starting new iteration]
            # del SCmatrixCollated, uniqueBarcodes, buildsPWDmatrixCollated, runName
            print("\nDone with dataset {}".format(datasetName))
        else:
            print("\n Could not load ANY dataset!\n")

        # [saves output files]

        # creates outputFileName root
        outputFileName = p["outputFolder"] + os.sep + datasetName + "_label:" + p["label"] + "_action:" + p["action"]
        
        # saves npy arrays
        np.save(outputFileName + "_ensembleContactProbability.npy", SCmatrixCollatedEnsemble)
        np.save(outputFileName + "_SCmatrixCollated.npy", joinsListArrays(SCmatrixCollated,axis=2))
        np.save(outputFileName + "_SClabeledCollated.npy", joinsListArrays(SClabeledCollated,axis=0))

        # saves lists
        with open(outputFileName + "_uniqueBarcodes.csv", "w", newline="") as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=" ", quotechar="|", quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(uniqueBarcodes)

        p["SClabeledCollated"] = []
        with open(outputFileName + "_parameters.json", "w") as f:
            json.dump(p, f, ensure_ascii=False, sort_keys=True, indent=4)

        with open(outputFileName + "_runName.csv", "w", newline="") as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=" ", quotechar="|", quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(runName)
        
        print("Finished execution")



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    