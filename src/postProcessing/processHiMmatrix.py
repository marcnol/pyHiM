#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:36:20 2020

@author: marcnol


This script takes JSON file with folders where datasets are
stored and processes multiple PWD matrices together.

$ processHiMmatrix.py -F rootFolder

outputs

SCmatrixCollated: 3D npy matrix. PWD matrix for single cells. Axes:0-1 barcodes, Axis:2, cellID
uniqueBarcodes: npy array. list of unique barcodes
SClabeledCollated: npy array. binary label indicating if cell is in pattern or not. Axis:0 cellID

"""

# =============================================================================
# IMPORTS
# =============================================================================q

import numpy as np
import os, sys
import json
from datetime import datetime
import argparse
import csv

#Olivier
csv.field_size_limit(sys.maxsize) 

from fileProcessing.fileManagement import writeString2File

from matrixOperations.HIMmatrixOperations import (
    loadsSCdata,
    plotsEnsemble3wayContactMatrix,
    loadsSCdataMATLAB,
    plotsSinglePWDmatrices,
    plotsInversePWDmatrix,
    plotsSingleContactProbabilityMatrix,
    plotsEnsembleContactProbabilityMatrix,
)

# to remove in a future version
import warnings

# warnings.filterwarnings("ignore")

# =============================================================================
# FUNCTIONS
# =============================================================================q


def joinsListArrays(ListArrays, axis=0):
    joinedArray = np.zeros(0)
    for iArray in ListArrays:
        if joinedArray.shape[0] == 0:
            joinedArray = iArray
        else:
            joinedArray = np.concatenate((joinedArray, iArray), axis=axis)
    return joinedArray

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument(
        "-P", "--parameters", help="Provide name of parameter files. folders2Load.json assumed as default",
    )
    parser.add_argument("-A", "--label", help="Add name of label (e.g. doc)")
    parser.add_argument("-W", "--action", help="Select: [all], [labeled] or [unlabeled] cells plotted ")
    parser.add_argument("--matlab", help="Use to load matlab formatted data", action="store_true")
    parser.add_argument("--saveMatrix", help="Use to load matlab formatted data", action="store_true")
    parser.add_argument("--getStructure", help="Use to save ShEc3D PDB structure", action="store_true")
    parser.add_argument("--pixelSize", help="pixelSize in um")
    parser.add_argument("--HiMnormalization", help="Normalization of contact matrix: nonNANs (default) or nCells")
    parser.add_argument("--d3", help="Use to load 3D maps", action="store_true")

    p={}

    args = parser.parse_args()
    if args.rootFolder:
        p["rootFolder"] = args.rootFolder
    else:
        p["rootFolder"] = "."
        # p["rootFolder"] = "/home/marcnol/grey/docPaper_fullDatasets/updatedDatasets/white_wt_docTAD_nc14"

    if args.parameters:
        p["parametersFileName"] = args.parameters
    else:
        p["parametersFileName"] = "folders2Load.json"

    if args.label:
        p["label"] = args.label
    else:
        p["label"] = "M"

    if args.action:
        p["action"] = args.action
    else:
        p["action"] = "all"

    if args.matlab:
        p["format"] = "matlab"
    else:
        p["format"] = "pyHiM"

    if args.saveMatrix:
        p["saveMatrix"] = True
    else:
        p["saveMatrix"] = False

    if args.getStructure:
        p["getStructure"] = True
    else:
        p["getStructure"] = False

    if args.pixelSize:
        p["pixelSize"] = float(args.pixelSize)
    else:
        p["pixelSize"] = 0.1

    if args.HiMnormalization:
        p["HiMnormalization"] = args.HiMnormalization
    else:
        p["HiMnormalization"] = "nonNANs"

    if args.d3:
        p["d3"] = args.d3
    else:
        p["d3"] = False

    return p
# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    # [parsing arguments]
    p=parseArguments()

    # [ initialises MD file]
    now = datetime.now()
    dateTime = now.strftime("%d%m%Y_%H%M%S")
    fileNameRoot = "processHiMmatrixAnalysis_"

    # [ Lists and loads datasets from different embryos]
    fileNameListDataJSON = p["rootFolder"] + os.sep + p["parametersFileName"]
    print("\n--------------------------------------------------------------------------")
    if os.path.exists(fileNameListDataJSON):
        with open(fileNameListDataJSON) as json_file:
            ListData = json.load(json_file)
        print("Loaded JSON file with {} datasets from {}\n".format(len(ListData), fileNameListDataJSON))
    else:
        print("File not found: {}".format(fileNameListDataJSON))
        sys.exit()

    # [ creates output folder]
    p["outputFolder"] = p["rootFolder"] + os.sep + "scHiMmatrices"
    if not os.path.exists(p["outputFolder"]):
        os.mkdir(p["outputFolder"])
        print("Folder created: {}".format(p["outputFolder"]))

    # [loops over lists of datafolders]
    for datasetName in list(ListData.keys()):

        # [loads SC matrices]
        if p["format"] == "pyHiM":
            print(">>> Loading pyHiM-formatted dataset")
            SCmatrixCollated, uniqueBarcodes, buildsPWDmatrixCollated, runName, SClabeledCollated = loadsSCdata(
                ListData, datasetName, p
            )
        elif p["format"] == "matlab":
            print(">>> Loading MATLAB-formatted dataset")
            SCmatrixCollated, uniqueBarcodes, runName, SClabeledCollated = loadsSCdataMATLAB(ListData, datasetName, p)

        fileNameMD = (
            p["rootFolder"] + os.sep + fileNameRoot + "_" + datasetName + "_Cells:" + p["action"] + "_" + dateTime + ".md"
        )
        writeString2File(fileNameMD, "# Post-processing of Hi-M matrices", "w")
        writeString2File(fileNameMD, "**dataset: {}** - **Cells: {}**".format(datasetName, p["action"]), "a")
        p["SClabeledCollated"] = SClabeledCollated

        if len(SCmatrixCollated) > 0:
            # [plots distance matrix for each dataset]
            writeString2File(fileNameMD, "## single dataset PWD matrices", "a")
            print(">>> Producing {} PWD matrices for dataset {}\n".format(len(SCmatrixCollated), datasetName))
            plotsSinglePWDmatrices(
                SCmatrixCollated,
                uniqueBarcodes,
                runName,
                ListData[datasetName],
                p,
                fileNameMD,
                datasetName=datasetName,
            )

            # plots histograms for each dataset
            # for iSCmatrixCollated, iuniqueBarcodes in zip(SCmatrixCollated, uniqueBarcodes):
            #     plotDistanceHistograms(iSCmatrixCollated, pixelSize, mode="KDE", limitNplots=15)

            # [plots inverse distance matrix for each dataset]
            writeString2File(fileNameMD, "## single dataset inverse PWD matrices", "a")
            print(">>> Producing {} inverse PWD matrices for dataset {}\n".format(len(SCmatrixCollated), datasetName))
            plotsInversePWDmatrix(
                SCmatrixCollated,
                uniqueBarcodes,
                runName,
                ListData[datasetName],
                p,
                fileNameMD,
                datasetName=datasetName,
            )

            # [Plots contact probability matrices for each dataset]
            writeString2File(fileNameMD, "## single dataset Contact Probability matrices", "a")
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
            SCmatrixCollatedEnsemble, commonSetUniqueBarcodes = plotsEnsembleContactProbabilityMatrix(
                SCmatrixCollated,
                uniqueBarcodes,
                runName,
                ListData[datasetName],
                p,
                fileNameMD=fileNameMD,
                datasetName=datasetName,
            )

            anchors = ListData[datasetName]["3wayContacts_anchors"]
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
        if 'SCmatrixCollatedEnsemble' in locals():
            np.save(outputFileName + "_ensembleContactProbability.npy", SCmatrixCollatedEnsemble)

        np.save(outputFileName + "_SCmatrixCollated.npy", joinsListArrays(SCmatrixCollated, axis=2))
        np.save(outputFileName + "_SClabeledCollated.npy", joinsListArrays(SClabeledCollated, axis=0))

        # saves lists
        if 'SCmatrixCollatedEnsemble' in locals():
            with open(outputFileName + "_uniqueBarcodes.csv", "w", newline="") as csvfile:
                spamwriter = csv.writer(csvfile, delimiter=" ", quotechar="|", quoting=csv.QUOTE_MINIMAL)
                spamwriter.writerow(commonSetUniqueBarcodes)

        p["SClabeledCollated"] = []
        with open(outputFileName + "_parameters.json", "w") as f:
            json.dump(p, f, ensure_ascii=False, sort_keys=True, indent=4)

        with open(outputFileName + "_runName.csv", "w", newline="") as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=" ", quotechar="|", quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(runName)

        print("Finished execution")
