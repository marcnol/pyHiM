#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 09:04:10 2020

@author: marcnol
 Produces 
"""


#%% imports and plotting settings
import os
import numpy as np
import argparse

# import matplotlib as plt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json, csv

# import scaleogram as scg

from matrixOperations.HIMmatrixOperations import plotDistanceHistograms, plotMatrix, listsSCtoKeep
from matrixOperations.HIMmatrixOperations import (
    analysisHiMmatrix,
    normalizeMatrix,
    shuffleMatrix,
    plotScalogram,
    calculatesEnsemblePWDmatrix,
)

#%% define and loads datasets


def parseArguments():
    # [parsing arguments]
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with dataset")
    parser.add_argument("-O", "--outputFolder", help="Folder for outputs")

    parser.add_argument(
        "-P", "--parameters", help="Provide name of parameter files. folders2Load.json assumed as default",
    )
    parser.add_argument("-A", "--label", help="Add name of label (e.g. doc)")
    parser.add_argument("-W", "--action", help="Select: [all], [labeled] or [unlabeled] cells plotted ")
    parser.add_argument("--fontsize", help="Size of fonts to be used in matrix")
    parser.add_argument("--axisLabel", help="Use if you want a label in x and y", action="store_true")
    parser.add_argument("--axisTicks", help="Use if you want axes ticks", action="store_true")
    parser.add_argument("--barcodes", help="Use if you want barcode images to be displayed", action="store_true")
    parser.add_argument(
        "--scalingParameter", help="Normalizing scaling parameter of colormap. Max will matrix.max()/scalingParameter"
    )
    parser.add_argument("--cScale", help="Colormap absolute scale")
    parser.add_argument("--plottingFileExtension", help="By default: svg. Other options: pdf, png")
    parser.add_argument(
        "--shuffle",
        help="Provide shuffle vector: 0,1,2,3... of the same size or smaller than the original matrix. No spaces! comma-separated!",
    )
    parser.add_argument("--scalogram", help="Use if you want scalogram image to be displayed", action="store_true")
    parser.add_argument("--inputMatrix", help="contact, PWD, or iPWD. Default is contact")
    parser.add_argument("--pixelSize", help="pixel size in um")
    parser.add_argument("--cmap", help="Colormap. Default: coolwarm")
    parser.add_argument(
        "--PWDmode", help="Mode used to calculate the mean distance. Can be either 'median' or KDE. Default: median"
    )

    args = parser.parse_args()

    runParameters = {}

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = "."
        # rootFolder='/home/marcnol/data'+os.sep+'Experiment_18'

    if args.outputFolder:
        outputFolder = args.outputFolder
    else:
        outputFolder = "none"

    if args.parameters:
        runParameters["parametersFileName"] = args.parameters
    else:
        runParameters["parametersFileName"] = "folders2Load.json"

    if args.label:
        runParameters["label"] = args.label
    else:
        runParameters["label"] = "doc"

    if args.action:
        runParameters["action"] = args.action
    else:
        runParameters["action"] = "labeled"

    if args.fontsize:
        runParameters["fontsize"] = args.fontsize
    else:
        runParameters["fontsize"] = 12

    if args.axisLabel:
        runParameters["axisLabel"] = args.axisLabel
    else:
        runParameters["axisLabel"] = False

    if args.axisTicks:
        runParameters["axisTicks"] = args.axisTicks
    else:
        runParameters["axisTicks"] = False

    if args.barcodes:
        runParameters["barcodes"] = args.barcodes
    else:
        runParameters["barcodes"] = False

    if args.scalingParameter:
        runParameters["scalingParameter"] = float(args.scalingParameter)
    else:
        runParameters["scalingParameter"] = 1.0

    if args.cScale:
        runParameters["cScale"] = float(args.cScale)
    else:
        runParameters["cScale"] = 0.0

    if args.plottingFileExtension:
        runParameters["plottingFileExtension"] = "." + args.plottingFileExtension
    else:
        runParameters["plottingFileExtension"] = ".png"

    if args.shuffle:
        runParameters["shuffle"] = args.shuffle
    else:
        runParameters["shuffle"] = 0

    if args.scalogram:
        runParameters["scalogram"] = args.scalogram
    else:
        runParameters["scalogram"] = False

    if args.inputMatrix:
        runParameters["inputMatrix"] = args.inputMatrix
    else:
        runParameters["inputMatrix"] = "contact"

    if args.pixelSize:
        runParameters["pixelSize"] = args.pixelSize
    else:
        runParameters["pixelSize"] = 0.1

    if args.cmap:
        runParameters["cmap"] = args.cmap
    else:
        runParameters["cmap"] = "coolwarm"

    if args.PWDmode:
        runParameters["PWDmode"] = args.PWDmode
    else:
        runParameters["PWDmode"] = "median"

    return rootFolder, outputFolder, runParameters


#%%

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    print(">>> Producing HiM matrix")
    rootFolder, outputFolder, runParameters = parseArguments()

    HiMdata = analysisHiMmatrix(runParameters, rootFolder)

    HiMdata.loadData()

    nCells = HiMdata.nCellsLoaded()

    HiMdata.retrieveSCmatrix()

    if runParameters["inputMatrix"] == "contact":
        # contact probability matrix
        matrix = HiMdata.data["ensembleContactProbability"]
        cScale = matrix.max() / runParameters["scalingParameter"]
    elif runParameters["inputMatrix"] == "PWD":
        # PWD matrix
        # SCmatrix = HiMdata.SCmatrixSelected
        # print("SC matrix size: {}".format(SCmatrix.shape))
        # matrix, keepPlotting = calculatesEnsemblePWDmatrix(SCmatrix, runParameters["pixelSize"],mode=runParameters["PWDmode"])
        SCmatrix = HiMdata.SCmatrixSelected
        cells2Plot = listsSCtoKeep(runParameters, HiMdata.data["SClabeledCollated"])
        print("N cells to plot: {}/{}".format(len(cells2Plot), SCmatrix.shape[2]))
        matrix, keepPlotting = calculatesEnsemblePWDmatrix(
            SCmatrix, runParameters["pixelSize"], cells2Plot, mode=runParameters["PWDmode"]
        )
        if runParameters["cScale"] == 0:
            cScale = matrix[~np.isnan(matrix)].max() / runParameters["scalingParameter"]
        else:
            cScale = runParameters["cScale"]
    elif runParameters["inputMatrix"] == "iPWD":
        SCmatrix = HiMdata.SCmatrixSelected
        cells2Plot = listsSCtoKeep(runParameters, HiMdata.data["SClabeledCollated"])
        print("N cells to plot: {}/{}".format(len(cells2Plot), SCmatrix.shape[2]))
        matrix, keepPlotting = calculatesEnsemblePWDmatrix(
            SCmatrix, runParameters["pixelSize"], cells2Plot, mode=runParameters["PWDmode"]
        )
        matrix = np.reciprocal(matrix)
        cScale = runParameters["cScale"]

    print("scalingParameters, scale={}, {}".format(runParameters["scalingParameter"], cScale))

    nCells = HiMdata.nCellsLoaded()

    nDatasets = len(HiMdata.data["runName"])

    if outputFolder == "none":
        outputFolder = HiMdata.dataFolder

    outputFileName = (
        outputFolder
        + os.sep
        + "Fig_HiMmatrix"
        + "_dataset1:"
        + HiMdata.datasetName
        + "_label:"
        + runParameters["label"]
        + "_action:"
        + runParameters["action"]
        + runParameters["plottingFileExtension"]
    )

    if runParameters["barcodes"]:
        fig1 = plt.figure(figsize=(10, 10), constrained_layout=False)
        gs1 = fig1.add_gridspec(nrows=19, ncols=22, left=0.05, right=0.95, wspace=0.05, hspace=0.05)
        f1 = fig1.add_subplot(gs1[0:-1, 5:-1])
        f2 = fig1.add_subplot(gs1[:-1, 3], sharey=f1)
        f3 = fig1.add_subplot(gs1[-1, 5:-1], sharex=f1)
        ATACseqMatrix = np.array(HiMdata.ListData[HiMdata.datasetName]["BarcodeColormap"]) / 10
        ATACseqMatrixV = np.copy(ATACseqMatrix).reshape((-1, 1))
        pos1 = f2.imshow(np.atleast_2d(ATACseqMatrixV), cmap="tab10")  # colormaps RdBu seismic
        f2.set_xticklabels(())
        f2.set_yticklabels(())
        pos1.set_clim(vmin=-1, vmax=1)

        pos2 = f3.imshow(np.atleast_2d(ATACseqMatrix), cmap="tab10")  # colormaps RdBu seismic
        f3.set_xticklabels(())
        f3.set_yticklabels(())
        pos2.set_clim(vmin=-1, vmax=1)

        barcodeLabels = np.arange(1, ATACseqMatrix.shape[0] + 1)
        for j in range(len(ATACseqMatrix)):
            text = f3.text(
                j,
                0,
                barcodeLabels[j],
                ha="center",
                va="center",
                color="w",
                fontsize=int((14.0 / 22.0) * float(runParameters["fontsize"])),
            )
            text = f2.text(
                0,
                j,
                barcodeLabels[j],
                ha="center",
                va="center",
                color="w",
                fontsize=int((14.0 / 22.0) * float(runParameters["fontsize"])),
            )

        colorbar = False
    else:
        fig1 = plt.figure(constrained_layout=True)
        spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
        f1 = fig1.add_subplot(spec1[0, 0])  # 16
        colorbar = True

    if runParameters["shuffle"] == 0:
        index = range(matrix.shape[0])
    else:
        index = [int(i) for i in runParameters["shuffle"].split(",")]
        matrix = shuffleMatrix(matrix, index)

    f1_ax1_im = HiMdata.plot2DMatrixSimple(
        f1,
        matrix,
        list(HiMdata.data["uniqueBarcodes"]),
        runParameters["axisLabel"],
        runParameters["axisLabel"],
        cmtitle="probability",
        cMin=0,
        cMax=cScale,
        cm=runParameters["cmap"],
        fontsize=runParameters["fontsize"],
        colorbar=colorbar,
        axisTicks=runParameters["axisTicks"],
        nCells=nCells,
        nDatasets=nDatasets,
        showTitle=True,
    )

    # HiMdata.update_clims(0, cScale, f1)
    print("Output written to {}".format(outputFileName))
    plt.savefig(outputFileName)
    np.save(outputFileName + ".npy", matrix)
    titleText = "N = {} | n = {}".format(nCells, nDatasets)
    print("Title: {}".format(titleText))
    print("Output figure: {}".format(outputFileName))

    # if runParameters["scalogram"]:
    #     outputFileNameScalogram = (
    #         outputFolder
    #         + os.sep
    #         + "Fig_HiMmatrix_scalogram"
    #         + "_dataset1:"
    #         + HiMdata.datasetName
    #         + "_label:"
    #         + runParameters["label"]
    #         + "_action:"
    #         + runParameters["action"]
    #         + runParameters["plottingFileExtension"]
    #     )

    #     plotScalogram(matrix, outputFileNameScalogram)

    print("\nDone\n\n")
