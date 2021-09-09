#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 09:18:15 2020

@author: marcnol
"""


#%% imports and plotting settings
import os
import numpy as np
import argparse

# import matplotlib as plt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json, csv
from matrixOperations.HIMmatrixOperations import plotDistanceHistograms, plotMatrix

# import scaleogram as scg

from matrixOperations.HIMmatrixOperations import (
    plotsEnsemble3wayContactMatrix,
    calculate3wayContactMatrix,
    getMultiContact,
    analysisHiMmatrix,
)

#%% define and loads datasets
def parseArguments():
    # [parsing arguments]
    parser = argparse.ArgumentParser()
    parser.add_argument("-F1", "--rootFolder1", help="Folder with dataset 1")
    parser.add_argument("-F2", "--rootFolder2", help="Folder with dataset 2")
    parser.add_argument("-O", "--outputFolder", help="Folder for outputs")
    parser.add_argument("-P", "--parameters", help="Provide name of parameter files. folders2Load.json assumed as default")
    parser.add_argument("-P2", "--parameters2", help="Provide name of parameter files for dataset 2. folders2Load.json assumed as default")
    parser.add_argument("-A1", "--label1", help="Add name of label for dataset 1 (e.g. doc)")
    parser.add_argument("-W1", "--action1", help="Select: [all], [labeled] or [unlabeled] cells plotted for dataset 1 ")
    parser.add_argument("-A2", "--label2", help="Add name of label for dataset 1  (e.g. doc)")
    parser.add_argument("-W2", "--action2", help="Select: [all], [labeled] or [unlabeled] cells plotted for dataset 1 ")
    parser.add_argument("--fontsize", help="Size of fonts to be used in matrix")
    # parser.add_argument("--axisLabel", help="Use if you want a label in x and y", action="store_true")
    # parser.add_argument("--axisTicks", help="Use if you want axes ticks", action="store_true")
    parser.add_argument("--scalingParameter", help="Scaling parameter of colormap")
    parser.add_argument("--colorbar", help="Use if you want a colorbar", action="store_true")
    parser.add_argument("--plottingFileExtension", help="By default: svg. Other options: pdf, png")
    parser.add_argument("--normalize", help="Matrices get normalized by their maximum", action="store_true")
    args = parser.parse_args()

    runParameters = {}
    runParameters["pixelSize"] = 0.1

    if args.rootFolder1:
        rootFolder1 = args.rootFolder1
    else:
        rootFolder1 = "."

    if args.rootFolder2:
        rootFolder2 = args.rootFolder2
        runParameters["run2Datasets"] = True
    else:
        rootFolder2 = "."
        runParameters["run2Datasets"] = False

    if args.outputFolder:
        outputFolder = args.outputFolder
    else:
        outputFolder = "none"

    if args.parameters:
        runParameters["parametersFileName"] = args.parameters
    else:
        runParameters["parametersFileName"] = "folders2Load.json"

    if args.parameters2:
        runParameters["parametersFileName2"] = args.parameters2
    else:
        runParameters["parametersFileName2"] = "folders2Load.json"

    if args.label1:
        runParameters["label1"] = args.label1
    else:
        runParameters["label1"] = "doc"

    if args.label2:
        runParameters["label2"] = args.label2
    else:
        runParameters["label2"] = "doc"

    if args.action1:
        runParameters["action1"] = args.action1
    else:
        runParameters["action1"] = "labeled"

    if args.action2:
        runParameters["action2"] = args.action2
    else:
        runParameters["action2"] = "unlabeled"

    if args.fontsize:
        runParameters["fontsize"] = args.fontsize
    else:
        runParameters["fontsize"] = 12

    if args.scalingParameter:
        runParameters["scalingParameter"] = float(args.scalingParameter)
    else:
        runParameters["scalingParameter"] = 1.0

    if args.colorbar:
        runParameters["colorbar"] = args.colorbar
    else:
        runParameters["colorbar"] = False

    if args.plottingFileExtension:
        runParameters["plottingFileExtension"] = "." + args.plottingFileExtension
    else:
        runParameters["plottingFileExtension"] = ".svg"

    if args.normalize:
        runParameters["normalize"] = args.normalize
    else:
        runParameters["normalize"] = False

    return rootFolder1, rootFolder2, outputFolder, runParameters


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print(">>> Producing HiM 3-way matrices")

    # [parsing arguments
    rootFolder1, rootFolder2, outputFolder, runParameters = parseArguments()

    print(">>> Loading first dataset from {}".format(rootFolder1))
    HiMdata1 = analysisHiMmatrix(runParameters, rootFolder1)
    HiMdata1.runParameters["action"] = HiMdata1.runParameters["action1"]
    HiMdata1.runParameters["label"] = HiMdata1.runParameters["label1"]
    HiMdata1.loadData()
    nCells1 = HiMdata1.nCellsLoaded()

    if outputFolder == "none":
        outputFolder = HiMdata1.dataFolder

    outputFileName = (
        outputFolder
        + os.sep
        + "Fig_3wayContacts"
        + "_dataset1:"
        + HiMdata1.datasetName
        + "_label1:"
        + runParameters["label1"]
        + "_action1:"
        + runParameters["action1"]
    )

    if runParameters["run2Datasets"]:
        print(">>> Loading second dataset from {}".format(rootFolder2))
        HiMdata2 = analysisHiMmatrix(runParameters, rootFolder2)
        HiMdata2.runParameters["action"] = HiMdata2.runParameters["action2"]
        HiMdata2.runParameters["label"] = HiMdata2.runParameters["label2"]
        HiMdata2.runParameters["parametersFileName"] = HiMdata2.runParameters["parametersFileName2"]
        HiMdata2.loadData()
        nCells2 = HiMdata2.nCellsLoaded()

        outputFileName = (
            outputFileName
            + "_dataset2:"
            + HiMdata2.datasetName
            + "_label2:"
            + runParameters["label2"]
            + "_action2:"
            + runParameters["action2"]
        )
    outputFileName += runParameters["plottingFileExtension"]

    # 3-way interaction matrices
    pixelSize = 0.1
    cMax = HiMdata1.data["ensembleContactProbability"].max() / runParameters["scalingParameter"]

    anchors = [int(i.split(":")[1]) for i in list(HiMdata1.dataFiles.keys()) if "anchor" in i]
    fig2 = plt.figure(constrained_layout=True)
    nCols = np.ceil(len(anchors) / 2).astype(int)
    nRows = 2
    spec2 = gridspec.GridSpec(ncols=nCols, nrows=nRows, figure=fig2)

    FigList, Yticks, Xticks = [], [], []
    for iRow in range(nRows):
        for iCol in range(nCols):
            FigList.append(fig2.add_subplot(spec2[iRow, iCol]))
            if iRow == nRows - 1:
                Xticks.append(False)
            else:
                Xticks.append(False)
            if iCol == 0:
                Yticks.append(False)
            else:
                Yticks.append(False)

    FigLabels = [i for i in list(HiMdata1.dataFiles.keys()) if "anchor" in i]
    # print("FigList:{}".format(FigLabels))
    legendList = [False] * len(anchors)
    legendList[0] = True

    for ifigure, iFigLabel, iyticks, ixticks in zip(FigList, FigLabels, Xticks, Yticks):
        if runParameters["run2Datasets"]:
            # mixed matrices from 2 datasets
            if runParameters["normalize"]:
                matrix = HiMdata1.data[iFigLabel] / HiMdata1.data[iFigLabel].max()
            else:
                matrix = HiMdata1.data[iFigLabel]

            for i in range(matrix.shape[0]):
                for j in range(0, i):
                    if runParameters["normalize"]:
                        matrix[i, j] = HiMdata2.data[iFigLabel][i, j] / HiMdata2.data[iFigLabel].max()
                    else:
                        matrix[i, j] = HiMdata2.data[iFigLabel][i, j]
        else:
            # only one matrix
            matrix = HiMdata1.data[iFigLabel]

        print("Dataset: {} | cScale= {}-{}".format(iFigLabel, 0, cMax))
        f2_ax1_im = HiMdata1.plot2DMatrixSimple(
            ifigure,
            matrix,
            list(HiMdata1.data["uniqueBarcodes"]),
            iyticks,
            ixticks,
            cmtitle="probability",
            cMin=0,
            cMax=cMax,
            fontsize=12,
        )
        
        matrixOutputFileName=outputFileName + "_" + iFigLabel + ".npy"
        print("Saving Matrix: {}".format(matrixOutputFileName))
        np.save(matrixOutputFileName,matrix)
        
    if runParameters["colorbar"]:
        cbar_ax = fig2.add_axes([0.995, 0.25, 0.02, 0.6])
        cbar = fig2.colorbar(f2_ax1_im, cax=cbar_ax, fraction=0.046, pad=0.04)
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs, fontsize=12)

    # for ifigure in FigList:
    #     HiMdata1.update_clims(0, cMax, ifigure)

    # update_clims(0, cMax, FigList)


    plt.savefig(outputFileName)
    print("Output figure: {}".format(outputFileName))

    print("\nDone\n\n")
