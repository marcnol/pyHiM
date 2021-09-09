#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 17:59:34 2020

@author: marcnol

plots 4M profiles given a list of anchors.

Can work with up to two datasets


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

from matrixOperations.HIMmatrixOperations import analysisHiMmatrix, plot1Dprofile2Datasets

#%% define and loads datasets


def parseArguments():
    # [parsing arguments]
    parser = argparse.ArgumentParser()
    parser.add_argument("-F1", "--rootFolder1", help="Folder with dataset 1")
    parser.add_argument("-F2", "--rootFolder2", help="Folder with dataset 2")
    parser.add_argument("-O", "--outputFolder", help="Folder for outputs")

    parser.add_argument(
        "-P", "--parameters", help="Provide name of parameter files. folders2Load.json assumed as default",
    )
    parser.add_argument("-A1", "--label1", help="Add name of label for dataset 1 (e.g. doc)")
    parser.add_argument("-W1", "--action1", help="Select: [all], [labeled] or [unlabeled] cells plotted for dataset 1 ")
    parser.add_argument("-A2", "--label2", help="Add name of label for dataset 1  (e.g. doc)")
    parser.add_argument("-W2", "--action2", help="Select: [all], [labeled] or [unlabeled] cells plotted for dataset 1 ")
    parser.add_argument("--fontsize", help="Size of fonts to be used in matrix")
    parser.add_argument("--axisLabel", help="Use if you want a label in x and y", action="store_true")
    parser.add_argument("--axisTicks", help="Use if you want axes ticks", action="store_true")
    parser.add_argument("--splines", help="Use if you want plot data using spline interpolations", action="store_true")
    parser.add_argument("--cAxis", help="absolute cAxis value for colormap")
    parser.add_argument("--plottingFileExtension", help="By default: svg. Other options: pdf, png")
    parser.add_argument("--legend", help="Use if you want to show legends", action="store_true")
    parser.add_argument(
        "--normalize", help="Matrix normalization factor: maximum, none, single value (normalize 2nd profile by)"
    )

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

    if args.axisLabel:
        runParameters["axisLabel"] = args.axisLabel
    else:
        runParameters["axisLabel"] = False

    if args.axisTicks:
        runParameters["axisTicks"] = args.axisTicks
    else:
        runParameters["axisTicks"] = True

    if args.splines:
        runParameters["splines"] = args.splines
    else:
        runParameters["splines"] = False

    if args.cAxis:
        runParameters["cAxis"] = float(args.cAxis)
    else:
        runParameters["cAxis"] = 0.8

    if args.plottingFileExtension:
        runParameters["plottingFileExtension"] = "." + args.plottingFileExtension
    else:
        runParameters["plottingFileExtension"] = ".svg"

    if args.legend:
        runParameters["legend"] = args.legend
    else:
        runParameters["legend"] = False

    if args.normalize:
        runParameters["normalize"] = args.normalize
    else:
        runParameters["normalize"] = "none"

    print("Input Folders:{}, {}".format(rootFolder1, rootFolder2))
    print("Input parameters:{}".format(runParameters))

    return rootFolder1, rootFolder2, outputFolder, runParameters


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    run2Datasets = False

    rootFolder1, rootFolder2, outputFolder, runParameters = parseArguments()
    print("RootFolders: \n{}\n{}".format(rootFolder1, rootFolder2))
    HiMdata1 = analysisHiMmatrix(runParameters, rootFolder1)
    HiMdata1.runParameters["action"] = HiMdata1.runParameters["action1"]
    HiMdata1.runParameters["label"] = HiMdata1.runParameters["label1"]
    HiMdata1.loadData()
    nCells = HiMdata1.nCellsLoaded()

    if outputFolder == "none":
        outputFolder = HiMdata1.dataFolder

    outputFileName = (
        outputFolder
        + os.sep
        + "Fig_4Mcontacts"
        + "_dataset1:"
        + HiMdata1.datasetName
        + "_label1:"
        + runParameters["label1"]
        + "_action1:"
        + runParameters["action1"]
    )

    if runParameters["run2Datasets"]:
        HiMdata2 = analysisHiMmatrix(runParameters, rootFolder2)
        HiMdata2.runParameters["action"] = HiMdata2.runParameters["action2"]
        HiMdata2.runParameters["label"] = HiMdata2.runParameters["label2"]
        HiMdata2.loadData()
        nCells2 = HiMdata2.nCellsLoaded()

        run2Datasets = True
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

    anchors = HiMdata1.ListData[HiMdata1.datasetName]["3wayContacts_anchors"]
    print("Anchors: {}".format(anchors))

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
                Yticks.append(True)
            else:
                Yticks.append(False)

    FigLabels = [i for i in list(HiMdata1.dataFiles.keys()) if "anchor" in i]
    legendList = [False] * len(anchors)
    if runParameters["legend"]:
        legendList[0] = True

    for anchor, ifigure, iFigLabel, yticks, xticks, legend in zip(
        anchors, FigList, FigLabels, Yticks, Xticks, legendList
    ):
        if not run2Datasets:
            HiMdata1.plot1Dprofile1Dataset(ifigure, anchor, iFigLabel, yticks, xticks)
        else:
            plot1Dprofile2Datasets(
                ifigure, HiMdata1, HiMdata2, runParameters, anchor, iFigLabel, yticks, xticks, legend
            )

    plt.savefig(outputFileName)
    print("Output figure: {}".format(outputFileName))
