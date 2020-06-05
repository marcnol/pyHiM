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
from alignBarcodesMasks import plotDistanceHistograms, plotMatrix
import scaleogram as scg
from HIMmatrixOperations import plotsEnsemble3wayContactMatrix, calculate3wayContactMatrix, getMultiContact

from HIMmatrixOperations import analysisHiMmatrix

#%% define and loads datasets


def parseArguments():
    # [parsing arguments]
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with dataset")
    parser.add_argument(
        "-P", "--parameters", help="Provide name of parameter files. folders2Load.json assumed as default",
    )
    parser.add_argument("-A", "--label", help="Add name of label (e.g. doc)")
    parser.add_argument("-W", "--action", help="Select: [all], [labeled] or [unlabeled] cells plotted ")
    parser.add_argument("--fontsize", help="Size of fonts to be used in matrix")
    parser.add_argument("--axisLabel", help="Use if you want a label in x and y", action="store_true")
    parser.add_argument("--axisTicks", help="Use if you want axes ticks", action="store_true")
    parser.add_argument("--scalingParameter", help="Scaling parameter of colormap")

    args = parser.parse_args()

    runParameters = {}
    runParameters["pixelSize"] = 0.1

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = "."
        # rootFolder='/home/marcnol/data'+os.sep+'Experiment_18'

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
        runParameters["action"] = "all"

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

    if args.scalingParameter:
        runParameters["scalingParameter"] = float(args.scalingParameter)
    else:
        runParameters["scalingParameter"] = 1.0

    return rootFolder, runParameters


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    print(">>> Producing HiM matrix")
    rootFolder, runParameters = parseArguments()

    HiMdata = analysisHiMmatrix(runParameters, rootFolder)

    HiMdata.loadData()

    # panel C: contact probability matrix

    cScale = HiMdata.data["ensembleContactProbability"].max() / runParameters["scalingParameter"]
    print("scalingParameters={}".format(runParameters["scalingParameter"]))
    nCells = HiMdata.data["SCmatrixCollated"].shape[2]
    plottingFileExtension = ".svg"
    outputFileName = (
        HiMdata.dataFolder
        + os.sep
        + "Fig_HiMmatrix"
        + "_label:"
        + runParameters["label"]
        + "_action:"
        + runParameters["action"]
        + plottingFileExtension
    )

    fig1 = plt.figure(constrained_layout=True)
    spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
    f1 = fig1.add_subplot(spec1[0, 0])  # 16

    f1_ax1_im = HiMdata.plot2DMatrixSimple(
        f1,
        HiMdata.data["ensembleContactProbability"],
        list(HiMdata.data["uniqueBarcodes"]),
        runParameters["axisLabel"],
        runParameters["axisLabel"],
        cmtitle="probability",
        cMin=0,
        cMax=cScale,
        fontsize=runParameters["fontsize"],
        colorbar=True,
        axisTicks=runParameters["axisTicks"],
    )
    plt.savefig(outputFileName)
    print("Output figure: {}".format(outputFileName))

    print("\nDone\n\n")
