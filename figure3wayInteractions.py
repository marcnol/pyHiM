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
from alignBarcodesMasks import plotDistanceHistograms, plotMatrix
import scaleogram as scg

from HIMmatrixOperations import plotsEnsemble3wayContactMatrix, calculate3wayContactMatrix, getMultiContact, analysisHiMmatrix

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
    # parser.add_argument("--axisLabel", help="Use if you want a label in x and y", action="store_true")
    # parser.add_argument("--axisTicks", help="Use if you want axes ticks", action="store_true")
    parser.add_argument("--scalingParameter", help="Scaling parameter of colormap")
    parser.add_argument("--colorbar", help="Use if you want a colorbar", action="store_true")
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

    # if args.axisLabel:
    #     runParameters["axisLabel"] = args.axisLabel
    # else:
    #     runParameters["axisLabel"] = False

    # if args.axisTicks:
    #     runParameters["axisTicks"] = args.axisTicks
    # else:
    #     runParameters["axisTicks"] = False

    if args.scalingParameter:
        runParameters["scalingParameter"] = float(args.scalingParameter)
    else:
        runParameters["scalingParameter"] = 1.0

    if args.colorbar:
        runParameters["colorbar"] = args.colorbar
    else:
        runParameters["colorbar"] = False

    return rootFolder, runParameters


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print(">>> Producing HiM 3-way matrices")

    # [parsing arguments

    rootFolder, runParameters = parseArguments()

    HiMdata = analysisHiMmatrix(runParameters, rootFolder)

    HiMdata.loadData()

    # panel D: 3-way interaction matrices

    pixelSize = 0.1
    cMax = HiMdata.data["ensembleContactProbability"].max() / runParameters["scalingParameter"]
    nCells = HiMdata.data["SCmatrixCollated"].shape[2]
    plottingFileExtension = ".svg"
    outputFileName = (
        HiMdata.dataFolder
        + os.sep
        + "Fig_3wayInteractions"
        + "_label:"
        + runParameters["label"]
        + "_action:"
        + runParameters["action"]
        + plottingFileExtension
    )

    fig2 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=3, nrows=2, figure=fig2)
    f2_P3 = fig2.add_subplot(spec2[0, 0])  # 16
    f2_P2 = fig2.add_subplot(spec2[0, 1])  # 10
    f2_P1 = fig2.add_subplot(spec2[0, 2])  # 6
    f2_Ea = fig2.add_subplot(spec2[1, 0])  # 13
    f2_Eb = fig2.add_subplot(spec2[1, 1])  # 9
    f2_Ec = fig2.add_subplot(spec2[1, 2])  # 4

    FigList = [f2_P1, f2_P2, f2_P3, f2_Ea, f2_Eb, f2_Ec]
    FigLabels = [i for i in list(HiMdata.dataFiles.keys()) if "anchor" in i]
    yticks = [False, False, True, True, False, False]
    xticks = [False, False, False, True, True, True]

    for ifigure, iFigLabel, iyticks, ixticks in zip(FigList, FigLabels, yticks, xticks):
        f2_ax1_im = HiMdata.plot2DMatrixSimple(
            ifigure,
            HiMdata.data[iFigLabel],
            list(HiMdata.data["uniqueBarcodes"]),
            iyticks,
            ixticks,
            cmtitle="probability",
            cMin=0,
            cMax=cMax,
            fontsize=12,
        )
    if runParameters["colorbar"]:
        cbar_ax = fig2.add_axes([0.995, 0.25, 0.02, 0.6])
        cbar = fig2.colorbar(f2_ax1_im, cax=cbar_ax, fraction=0.046, pad=0.04)
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs, fontsize=12)
    # update_clims(0, cMax, FigList)

    plt.savefig(outputFileName)
    print("Output figure: {}".format(outputFileName))

    print("\nDone\n\n")
