#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 09:24:51 2020

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
from HIMmatrixOperations import plotsEnsemble3wayContactMatrix, calculate3wayContactMatrix, getMultiContact

from HIMmatrixOperations import analysisHiMmatrix

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
    parser.add_argument("--cAxis", help="absolute cAxis value for colormap")

    args = parser.parse_args()

    runParameters = {}
    runParameters["pixelSize"] = 0.1

    if args.rootFolder1:
        rootFolder1 = args.rootFolder1
    else:
        rootFolder1 = "."

    if args.rootFolder2:
        rootFolder2 = args.rootFolder2
    else:
        rootFolder2 = "."

    if args.outputFolder:
        outputFolder = args.outputFolder
    else:
        outputFolder = "."

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
        runParameters["axisTicks"] = False

    if args.cAxis:
        runParameters["cAxis"] = float(args.cAxis)
    else:
        runParameters["cAxis"] = 2.0

    print("Input Folders:{}, {}".format(rootFolder1, rootFolder2))
    print("Input parameters:{}".format(runParameters))

    return rootFolder1, rootFolder2, outputFolder, runParameters


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    rootFolder1, rootFolder2, outputFolder, runParameters = parseArguments()

    HiMdata1 = analysisHiMmatrix(runParameters, rootFolder1)
    HiMdata1.runParameters["action"] = HiMdata1.runParameters["action1"]
    HiMdata1.runParameters["label"] = HiMdata1.runParameters["label1"]
    HiMdata1.loadData()

    HiMdata2 = analysisHiMmatrix(runParameters, rootFolder2)
    HiMdata2.runParameters["action"] = HiMdata2.runParameters["action2"]
    HiMdata2.runParameters["label"] = HiMdata2.runParameters["label2"]
    HiMdata2.loadData()

    # cScale1 = HiMdata1.data['ensembleContactProbability'].max() / runParameters['cAxis']
    # cScale2 = HiMdata2.data['ensembleContactProbability'].max() / runParameters['scalingParameter']
    # print('scalingParameters={}'.format(runParameters["scalingParameter"] ))

    plottingFileExtension = ".svg"
    outputFileName = (
        outputFolder
        + os.sep
        + "Fig_ratio2HiMmatrices"
        + "_dataset1:"
        + HiMdata1.datasetName
        + "_label1:"
        + runParameters["label1"]
        + "_action1:"
        + runParameters["action1"]
        + "_dataset2:"
        + HiMdata2.datasetName
        + "_label2:"
        + runParameters["label2"]
        + "_action2:"
        + runParameters["action2"]
        + plottingFileExtension
    )

    fig1 = plt.figure(constrained_layout=True)
    spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
    f1 = fig1.add_subplot(spec1[0, 0])  # 16

    if HiMdata1.data["ensembleContactProbability"].shape == HiMdata2.data["ensembleContactProbability"].shape:
        matrix = np.log(HiMdata1.data["ensembleContactProbability"] / HiMdata2.data["ensembleContactProbability"])

        # HiMdata1.plot2DMatrixSimple(f1, HiMdata1.data['ensembleContactProbability'],\
        #                    list(HiMdata1.data['uniqueBarcodes']),\
        #                    runParameters['axisLabel'],\
        #                    runParameters['axisLabel'],\
        #                    cmtitle='probability',
        #                    cMin=0, cMax=cScale1,\
        #                    fontsize=runParameters['fontsize'],\
        #                    colorbar=True,\
        #                    axisTicks=runParameters["axisTicks"],\
        #                    cm='coolwarm')

        # HiMdata1.plot2DMatrixSimple(f1, HiMdata2.data['ensembleContactProbability'],\
        #                     list(HiMdata2.data['uniqueBarcodes']),\
        #                     runParameters['axisLabel'],\
        #                     runParameters['axisLabel'],\
        #                     cmtitle='probability',
        #                     cMin=0, cMax=cScale2,\
        #                     fontsize=runParameters['fontsize'],\
        #                     colorbar=True,\
        #                     axisTicks=runParameters["axisTicks"],\
        #                     cm='coolwarm')

        f1_ax1_im = HiMdata1.plot2DMatrixSimple(
            f1,
            matrix,
            list(HiMdata1.data["uniqueBarcodes"]),
            runParameters["axisLabel"],
            runParameters["axisLabel"],
            cmtitle="log(ratio)",
            cMin=-runParameters["cAxis"],
            cMax=runParameters["cAxis"],
            fontsize=runParameters["fontsize"],
            colorbar=True,
            axisTicks=runParameters["axisTicks"],
            cm="RdBu",
        )
        plt.savefig(outputFileName)
        print("Output figure: {}".format(outputFileName))

    else:
        print("Error: matrices do not have the same dimensions!")
