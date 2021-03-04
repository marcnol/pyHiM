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
from matrixOperations.alignBarcodesMasks import plotDistanceHistograms, plotMatrix

# import scaleogram as scg
from matrixOperations.HIMmatrixOperations import (
    plotsEnsemble3wayContactMatrix,
    calculate3wayContactMatrix,
    getMultiContact,
    shuffleMatrix,
    analysisHiMmatrix,
)


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
    parser.add_argument("--ratio", help="Does ratio between matrices. Default: difference", action="store_true")
    parser.add_argument("--cAxis", help="absolute cAxis value for colormap")
    parser.add_argument("--plottingFileExtension", help="By default: svg. Other options: pdf, png")
    parser.add_argument(
        "--shuffle1",
        help="Provide shuffle vector: 0,1,2,3... of the same size or smaller than the original matrix. No spaces! comma-separated!",
    )
    parser.add_argument(
        "--shuffle2",
        help="Provide shuffle vector: 0,1,2,3... of the same size or smaller than the original matrix. No spaces! comma-separated!",
    )
    parser.add_argument(
        "--cMinMax",
        help="Provide min and max value for the colormap. Comma-separated, no spaces: 0,0.5 Overwrites --cAxis.",
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
        runParameters["label2"] = "NE"

    if args.action1:
        runParameters["action1"] = args.action1
    else:
        runParameters["action1"] = "labeled"

    if args.action2:
        runParameters["action2"] = args.action2
    else:
        runParameters["action2"] = "labeled"

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

    if args.ratio:
        runParameters["ratio"] = args.ratio
    else:
        runParameters["ratio"] = False

    if args.cAxis:
        runParameters["cAxis"] = float(args.cAxis)
    else:
        runParameters["cAxis"] = 0.6

    if args.plottingFileExtension:
        runParameters["plottingFileExtension"] = "." + args.plottingFileExtension
    else:
        runParameters["plottingFileExtension"] = ".svg"

    if args.shuffle1:
        runParameters["shuffle1"] = args.shuffle1
    else:
        runParameters["shuffle1"] = 0

    if args.shuffle2:
        runParameters["shuffle2"] = args.shuffle2
    else:
        runParameters["shuffle2"] = 0

    if args.cMinMax:
        runParameters["cMinMax"] = args.cMinMax
    else:
        runParameters["cMinMax"] = 0

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
    nCells = HiMdata1.nCellsLoaded()

    HiMdata2 = analysisHiMmatrix(runParameters, rootFolder2)
    HiMdata2.runParameters["action"] = HiMdata2.runParameters["action2"]
    HiMdata2.runParameters["label"] = HiMdata2.runParameters["label2"]
    HiMdata2.loadData()
    nCells2 = HiMdata2.nCellsLoaded()

    # cScale1 = HiMdata1.data['ensembleContactProbability'].max() / runParameters['cAxis']
    # cScale2 = HiMdata2.data['ensembleContactProbability'].max() / runParameters['scalingParameter']
    # print('scalingParameters={}'.format(runParameters["scalingParameter"] ))

    if outputFolder == "none":
        outputFolder = HiMdata1.dataFolder

    outputFileName1 = (
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
        + runParameters["plottingFileExtension"]
    )

    outputFileName2 = (
        outputFolder
        + os.sep
        + "Fig_mixedHiMmatrices"
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
        + runParameters["plottingFileExtension"]
    )

    if HiMdata1.data["ensembleContactProbability"].shape == HiMdata2.data["ensembleContactProbability"].shape:
        ### Fig1: difference or ratio of the two matrices
        fig1 = plt.figure(constrained_layout=True)
        spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
        f1 = fig1.add_subplot(spec1[0, 0])  # 16
        m1 = HiMdata1.data["ensembleContactProbability"]
        m2 = HiMdata2.data["ensembleContactProbability"]

        if runParameters["shuffle1"] != 0:
            index1 = [int(i) for i in runParameters["shuffle1"].split(",")]
            m1 = shuffleMatrix(m1, index1)

        if runParameters["shuffle2"] != 0:
            index2 = [int(i) for i in runParameters["shuffle2"].split(",")]
            m2 = shuffleMatrix(m2, index2)

        m1 = m1 / m1.max()
        m2 = m2 / m2.max()

        if runParameters["ratio"] == True:
            matrix = np.log(m1 / m2)
            cmtitle = "log(ratio)"
        else:
            matrix = m1 - m2
            cmtitle = "difference"

        f1_ax1_im = HiMdata1.plot2DMatrixSimple(
            f1,
            matrix,
            list(HiMdata1.data["uniqueBarcodes"]),
            runParameters["axisLabel"],
            runParameters["axisLabel"],
            cmtitle=cmtitle,
            cMin=-runParameters["cAxis"],
            cMax=runParameters["cAxis"],
            fontsize=runParameters["fontsize"],
            colorbar=True,
            axisTicks=runParameters["axisTicks"],
            cm="RdBu",
        )
        plt.savefig(outputFileName1)
        print("Output figure: {}".format(outputFileName1))
        # plt.close()

        ### Fig2: "mixed matrix"
        fig2 = plt.figure(constrained_layout=True)
        spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
        f2 = fig2.add_subplot(spec2[0, 0])  # 16

        # load data once more
        matrix1 = HiMdata1.data["ensembleContactProbability"]
        matrix2 = HiMdata2.data["ensembleContactProbability"]

        if runParameters["shuffle1"] != 0:
            index1 = [int(i) for i in runParameters["shuffle1"].split(",")]
            matrix1 = shuffleMatrix(matrix1, index1)

        if runParameters["shuffle2"] != 0:
            index2 = [int(i) for i in runParameters["shuffle2"].split(",")]
            matrix2 = shuffleMatrix(matrix2, index2)

        for i in range(matrix1.shape[0]):
            for j in range(0, i):
                matrix1[i, j] = matrix2[i, j]

        if runParameters["cMinMax"] == 0:
            cMin = 0
            cMax = runParameters["cAxis"]
        else:
            index = [float(i) for i in runParameters["cMinMax"].split(",")]
            cMin = index[0]
            cMax = index[1]

        HiMdata1.plot2DMatrixSimple(
            f2,
            matrix1,
            list(HiMdata1.data["uniqueBarcodes"]),
            runParameters["axisLabel"],
            runParameters["axisLabel"],
            cmtitle="probability",
            cMin=cMin,
            cMax=cMax,
            fontsize=runParameters["fontsize"],
            colorbar=True,
            axisTicks=runParameters["axisTicks"],
            cm="coolwarm",
        )
        plt.savefig(outputFileName2)
        print("Output figure: {}".format(outputFileName2))

        # save also the npy
        outputFileName3 = outputFileName2.replace(runParameters["plottingFileExtension"], ".npy")
        np.save(outputFileName3, matrix1)
    else:
        print("Error: matrices do not have the same dimensions!")
