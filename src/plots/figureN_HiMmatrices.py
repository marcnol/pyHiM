#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 17:01:30 2020

plots N Hi-M matrices in a subplot
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
from matrixOperations.HIMmatrixOperations import calculateContactProbabilityMatrix

# import scaleogram as scg

from matrixOperations.HIMmatrixOperations import (
    analysisHiMmatrix,
    normalizeMatrix,
    shuffleMatrix,
    plotScalogram,
    listsSCtoKeep,
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
    parser.add_argument("--scalingParameter", help="Scaling parameter of colormap")
    parser.add_argument("--plottingFileExtension", help="By default: svg. Other options: pdf, png")
    parser.add_argument(
        "--shuffle",
        help="Provide shuffle vector: 0,1,2,3... of the same size or smaller than the original matrix. No spaces! comma-separated!",
    )
    parser.add_argument("--scalogram", help="Use if you want scalogram image to be displayed", action="store_true")
    parser.add_argument("--type", help="Provide one of the following: PWD, contact, iPWD")
    parser.add_argument("--pixelSize", help="Provide pixelSize in um")
    parser.add_argument("--cAxis", help="absolute cAxis value for colormap")
    parser.add_argument("--ratio", help="Does ratio between matrices. Default: difference", action="store_true")
    parser.add_argument("--normalizeMatrix", help="Normalizes matrices by maximum. Default: True", action="store_true")

    args = parser.parse_args()

    runParameters = {}

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        # rootFolder = "."
        # rootFolder='/home/marcnol/data'+os.sep+'Experiment_18'
        rootFolder = "/mnt/PALM_dataserv/DATA/gurgo/Quarantaine/Analysis_embryos_cycle_14_16_2020/mixed_embryos_data/26_06_2020_analysis_T=2µm/plotSegments"

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
        runParameters["label"] = "M"

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

    if args.barcodes:
        runParameters["barcodes"] = args.barcodes
    else:
        runParameters["barcodes"] = False

    if args.scalingParameter:
        runParameters["scalingParameter"] = float(args.scalingParameter)
    else:
        runParameters["scalingParameter"] = 1.0

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

    if args.type:
        runParameters["type"] = args.type
    else:
        runParameters["type"] = "contact"

    if args.pixelSize:
        runParameters["pixelSize"] = float(args.pixelSize)
    else:
        runParameters["pixelSize"] = 1

    if args.cAxis:
        runParameters["cAxis"] = [float(i) for i in args.cAxis.split(",")]
    else:
        runParameters["cAxis"] = [0.1]

    if args.ratio:
        runParameters["ratio"] = args.ratio
    else:
        runParameters["ratio"] = False

    if args.normalizeMatrix:
        runParameters["normalizeMatrix"] = args.normalizeMatrix
    else:
        runParameters["normalizeMatrix"] = False

    runParameters["outputFolder"] = outputFolder
    runParameters["rootFolder"] = rootFolder

    return runParameters


# =============================================================================
# FUNCTIONS
# =============================================================================
def plotTADs(ListData, runParameters):

    if len(runParameters["cAxis"]) == 2:
        cScale = runParameters["cAxis"][1]
    else:
        cScale = runParameters["cAxis"][0]
    print("--------\nClim used: {}\n--------------\n".format(cScale))
    fontsize = runParameters["fontsize"]

    for idataSet in dataSets:

        if runParameters["type"] == "contact" and "TAD2plot" in ListData[idataSet].keys():

            Samples = ListData[idataSet]["Folders"]
            cm = ListData[idataSet]["ContactProbability_cm"]

            TAD2plot = ListData[idataSet]["TAD2plot"]
            segmentLabels = ListData[idataSet]["segmentLabels"]
            segment2plot = ListData[idataSet]["segment2plot"]
            Nplots = len(Samples)

            HiMdata = analysisHiMmatrix(runParameters, os.path.dirname(Samples[segment2plot]))
            HiMdata.loadData()

            m1 = HiMdata.data["ensembleContactProbability"]
            if runParameters["normalizeMatrix"]:
                m1 = m1 / m1.max()

            submatrixReference = m1[TAD2plot[0] : TAD2plot[1], TAD2plot[0] : TAD2plot[1]]

            numberBarcodes = HiMdata.data["ensembleContactProbability"].shape[0]
            numberSegments = len(Samples)

            matrixSegmentAnchor = np.zeros((numberBarcodes, numberSegments))

            fig3 = plt.figure(constrained_layout=False, figsize=(5 * Nplots, 5), dpi=300, facecolor="w", edgecolor="k")
            nCols, nRows = Nplots, 1
            spec2 = gridspec.GridSpec(ncols=nCols, nrows=nRows, figure=fig3)

            FigList, Yticks, Xticks = [], [], []
            for iRow in range(nRows):
                for iCol in range(nCols):
                    FigList.append(fig3.add_subplot(spec2[iRow, iCol]))
                    if iRow == nRows - 1:
                        Xticks.append(False)
                    else:
                        Xticks.append(False)
                    if iCol == 0:
                        Yticks.append(True)
                    else:
                        Yticks.append(False)

            FigLabels = [isample.split(os.sep)[-2] for isample in Samples]
            legendList = [False] * len(Samples)
            colorbar = [False] * len(Samples)
            i = 0
            for isample, ifigure, iFigLabel, yticks, xticks, legend, icolorbar in zip(
                Samples, FigList, FigLabels, Yticks, Xticks, legendList, colorbar
            ):

                HiMdata = analysisHiMmatrix(runParameters, os.path.dirname(isample))
                HiMdata.loadData()

                subMatrix = HiMdata.data["ensembleContactProbability"][
                    TAD2plot[0] : TAD2plot[1], TAD2plot[0] : TAD2plot[1]
                ]

                if runParameters["normalizeMatrix"]:
                    subMatrix = subMatrix / subMatrix.max()

                if "ContactProbability_cm" in ListData[idataSet].keys():
                    colormap = ListData[idataSet]["ContactProbability_cm"]

                if runParameters["ratio"] == True:
                    subMatrixNormalized = np.log(submatrixReference / subMatrix)
                    cmtitle = "log(ratio)"
                else:
                    subMatrixNormalized = submatrixReference - subMatrix
                    cmtitle = "difference"

                print("scalingParameters, scale={}, {}".format(runParameters["scalingParameter"], cScale))

                nCells = HiMdata.nCellsLoaded()

                nDatasets = len(HiMdata.data["runName"])

                f2_ax1_im = HiMdata.plot2DMatrixSimple(
                    ifigure,
                    subMatrixNormalized,
                    list(HiMdata.data["uniqueBarcodes"]),
                    runParameters["axisLabel"],
                    runParameters["axisLabel"],
                    cmtitle=segmentLabels[i],
                    cMin=-cScale,
                    cMax=cScale,
                    fontsize=runParameters["fontsize"],
                    colorbar=icolorbar,
                    axisTicks=runParameters["axisTicks"],
                    nCells=nCells,
                    nDatasets=nDatasets,
                    showTitle=True,
                    figTitle=iFigLabel,
                    cm=colormap,
                )

                del HiMdata, subMatrixNormalized
                i += 1
                f2_ax1_im.set_clim(vmin=-cScale, vmax=cScale)
                # print("\n\n======--==--=={}\n\n======--==--==".format(cScale))

            # colorbar=True
            # cbar = fig3.colorbar(f2_ax1_im, ax=ifigure, fraction=0.046, pad=0.04)
            # cbar.minorticks_on()
            # cbar.set_label("difference",fontsize=float(fontsize)*0.85)
            # f2_ax1_im.set_clim(vmin=-cScale, vmax=cScale)

            outputFileName2 = runParameters["outputFileName"].replace("Fig_HiMmatrix", "Fig_TAD")
            print("Output written to {}".format(outputFileName2))
            fig3.savefig(outputFileName2)


def makesplotHiMLineProfile(matrixSegmentAnchor, uniqueBarcodes, segmentLabels, cScale=0.3, cm="RdBu", fontsize=8):

    numberSegments = matrixSegmentAnchor.shape[1]

    fig1 = plt.figure(constrained_layout=True)
    spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
    f1 = fig1.add_subplot(spec1[0, 0])  # 16
    pos = f1.imshow(matrixSegmentAnchor, cmap=cm)
    f1.set_ylabel("barcode #", fontsize=1.5 * float(fontsize))
    f1.set_xlabel("segment ID", fontsize=1.5 * float(fontsize))

    f1.set_xticks(np.arange(numberSegments))
    f1.set_xticklabels(segmentLabels, fontsize=fontsize)
    f1.set_yticks(np.arange(len(uniqueBarcodes)))
    f1.set_yticklabels(uniqueBarcodes, fontsize=fontsize)

    colorbar = True
    cbar = fig1.colorbar(pos, ax=f1, fraction=0.046, pad=0.04)
    cbar.minorticks_on()
    # cbar.set_label("difference",fontsize=float(fontsize)*0.85)
    pos.set_clim(vmin=-cScale, vmax=cScale)

    return fig1


def plotHiMLineProfile(ListData, runParameters):
    if len(runParameters["cAxis"]) == 2:
        cScale = runParameters["cAxis"][1]
    else:
        cScale = runParameters["cAxis"][0]
    print("--------\nClim used: {}\n--------------\n".format(cScale))
    fontsize = runParameters["fontsize"]

    for idataSet in dataSets:
        if runParameters["type"] == "contact" and "plotSegment_anchor" in ListData[idataSet].keys():

            Samples = ListData[idataSet]["Folders"]

            plotSegment_anchor = ListData[idataSet]["plotSegment_anchor"]
            segmentLabels = ListData[idataSet]["segmentLabels"]
            segment2plot = ListData[idataSet]["segment2plot"]
            cm = ListData[idataSet]["ContactProbability_cm"]

            HiMdata = analysisHiMmatrix(runParameters, os.path.dirname(Samples[segment2plot]))
            HiMdata.loadData()
            # m1=HiMdata.data["ensembleContactProbability"]
            m1, _ = calculateContactProbabilityMatrix(
                HiMdata.data["SCmatrixCollated"],
                list(HiMdata.data["uniqueBarcodes"]),
                pixelSize=runParameters["pixelSize"],
                threshold=0.25,
                norm="nonNANs",
            )

            if runParameters["normalizeMatrix"]:
                m1 = m1 / m1.max()
            contactsAnchor = m1[plotSegment_anchor, :]

            numberBarcodes = HiMdata.data["ensembleContactProbability"].shape[0]
            numberSegments = len(Samples)

            matrixSegmentAnchor = np.zeros((numberBarcodes, numberSegments))

            for iSample, sample in enumerate(Samples):

                HiMdata = analysisHiMmatrix(runParameters, os.path.dirname(sample))
                HiMdata.loadData()

                # matrix=HiMdata.data["ensembleContactProbability"]
                matrix, _ = calculateContactProbabilityMatrix(
                    HiMdata.data["SCmatrixCollated"],
                    list(HiMdata.data["uniqueBarcodes"]),
                    pixelSize=runParameters["pixelSize"],
                    threshold=0.25,
                    norm="nonNANs",
                )
                if runParameters["normalizeMatrix"]:
                    matrix = matrix / matrix.max()

                if runParameters["ratio"] == True:
                    matrixSegmentAnchor[:, iSample] = np.log(contactsAnchor / matrix[plotSegment_anchor, :])
                    cmtitle = "log(ratio)"
                else:
                    matrixSegmentAnchor[:, iSample] = contactsAnchor - matrix[plotSegment_anchor, :]
                    cmtitle = "difference"

            uniqueBarcodes = list(HiMdata.data["uniqueBarcodes"])

            fig1 = makesplotHiMLineProfile(
                matrixSegmentAnchor, uniqueBarcodes, segmentLabels, cScale=cScale, cm=cm, fontsize=fontsize
            )
            outputFileName1 = runParameters["outputFileName"].replace("Fig_HiMmatrix", "Fig_Segment")
            print("Output written to {}".format(outputFileName1))
            fig1.savefig(outputFileName1)

            if "barcodes2plot" in ListData[idataSet].keys():
                barcodes2plot = ListData[idataSet]["barcodes2plot"]
                fig2 = makesplotHiMLineProfile(
                    matrixSegmentAnchor[np.arange(barcodes2plot[0], barcodes2plot[1])],
                    uniqueBarcodes[barcodes2plot[0] : barcodes2plot[1]],
                    segmentLabels,
                    cScale=cScale,
                    cm=cm,
                    fontsize=fontsize,
                )
                outputFileName2 = runParameters["outputFileName"].replace("Fig_HiMmatrix", "Fig_Segment_subMatrix")
                print("Output written to {}".format(outputFileName2))
                fig2.savefig(outputFileName2)


def plotMultipleHiMmatrices(ListData, runParameters):
    for idataSet in dataSets:

        Samples = ListData[idataSet]["Folders"]

        Nplots = len(Samples)

        fig2 = plt.figure(constrained_layout=False, figsize=(5 * Nplots, 5), dpi=300, facecolor="w", edgecolor="k")
        # nCols=np.ceil(len(anchors)/2).astype(int)
        nCols = Nplots
        nRows = 1
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

        FigLabels = [isample.split(os.sep)[-2] for isample in Samples]
        legendList = [False] * len(Samples)
        colorbar = [False] * len(Samples)
        # colorbar[-1]=True

        for isample, ifigure, iFigLabel, yticks, xticks, legend, icolorbar in zip(
            Samples, FigList, FigLabels, Yticks, Xticks, legendList, colorbar
        ):

            HiMdata = analysisHiMmatrix(runParameters, os.path.dirname(isample))
            HiMdata.loadData()

            if runParameters["type"] == "contact":
                matrix = HiMdata.data["ensembleContactProbability"]
                # matrix=normalizeMatrix(matrix)
                cScale = matrix.max() / runParameters["scalingParameter"]
                if "ContactProbability_cm" in ListData[idataSet].keys():
                    colormap = ListData[idataSet]["ContactProbability_cm"]

            elif runParameters["type"] == "PWD":
                matrixSC = HiMdata.data["SCmatrixCollated"]
                cells2Plot = listsSCtoKeep(runParameters, HiMdata.data["SClabeledCollated"])
                matrix = runParameters["pixelSize"] * np.nanmedian(matrixSC[:, :, cells2Plot], axis=2)
                cScale = 3 * np.nanmedian(matrix) / runParameters["scalingParameter"]
                if "PWD_cm" in ListData[idataSet].keys():
                    colormap = ListData[idataSet]["PWD_cm"]
                del matrixSC

            elif runParameters["type"] == "iPWD":
                matrixSC = HiMdata.data["SCmatrixCollated"]
                cells2Plot = listsSCtoKeep(runParameters, HiMdata.data["SClabeledCollated"])
                matrixPWD = runParameters["pixelSize"] * np.nanmedian(matrixSC[:, :, cells2Plot], axis=2)
                matrix = np.reciprocal(matrixPWD)
                cScale = 3 * np.reciprocal(np.nanmedian(matrix)) / runParameters["scalingParameter"]
                if "iPWD_cm" in ListData[idataSet].keys():
                    colormap = ListData[idataSet]["iPWD_cm"]
                del matrixPWD, matrixSC

            print("scalingParameters, scale={}, {}".format(runParameters["scalingParameter"], cScale))

            nCells = HiMdata.nCellsLoaded()

            nDatasets = len(HiMdata.data["runName"])

            if runParameters["shuffle"] == 0:
                index = range(matrix.shape[0])
            else:
                index = [int(i) for i in runParameters["shuffle"].split(",")]
                matrix = shuffleMatrix(matrix, index)

            f2_ax1_im = HiMdata.plot2DMatrixSimple(
                ifigure,
                matrix,
                list(HiMdata.data["uniqueBarcodes"]),
                runParameters["axisLabel"],
                runParameters["axisLabel"],
                cmtitle=runParameters["type"],
                cMin=0,
                cMax=cScale,
                fontsize=runParameters["fontsize"],
                colorbar=icolorbar,
                axisTicks=runParameters["axisTicks"],
                nCells=nCells,
                nDatasets=nDatasets,
                showTitle=True,
                figTitle=iFigLabel,
                cm=colormap,
            )

            del HiMdata, matrix

        cbar_ax = fig2.add_axes([0.92, 0.20, 0.005, 0.6])
        cbar = fig2.colorbar(f2_ax1_im, cax=cbar_ax, fraction=0.046, pad=0.04)
        ticklabs = cbar.ax.get_yticklabels()
        ticklabs1 = ["{:04.2f}".format(i * cScale / (len(ticklabs) - 1)) for i in range(len(ticklabs))]
        cbar.ax.set_yticklabels(ticklabs1, fontsize=runParameters["fontsize"])
        cbar.set_label(runParameters["type"], fontsize=1.2 * float(runParameters["fontsize"]))

        # HiMdata.update_clims(0, cScale, f1)
        print("Output written to {}".format(runParameters["outputFileName"]))
        plt.savefig(runParameters["outputFileName"])
        titleText = "N = {} | n = {}".format(nCells, nDatasets)
        print("Title: {}".format(titleText))
        print("Output figure: {}".format(runParameters["outputFileName"]))


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    print(">>> Producing HiM matrix")
    runParameters = parseArguments()

    # loads datasets: parameter files
    fileNameListDataJSON = runParameters["rootFolder"] + os.sep + runParameters["parametersFileName"]
    with open(fileNameListDataJSON) as json_file:
        ListData = json.load(json_file)

    dataSets = list(ListData.keys())
    if runParameters["outputFolder"] == "none":
        runParameters["outputFolder"] = runParameters["rootFolder"]

    runParameters["outputFileName"] = (
        runParameters["outputFolder"]
        + os.sep
        + "Fig_HiMmatrix"
        + "_label:"
        + runParameters["label"]
        + "_action:"
        + runParameters["action"]
        + runParameters["plottingFileExtension"]
    )

    plotMultipleHiMmatrices(ListData, runParameters)

    plotHiMLineProfile(ListData, runParameters)

    plotTADs(ListData, runParameters)

    print("\nDone\n\n")
