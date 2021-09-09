#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 15:36:31 2020

@author: marcnol

Figure 1

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


# font = {'family' : 'normal',
#         'weight' : 'normal',
#         'size'   : 18}

# plt.rc('font', **font)

plottingFileExtension = ".pdf"


class analysisHiMmatrix:
    def __init__(self, runParameters, rootFolder="."):
        self.dataFolder = rootFolder + os.sep + "scHiMmatrices"
        self.runParameters = runParameters
        self.rootFolder = rootFolder

    def loadData(self):

        # loads datasets
        fileNameListDataJSON = self.rootFolder + os.sep + self.runParameters["parametersFileName"]

        with open(fileNameListDataJSON) as json_file:
            ListData = json.load(json_file)

        datasetName = list(ListData.keys())[0]
        print("Dataset: {}".format(datasetName))

        outputFileName = self.dataFolder + os.sep + datasetName + "_Cells:" + self.runParameters["action"]

        fileNameParametersJSON = outputFileName + "_parameters.json"
        with open(fileNameParametersJSON) as json_file:
            p = json.load(json_file)

        dataFiles = {}
        dataFiles["ensembleContactProbability"] = "_ensembleContactProbability.npy"
        # dataFiles['uniqueBarcodes']="_uniqueBarcodes.npy"
        dataFiles["SCmatrixCollated"] = "_SCmatrixCollated.npy"
        dataFiles["SClabeledCollated"] = "_SClabeledCollated.npy"

        dataFiles["anchor:4"] = "_anchor:4_ensemble3wayContacts.npy"
        dataFiles["anchor:6"] = "_anchor:6_ensemble3wayContacts.npy"
        dataFiles["anchor:9"] = "_anchor:9_ensemble3wayContacts.npy"
        dataFiles["anchor:10"] = "_anchor:10_ensemble3wayContacts.npy"
        dataFiles["anchor:13"] = "_anchor:13_ensemble3wayContacts.npy"
        dataFiles["anchor:16"] = "_anchor:16_ensemble3wayContacts.npy"

        data = {}

        for idataFile in dataFiles.keys():
            print("Loaded: {}".format(idataFile))
            data[idataFile] = np.load(outputFileName + dataFiles[idataFile]).squeeze()

        runName = loadList(outputFileName + "_runName.csv")
        data["uniqueBarcodes"] = loadList(outputFileName + "_uniqueBarcodes.csv")

        print("Loaded: {}".format("runName"))
        data["runName"] = runName

        self.data = data
        self.dataFiles = dataFiles

    # functions

    def plot2DMatrixSimple(
        self,
        ifigure,
        matrix,
        uniqueBarcodes,
        yticks,
        xticks,
        cmtitle="probability",
        cMin=0,
        cMax=1,
        cm="coolwarm",
        fontsize=12,
        colorbar=False,
        axisTicks=False,
    ):

        pos = ifigure.imshow(matrix, cmap=cm)  # colormaps RdBu seismic
        # plots figure
        if xticks:
            ifigure.set_xlabel("barcode #", fontsize=fontsize)
            if not axisTicks:
                ifigure.set_xticklabels(())
            else:
                print("barcodes:{}".format(uniqueBarcodes))
                # ifigure.set_xticks(np.arange(matrix.shape[0]),uniqueBarcodes)
                ifigure.set_xticklabels(uniqueBarcodes)

        else:
            ifigure.set_xticklabels(())
        if yticks:
            ifigure.set_ylabel("barcode #", fontsize=fontsize)
            if not axisTicks:
                ifigure.set_yticklabels(())
            else:
                # ifigure.set_yticks(np.arange(matrix.shape[0]), uniqueBarcodes)
                ifigure.set_yticklabels(uniqueBarcodes)
        else:
            ifigure.set_yticklabels(())

        for xtick, ytick in zip(ifigure.xaxis.get_majorticklabels(), ifigure.yaxis.get_majorticklabels()):
            xtick.set_fontsize(fontsize)
            ytick.set_fontsize(fontsize)

        if colorbar:
            cbar = plt.colorbar(pos, fraction=0.046, pad=0.04)
            cbar.minorticks_on()
            cbar.set_label(cmtitle)
            pos.set_clim(vmin=cMin, vmax=cMax)
        return pos

    def update_clims(self, cMin, cMax, axes):
        for ax in axes:
            ax.set_clim(vmin=cMin, vmax=cMax)


def loadList(fileName):
    with open(fileName, newline="") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ", quotechar="|")
        runName = []
        for row in spamreader:
            print(", ".join(row))
            if len(runName) > 0:
                runName.append(row)
            else:
                runName = row

    return runName
