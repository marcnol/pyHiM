#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Nov 10 2020

@author: markus

plot a mixed matrix for the "brightest spots" subset

"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

# from matrixOperations.alignBarcodesMasks import plotMatrix


#%% define and loads datasets
def parseArguments():
    # [parsing arguments]
    parser = argparse.ArgumentParser()
    parser.add_argument("-F1", "--rootFolder1", help="Folder with dataset 1")
    parser.add_argument("-F2", "--rootFolder2", help="Folder with dataset 2")
    parser.add_argument("-O", "--outputFolder", help="Folder for outputs")

    args = parser.parse_args()

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
        outputFolder = "none"

    return rootFolder1, rootFolder2, outputFolder


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    # %% global settings
    p = {}
    p["pixelSize"] = 0.1

    rootFolder1, rootFolder2, outputFolder = parseArguments()  # F1 all, F2 brightest

    # %% define paths
    barcodeFile = "scHiMmatrices/wt_docTAD_nc14_label:docHigh_action:labeled_uniqueBarcodes.csv"

    dataFile1 = "scHiMmatrices/wt_docTAD_nc14_label:doc_action:labeled_ensembleContactProbability.npy"
    dataFile2 = "scHiMmatrices/wt_docTAD_nc14_label:docHigh_action:labeled_ensembleContactProbability.npy"

    SClabeled1 = "scHiMmatrices/wt_docTAD_nc14_label:doc_action:labeled_SClabeledCollated.npy"
    SClabeled2 = "scHiMmatrices/wt_docTAD_nc14_label:docHigh_action:labeled_SClabeledCollated.npy"

    # %% load the data

    # load ensemble contact probability maps
    fn1 = os.path.join(rootFolder1, dataFile1)
    fn2 = os.path.join(rootFolder2, dataFile2)

    contactMap1 = np.load(fn1)
    contactMap2 = np.load(fn2)

    # load array with info with cells were "labeled", i.e. selected for plotting
    fn1_SClabeled = os.path.join(rootFolder1, SClabeled1)
    fn2_SClabeled = os.path.join(rootFolder2, SClabeled2)

    SClabeled1 = np.load(fn1_SClabeled)
    SClabeled2 = np.load(fn2_SClabeled)

    numCells1 = np.sum(SClabeled1 == 1)
    numCells2 = np.sum(SClabeled2 == 1)

    # load barcodes
    fnBarcodes = os.path.join(rootFolder2, barcodeFile)
    commonSetUniqueBarcodes = np.loadtxt(fnBarcodes).astype(int)

    # %% mix it

    contactMap_mix = contactMap1.copy()
    sel = np.tril(np.full(contactMap1.shape, True), -1)  # -1 to exclude the diagonal
    contactMap_mix[sel] = contactMap2[sel]

    contactMap1_lin = contactMap1[sel]
    contactMap2_lin = contactMap2[sel]

    corrcoef_matrix = np.corrcoef(contactMap1_lin, contactMap2_lin)
    print("Pearson correlation coefficient = {}".format(corrcoef_matrix[0, 1]))

    # %% plots results
    mode, clim, cMin = "counts", 0.6, 0.0  # mode = counts -> just use matrix as is, no scaling
    minVal = np.nanmin(contactMap_mix)
    maxVal = np.nanmax(contactMap_mix)
    print("Matrix to be plotted: min {}, max {}".format(minVal, maxVal))

    figtitle = "Contact probability map"
    nCells = "{}\{}".format(numCells2, numCells1)
    numberROIs = 0
    uniqueBarcodes = range(1, contactMap_mix.shape[0] + 1)
    cmtitle = "Contact probability"
    matrixPlot = contactMap_mix

    fig = plt.figure(figsize=(6, 6))
    pos = plt.imshow(matrixPlot, cmap="coolwarm")  # colormaps RdBu seismic
    plt.xlabel("Barcode #")
    plt.ylabel("Barcode #")
    plt.title(
        figtitle + " | " + str(matrixPlot.shape[0]) + " barcodes | n=" + str(nCells) + " | ROIs=" + str(numberROIs)
    )
    plt.xticks(np.arange(matrixPlot.shape[0]), uniqueBarcodes)
    plt.yticks(np.arange(matrixPlot.shape[0]), uniqueBarcodes)
    cbar = plt.colorbar(pos, fraction=0.046, pad=0.04)
    # cbar.minorticks_on()
    cbar.set_label(cmtitle)
    plt.clim(cMin, clim)

    # save to PNG
    saveExt = "png"  # png, svg, jpg
    fnSave = os.path.join(outputFolder, "label:doc_action:labeled_label:docHigh_action:labeled." + saveExt)
    fig.savefig(fnSave)

    # save also the npy
    fnSave2 = os.path.join(outputFolder, "label:doc_action:labeled_label:docHigh_action:labeled." + "npy")
    np.save(fnSave2, matrixPlot)

    # %% plot 4M profiles
    anchors = [17]  # 1-based counting, bin number, not RT

    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]
    lwbase = plt.rcParams["lines.linewidth"]
    thin, thick = lwbase / 2, lwbase * 8

    for anchor in anchors:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))

        profile1 = contactMap1[:, anchor - 1]
        profile2 = contactMap2[:, anchor - 1]

        ax.plot(profile1, ls="-", color="magenta")
        ax.plot(profile2, ls="--", color="darkmagenta")

        ax.axvline(x=anchor - 1, color=colors[4], lw=thick, alpha=0.5)
        ax.set_ylim(0, 1)
        ax.set_xlim(0, profile1.shape[0] - 1)

        ax.set_xticks(np.arange(profile1.shape[0]))
        # ax.set_xticklabels(commonSetUniqueBarcodes)
        ax.set_xticklabels(np.arange(profile1.shape[0]) + 1)

        # ax.set_yticks([0, 0.5, 1])

        ax.set_xlabel("Barcode #")
        ax.set_ylabel("Probability")

        ax.legend(("Dorsal ectoderm", "Subset transcription\nhotspots"))

        saveExt = "png"  # png, svg, jpg
        outputFileName = os.path.join(outputFolder, "4M_doc_subset_anchor_{}.{}".format(anchor, saveExt))
        fig.savefig(outputFileName)
