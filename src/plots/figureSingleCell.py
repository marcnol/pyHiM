#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 10:41:05 2020

@author: marcnol

produces movies and structures from single cell PWD matrices

"""

#%% imports and plotting settings
import os
import numpy as np
import argparse

import cv2

# from mayavi.mlab import *
# import matplotlib as plt
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from sklearn.neighbors import KernelDensity
from sklearn import manifold
from sklearn.model_selection import GridSearchCV, LeaveOneOut

from matrixOperations.HIMmatrixOperations import getRgFromPWD, getDetectionEffBarcodes, getBarcodesPerCell, kdeFit
from matrixOperations.HIMmatrixOperations import analysisHiMmatrix, getsCoordinatesFromPWDmatrix, sortsCellsbyNumberPWD
from matrixOperations.HIMmatrixOperations import plotDistanceHistograms, write_XYZ_2_pdb


font = {"family": "DejaVu Sans", "weight": "normal", "size": 22}

matplotlib.rc("font", **font)

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
    parser.add_argument("--nRows", help="The number of cells is set by nRows**2")
    parser.add_argument("--pixelSize", help="Pixel Size in um")
    parser.add_argument("--maxDistance", help="Maximum distance for histograms, in um")

    parser.add_argument("--plottingFileExtension", help="By default: svg. Other options: pdf, png")
    parser.add_argument(
        "--shuffle",
        help="Provide shuffle vector: 0,1,2,3... of the same size or smaller than the original matrix. No spaces! comma-separated!",
    )
    parser.add_argument(
        "--ensembleMatrix",
        help="Use if you want ensembleMatrix to be plotted alongside sc matrices",
        action="store_true",
    )
    parser.add_argument("--video", help="Use if you want to output video", action="store_true")
    parser.add_argument(
        "--videoAllcells",
        help="Use if you want all nRows**2 single cells to take part of the video",
        action="store_true",
    )
    parser.add_argument(
        "--plotHistogramMatrix",
        help="Use if you want to plot the PWD histograms for all bin combinations. This is slow!",
        action="store_true",
    )
    parser.add_argument("--minNumberPWD", help="Minimum number of PWD to calculate Rg")
    parser.add_argument("--threshold", help="Maximum accepted PWD to calculate Rg, in px")

    args = parser.parse_args()

    runParameters = {}

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        # rootFolder = "."
        # rootFolder='/home/marcnol/data'+os.sep+'Experiment_18'
        # rootFolder = "/mnt/grey/DATA/docPaper_fullDatasets/updatedDatasets/wt_docTAD_nc14"
        rootFolder = "/home/marcnol/data/updatedDatasets/wt_docTAD_nc14"

    if args.outputFolder:
        outputFolder = args.outputFolder
    else:
        outputFolder = rootFolder + os.sep + "figureSingleCell"

    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
        print("Folder created: {}".format(outputFolder))

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

    if args.pixelSize:
        runParameters["pixelSize"] = float(args.pixelSize)
    else:
        runParameters["pixelSize"] = 0.1

    if args.maxDistance:
        runParameters["maxDistance"] = float(args.maxDistance)
    else:
        runParameters["maxDistance"] = 4.0

    if args.threshold:
        runParameters["threshold"] = float(args.threshold)
    else:
        runParameters["threshold"] = 8

    if args.minNumberPWD:
        runParameters["minNumberPWD"] = args.minNumberPWD
    else:
        runParameters["minNumberPWD"] = 6

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

    if args.plottingFileExtension:
        runParameters["plottingFileExtension"] = "." + args.plottingFileExtension
    else:
        runParameters["plottingFileExtension"] = ".png"

    if args.shuffle:
        runParameters["shuffle"] = args.shuffle
    else:
        runParameters["shuffle"] = 0

    if args.nRows:
        runParameters["nRows"] = int(args.nRows)
    else:
        runParameters["nRows"] = int(10)

    if args.ensembleMatrix:
        runParameters["ensembleMatrix"] = args.ensembleMatrix
    else:
        runParameters["ensembleMatrix"] = False

    if args.video:
        runParameters["video"] = args.video
    else:
        runParameters["video"] = False

    if args.videoAllcells:
        runParameters["videoAllcells"] = args.videoAllcells
    else:
        runParameters["videoAllcells"] = False

    if args.plotHistogramMatrix:
        runParameters["plotHistogramMatrix"] = args.plotHistogramMatrix
    else:
        runParameters["plotHistogramMatrix"] = False

    return rootFolder, outputFolder, runParameters


def returnCellsHighestNumberPWD(sortedValues, n):
    cellID = list()
    Npwd = list()

    for i in range(1, n + 1):
        cellID.append(sortedValues[-i][0])
        Npwd.append(sortedValues[-i][1])

    return cellID, Npwd


def visualize3D(coordinates, colors=[], cmap="hsv", title=[], output="visualize3D.png"):

    fig = plt.figure()
    fig.set_size_inches((10, 10))

    ax = plt.axes(projection="3d")

    xdata, ydata, zdata, barcodeID = [], [], [], []
    for i, r in enumerate(coordinates):
        xdata.append(r[0])
        ydata.append(r[1])
        zdata.append(r[2])
        barcodeID.append(i)

    if len(colors) == 0:
        colors = np.arange(coordinates.shape[0])
        colors = colors.flatten()

    ax.plot3D(xdata, ydata, zdata, "black")
    pos = ax.scatter3D(
        xdata,
        ydata,
        zdata,
        c=colors,
        s=200,
        cmap=cmap,
        marker="o",
        edgecolors="k",
        linewidths=2,
        alpha=0.7,
        vmin=0,
        vmax=1,
    )  #'terrain_r'
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(title)
    ax.set_xticklabels(())
    ax.set_yticklabels(())
    ax.set_zticklabels(())
    fig.savefig(output)
    plt.close(fig)


# def visualize3D_mayavi(coordinates,colors=[],cmap='hsv',title=[],output='visualize3D.png'):

#     fig = plt.figure()
#     fig.set_size_inches((10, 10))

#     ax = plt.axes(projection='3d')

#     ax = plt.axes(projection='3d')

#     xdata, ydata, zdata, barcodeID = [],[],[], []
#     for i,r in enumerate(coordinates):
#         xdata.append(r[0])
#         ydata.append(r[1])
#         zdata.append(r[2])
#         barcodeID.append(i)

#     if len(colors)==0:
#         colors=np.arange(coordinates.shape[0])
#         colors=colors.flatten()

#     # points3d()
#     points3d(xdata, ydata, zdata, np.ones(coordinates.shape[0]), colormap="hsv", scale_factor=.1)
#     plot3d(xdata, ydata, zdata,tube_radius=0.01, colormap='Spectral')

#     ax.plot3D(xdata, ydata, zdata, 'black')
#     pos=ax.scatter3D(xdata, ydata, zdata, c=colors,
#                       s=200,
#                       cmap=cmap,
#                       marker='o',
#                       edgecolors = 'k',
#                       linewidths=2,
#                       alpha=0.7) #'terrain_r'
#     ax.set_xlabel('x')
#     ax.set_ylabel('y')
#     ax.set_zlabel('z')
#     ax.set_title(title)
#     fig.savefig(output)
#     plt.close(fig)


def visualize2D(coordinateList, colors=[], cmap="hsv", titles=[], output="visualize2D.png"):

    nRows = len(coordinateList)

    fig, allAxes = plt.subplots(1, nRows)
    fig.set_size_inches((10 * nRows, 10))
    axes = allAxes.ravel()

    for coordinates, ax, title in zip(coordinateList, axes, titles):
        if len(colors) == 0:
            colors = np.arange(coordinates.shape[0])
            colors = colors.flatten()

        xdata, ydata, zdata, barcodeID = [], [], [], []
        for i, r in enumerate(coordinates):
            xdata.append(r[0])
            ydata.append(r[1])
            zdata.append(r[2])
            barcodeID.append(i)

        ax.plot(xdata, ydata, "black")
        pos = ax.scatter(
            xdata,
            ydata,
            s=200,
            c=colors,
            cmap=cmap,
            marker="o",
            edgecolors="k",
            linewidths=2,
            alpha=0.7,
            vmax=1,
            vmin=0,
        )  #'terrain_r'
        ax.set_xlabel("x")
        ax.set_ylabel("y")

        ax.set_title(title)

        fig.savefig(output)
        plt.close(fig)


def plotTrajectories(HiMdata, runParameters, outputFileNameRoot, cellID, mode="matplotlib"):

    PWDmatrix = SCmatrix[:, :, cellID]
    EnsembleMatrix = 1 / HiMdata.data["ensembleContactProbability"]

    ATACseqMatrix = np.array(HiMdata.ListData[HiMdata.datasetName]["BarcodeColormap"]) / 10
    colors = np.atleast_2d(ATACseqMatrix).flatten()

    singleCellTitle = "Cell #" + str(cellID)
    ensembleTitle = "Ensemble"

    # removes nans
    for i in range(PWDmatrix.shape[0]):
        PWDmatrix[i, i] = 0
        EnsembleMatrix[i, i] = 0

    # gets coordinates and saves in PDB format
    EnsembleMatrix[np.isnan(EnsembleMatrix)] = 0  # removes NaNs from matrix
    coordinatesEnsemble = getsCoordinatesFromPWDmatrix(EnsembleMatrix)

    PWDmatrix[np.isnan(PWDmatrix)] = 0  # removes NaNs from matrix
    coordinates = runParameters["pixelSize"] * getsCoordinatesFromPWDmatrix(PWDmatrix)
    outputFileNamePDB = outputFileNameRoot + "_SingleCellTrajectory:" + str(cellID) + ".pdb"
    write_XYZ_2_pdb(outputFileNamePDB, coordinates)

    # makes plots
    cmap = "tab10"
    output = outputFileNameRoot + "_2DsingleCell:" + str(cellID) + runParameters["plottingFileExtension"]
    visualize2D(
        [coordinates, coordinatesEnsemble],
        colors=colors,
        cmap=cmap,
        titles=[singleCellTitle, ensembleTitle],
        output=output,
    )

    output = outputFileNameRoot + "_3DensembleMatrix" + runParameters["plottingFileExtension"]
    visualize3D(coordinatesEnsemble, colors=colors, cmap=cmap, title=ensembleTitle, output=output)

    output = outputFileNameRoot + "_3DsingleCell:" + str(cellID) + runParameters["plottingFileExtension"]

    if mode == "matplotlib":
        visualize3D(coordinates, colors=colors, cmap=cmap, title=singleCellTitle, output=output)

    # else:
    #     visualize3D_mayavi(coordinates,colors=colors,cmap='rainbow',title=singleCellTitle, output=output)


def plotSCmatrix(HiMdata, cellID, outputFileNameRoot="SCmatrix.png", ensembleMatrix=False, searchPattern="_scMatrix:"):
    datasetName = list(HiMdata.ListData.keys())[0]

    vmax = HiMdata.ListData[datasetName]["iPWD_clim"]
    cmap = HiMdata.ListData[datasetName]["iPWD_cm"]
    SCmatrix = HiMdata.SCmatrixSelected
    singleCellTitle = "Cell #" + str(cellID)

    if ensembleMatrix:
        ensembleTitle = "Ensemble"
        EnsembleMatrix = 1 / HiMdata.data["ensembleContactProbability"]

    PWDmatrix = SCmatrix[:, :, cellID]

    if ensembleMatrix:
        fig, allAxes = plt.subplots(1, 2)
        ax = allAxes.ravel()
        fig.set_size_inches((20, 10))
    else:
        fig, allAxes = plt.subplots(1, 1)
        fig.set_size_inches((10, 10))
        ax = [allAxes]

    uniqueBarcodes = 1 + np.arange(PWDmatrix.shape[0])

    p1 = ax[0].imshow(1 / PWDmatrix, cmap=cmap, vmin=0, vmax=vmax)
    fig.colorbar(p1, ax=ax[0], fraction=0.046, pad=0.04)
    ax[0].set_title(singleCellTitle)
    plt.xticks(np.arange(PWDmatrix.shape[0]), uniqueBarcodes)
    plt.yticks(np.arange(PWDmatrix.shape[0]), uniqueBarcodes)

    if ensembleMatrix:
        p2 = ax[1].imshow(1 / EnsembleMatrix[:, :], cmap=cmap, vmin=0, vmax=vmax)
        fig.colorbar(p2, ax=ax[1], fraction=0.046, pad=0.04)
        ax[1].set_title(ensembleTitle)

    output = outputFileNameRoot + searchPattern + str(cellID) + runParameters["plottingFileExtension"]
    plt.savefig(output)
    plt.close(fig)


def plotsSubplotSCmatrices(HiMdata, nRows, output="subplotMatrices.png"):

    datasetName = list(HiMdata.ListData.keys())[0]

    SCmatrix, sortedValues, nCells = sortsCellsbyNumberPWD(HiMdata)

    # displays plots
    Ncells2Process = nRows ** 2
    cellID, Npwd = returnCellsHighestNumberPWD(sortedValues, Ncells2Process)

    fig, allAxes = plt.subplots(nRows, nRows)
    fig.set_size_inches((50, 50))
    ax = allAxes.ravel()

    cmap = HiMdata.ListData[datasetName]["ContactProbability_cm"]
    vmax = HiMdata.ListData[datasetName]["iPWD_clim"]

    iplot = 0
    for iCell in cellID:
        pos = ax[iplot].imshow(1 / SCmatrix[:, :, iCell], cmap=cmap, vmin=0, vmax=vmax)
        ax[iplot].set_xticklabels(())
        ax[iplot].set_yticklabels(())
        ax[iplot].set_axis_off()
        ax[iplot].set_title(str(iCell))

        iplot += 1

    plt.savefig(output)
    plt.close(fig)
    return cellID, SCmatrix


def makesVideo(folder, video_name, searchPattern):

    images = [img for img in os.listdir(os.path.dirname(folder)) if img.endswith(".png") and searchPattern in img]
    if len(images) > 0:
        frame = cv2.imread(os.path.join(os.path.dirname(folder), images[0]))
        height, width, layers = frame.shape

        video = cv2.VideoWriter(video_name, 0, 5, (width, height))

        for image in images:
            video.write(cv2.imread(os.path.join(os.path.dirname(folder), image)))

        cv2.destroyAllWindows()
        video.release()
    else:
        print("Sorry, no images found fitting the pattern {} in this folder: {}".format(searchPattern, folder))


def plotsBarcodesPerCell(SCmatrix, runParameters, outputFileNameRoot="./"):

    numBarcodes = getBarcodesPerCell(SCmatrix)
    maxNumberBarcodes = SCmatrix.shape[0]

    fig, ax = plt.subplots()
    fig.set_size_inches((10, 10))
    ax.hist(numBarcodes, bins=range(2, maxNumberBarcodes + 1))
    ax.set_xlabel("number of barcodes")
    ax.set_ylabel("counts")

    output = outputFileNameRoot + "_SChistBarcodesPerCell" + runParameters["plottingFileExtension"]
    plt.savefig(output)
    plt.close(fig)


def plotsBarcodesEfficiencies(SCmatrix, runParameters, uniqueBarcodes, outputFileNameRoot="./"):

    eff = getDetectionEffBarcodes(SCmatrix)

    fig, ax = plt.subplots()
    fig.set_size_inches((10, 10))
    ax.bar(uniqueBarcodes, eff)
    ax.set_xlabel("barcode ID")
    ax.set_ylabel("efficiency")
    ax.set_xticks(np.arange(len(eff)))
    ax.set_xticklabels(uniqueBarcodes)

    output = outputFileNameRoot + "_SCBarcodesEfficiency" + runParameters["plottingFileExtension"]
    plt.savefig(output)
    plt.close(fig)


def plotsRgvalues(
    HiMdata,
    nRows,
    runParameters,
    outputFileName="./RgValues.png",
    minNumberPWD=6,
    threshold=6,
    bandwidths=10 ** np.linspace(-1.5, 0, 20),
):

    print("Threshold = {} px | min number PWDs = {}".format(threshold, minNumberPWD))

    SCmatrix, sortedValues, nCells = sortsCellsbyNumberPWD(HiMdata)
    Ncells2Process = nRows ** 2
    selectedCellsIDs, Npwd = returnCellsHighestNumberPWD(sortedValues, Ncells2Process)

    # calculates Rg for all cells
    RgList = list()
    for cellID in selectedCellsIDs:
        RgList.append(
            runParameters["pixelSize"]
            * getRgFromPWD(SCmatrix[:, :, cellID], minNumberPWD=minNumberPWD, threshold=threshold)
        )

    RgListArray = np.array(RgList)
    maxRange = 1.5 * RgListArray.max()

    # plots single Rg as bars
    fig, ax = plt.subplots(
        1, 1, figsize=(10, 10), sharex=True, sharey=True, subplot_kw={"xlim": (0, maxRange), "ylim": (-0.02, 1.2)}
    )
    ax.plot(RgListArray, np.full_like(RgListArray, -0.01), "|k", markeredgewidth=0.5)

    # calculates and displays median
    mean = np.nanmedian(RgListArray)
    print("Median Rg = {}".format(mean))
    ax.set_ylabel("counts")
    ax.set_xlabel("Rg, " + r"$\mu$m")
    x_d = np.linspace(0, maxRange, 100)

    # KDE fit
    # fist finds best bandwidth
    print("Calculating optimal KDE bandwidth...")

    grid = GridSearchCV(KernelDensity(kernel="gaussian"), {"bandwidth": bandwidths}, cv=LeaveOneOut())
    grid.fit(RgListArray[:, None])
    bandwidth = grid.best_params_["bandwidth"]
    print("bandwidth = {}".format(bandwidth))

    # calculates KDE with optimal bandwidth
    logprob, kde = kdeFit(RgListArray, x_d, bandwidth=bandwidth)
    kde_params = kde.get_params()
    maxlogprob = logprob.max()
    ax.fill_between(x_d, np.exp(logprob) / np.exp(maxlogprob), alpha=0.3)

    mean = x_d[np.argmax(logprob, axis=0)]
    print("KDE max Rg = {}".format(mean))

    ax.axvline(x=mean, color="black", linestyle=(0, (5, 5)))

    plt.savefig(outputFileName)
    plt.close(fig)
    return RgList


def makesPlotHistograms(
    HiMdata,
    runParameters,
    outputFileName="./HiMhistograms.png",
    mode="KDE",
    kernelWidth=0.25,
    optimizeKernelWidth=False,
):
    SCmatrix, sortedValues, nCells = sortsCellsbyNumberPWD(HiMdata)

    plotDistanceHistograms(
        SCmatrix,
        runParameters["pixelSize"],
        outputFileName,
        mode="KDE",
        kernelWidth=0.25,
        optimizeKernelWidth=False,
        maxDistance=runParameters["maxDistance"],
    )


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

    nDatasets = len(HiMdata.data["runName"])

    if outputFolder == "none":
        outputFolder = HiMdata.dataFolder

    outputFileNameRoot = (
        outputFolder
        + os.sep
        + "Fig_SCmatrices"
        + "_dataset1:"
        + HiMdata.datasetName
        + "_label:"
        + runParameters["label"]
        + "_action:"
        + runParameters["action"]
    )
    datasetName = list(HiMdata.ListData.keys())[0]
    print("Data output: {}".format(outputFileNameRoot))

    # "makes subplots with sc 1/PWD matrices"
    print("\n>>>Plotting subplots with 1/PWD matrices<<<\n")
    nRows = runParameters["nRows"]
    output = outputFileNameRoot + "_scMatrices" + runParameters["plottingFileExtension"]
    cellID_most_PWDs, SCmatrix = plotsSubplotSCmatrices(HiMdata, nRows, output=output)

    # "calculates the number of barcodes per cell and makes histograms"
    print("\n>>>Calculating distribution of barcodes<<<\n")
    plotsBarcodesPerCell(SCmatrix, runParameters, outputFileNameRoot=outputFileNameRoot)

    # "calculates the detection efficiency for each barcode"
    print("\n>>>Calculating detection efficiency distribution<<<\n")
    plotsBarcodesEfficiencies(
        SCmatrix, runParameters, list(HiMdata.data["uniqueBarcodes"]), outputFileNameRoot=outputFileNameRoot
    )

    # "calculates the Rg for each cell from the PWD sc matrix"
    print("\n>>>Calculating Rg distributions<<<\n")
    output = outputFileNameRoot + "_RgValues" + runParameters["plottingFileExtension"]
    RgList = plotsRgvalues(
        HiMdata,
        nRows,
        runParameters,
        outputFileName=output,
        minNumberPWD=int(runParameters["minNumberPWD"]),
        threshold=float(runParameters["threshold"]),
        bandwidths=10 ** np.linspace(-1, 0, 20),
    )

    # plots distance histograms
    if runParameters["plotHistogramMatrix"]:
        print("\n>>>Plotting distance histograms<<<\n")
        output = outputFileNameRoot + "_HistogramPWDs" + runParameters["plottingFileExtension"]
        makesPlotHistograms(
            HiMdata, runParameters, outputFileName=output, mode="KDE", kernelWidth=0.25, optimizeKernelWidth=False
        )

    # "plots trajectories for selected cells"
    print("\n>>>Plotting trajectories for selected cells<<<\n")
    if "CellIDs" in HiMdata.ListData[datasetName].keys():
        CellIDs = HiMdata.ListData[datasetName]["CellIDs"]
        print("CellIDs to process: {}".format(CellIDs))
        for cellID in CellIDs:
            if cellID < HiMdata.SCmatrixSelected.shape[2]:
                #  Plots sc 1/PWD matrix and ensemble 1/PWD matrix together
                plotSCmatrix(HiMdata, cellID, outputFileNameRoot, ensembleMatrix=runParameters["ensembleMatrix"])

                # plots trajectories
                plotTrajectories(HiMdata, runParameters, outputFileNameRoot, cellID)

    # "makes video of SC matrix for selected cells"
    if runParameters["video"]:
        print("\n>>>Making video<<<\n")
        if runParameters["videoAllcells"]:
            searchPattern = "_scMostPWD:"
            for cellID in cellID_most_PWDs:
                plotSCmatrix(
                    HiMdata,
                    cellID,
                    outputFileNameRoot,
                    ensembleMatrix=runParameters["ensembleMatrix"],
                    searchPattern=searchPattern,
                )
        else:
            searchPattern = "_scMatrix:"

        video_name = os.path.dirname(outputFileNameRoot) + os.sep + "video.avi"

        makesVideo(outputFileNameRoot, video_name, searchPattern)

    print("\nDone\n\n")
