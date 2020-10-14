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
import numpy.linalg as npl
from sklearn import manifold

import networkx as nx
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
    
# import matplotlib as plt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json, csv
from matrixOperations.HIMmatrixOperations import plotDistanceHistograms, plotMatrix

from matrixOperations.HIMmatrixOperations import analysisHiMmatrix, normalizeMatrix, shuffleMatrix, plotScalogram

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

    args = parser.parse_args()

    runParameters = {}
    runParameters["pixelSize"] = 0.1

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = "."
        # rootFolder='/home/marcnol/data'+os.sep+'Experiment_18'

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
        runParameters["label"] = "doc"

    if args.action:
        runParameters["action"] = args.action
    else:
        runParameters["action"] = "labeled"

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

    return rootFolder, outputFolder, runParameters


def returnCellsHighestNumberPWD(sortedValues, n):
    cellID = list()
    Npwd = list()
    
    for i in range(1,n+1):
        cellID.append(sortedValues[-i][0])
        Npwd.append(sortedValues[-i][1])
    
    return cellID, Npwd

def visualize3D(coordinates,colors=[],cmap='hsv',title=[]):
   
    fig = plt.figure()
    fig.set_size_inches((10, 10))
    
    ax = plt.axes(projection='3d')
        
    ax = plt.axes(projection='3d')
 
    xdata, ydata, zdata, barcodeID = [],[],[], [] 
    for i,r in enumerate(coordinates):
        xdata.append(r[0])    
        ydata.append(r[1])    
        zdata.append(r[2])    
        barcodeID.append(i)

    if len(colors)==0:
        colors=np.arange(coordinates.shape[0])
        colors=colors.flatten()
            
    ax.plot3D(xdata, ydata, zdata, 'black')
    pos=ax.scatter3D(xdata, ydata, zdata, c=colors,
                     s=200,
                     cmap=cmap,
                     marker='o',
                     edgecolors = 'k',
                     linewidths=2,
                     alpha=0.7) #'terrain_r'
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(title)

def visualize2D(coordinateList,colors=[],cmap='hsv',titles=[]):
    
    nRows=len(coordinateList)

    
    fig, allAxes = plt.subplots(1,nRows)
    fig.set_size_inches((10*nRows, 10))
    axes=allAxes.ravel()

    for coordinates,ax,title in zip(coordinateList,axes,titles):
        if len(colors)==0:
            colors=np.arange(coordinates.shape[0])
            colors=colors.flatten()


            
        xdata, ydata, zdata, barcodeID = [],[],[], [] 
        for i,r in enumerate(coordinates):
            xdata.append(r[0])    
            ydata.append(r[1])    
            zdata.append(r[2])    
            barcodeID.append(i)
            
        ax.plot(xdata, ydata, 'black')
        pos=ax.scatter(xdata, ydata, s=200, c=colors,
                         cmap=cmap,
                         marker='o',
                         edgecolors = 'k',
                         linewidths=2,
                         alpha=0.7) #'terrain_r'
        ax.set_xlabel('')
        ax.set_ylabel('')

        ax.set_title(title)
            
def getsCoordinatesFromPWDmatrix(matrix):
    ## multi-dimensional scaling to get coordinates from PWDs
    # make sure meanSCmatrix is symmetric
    matrix = 0.5 * (matrix + np.transpose(matrix))
    # run metric mds
    verbosity = 0  # default: 0, quite verbose: 2
    mds = manifold.MDS(
        n_components=3,
        metric=True,
        n_init=20,
        max_iter=3000,
        verbose=verbosity,
        eps=1e-9,
        n_jobs=1,
        random_state=1,
        dissimilarity="precomputed",  # euclidean | precomputed
    )
    
    XYZ = mds.fit(matrix).embedding_


    return XYZ
 

#%%
# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    print(">>> Producing HiM matrix")
    rootFolder, outputFolder, runParameters = parseArguments()
    
    rootFolder = "/mnt/grey/DATA/docPaper_fullDatasets/updatedDatasets/wt_docTAD_nc14"
    outputFolder = "/home/marcnol/Downloads/test"
    HiMdata = analysisHiMmatrix(runParameters, rootFolder)

    HiMdata.loadData()

    # panel C: contact probability matrix

    # matrix = HiMdata.data["ensembleContactProbability"]
    # matrix=normalizeMatrix(matrix)

    # cScale = matrix.max() / runParameters["scalingParameter"]
    # print("scalingParameters, scale={}, {}".format(runParameters["scalingParameter"], cScale))

    nCells = HiMdata.nCellsLoaded()

    nDatasets = len(HiMdata.data["runName"])

    if outputFolder == "none":
        outputFolder = HiMdata.dataFolder

    outputFileName = (
        outputFolder
        + os.sep
        + "Fig_SCmatrices"
        + "_dataset1:"
        + HiMdata.datasetName
        + "_label:"
        + runParameters["label"]
        + "_action:"
        + runParameters["action"]
        + runParameters["plottingFileExtension"]
    )

    SCmatrix = HiMdata.data["SCmatrixCollated"]
    nCells=SCmatrix.shape[2]
    print("Number of cells loaded: {}".format(nCells))
    
  

    
    #%%
    # finds the number of barcodes detected per cell.
    cmap = 'coolwarm'
    nBarcodePerCell = list()
    values = list()
    dtype = [('cellID', int), ('nPWD', int)]
    
    for iCell in range(nCells):
        SCmatrixCell = SCmatrix[:,:,iCell]
        nPWD = int(np.count_nonzero(~np.isnan(SCmatrixCell))/2)
        nBarcodePerCell.append(nPWD)
        values.append((iCell,nPWD))        
        
    # plt.hist(nBarcodePerCell)

    valuesArray = np.array(values, dtype=dtype)       # create a structured array
    sortedValues = np.sort(valuesArray, order='nPWD')                        
    
    nRows=10
    Ncells2Process = nRows**2
    
    cellID,Npwd = returnCellsHighestNumberPWD(sortedValues, Ncells2Process)
    
    fig, allAxes = plt.subplots(nRows,nRows)
    fig.set_size_inches((50, 50))
    ax=allAxes.ravel()

    iplot=0    
    for iCell in cellID:
        pos = ax[iplot].imshow(1/SCmatrix[:,:,iCell],cmap=cmap,vmin=0,vmax=1.2)
        ax[iplot].set_xticklabels(())
        ax[iplot].set_yticklabels(())
        ax[iplot].set_axis_off()
        ax[iplot].set_title(str(iCell))

        iplot+=1

    plt.savefig(outputFileName)
    
    #%%
    #  Plots SC matrix and ensemble matrix together

    vmax=1.2
    cellID=14433
    singleCellTitle="Cell #"+str(cellID)
    ensembleTitle = "Ensemble"
    
    EnsembleMatrix= 1/HiMdata.data["ensembleContactProbability"]

    PWDmatrix=SCmatrix[:,:,cellID]
    
    fig, allAxes = plt.subplots(1,2)
    fig.set_size_inches((10,10))
    
    ax=allAxes.ravel()
    p1=ax[0].imshow(1/PWDmatrix, cmap='Reds',vmin=0,vmax=vmax)
    fig.colorbar(p1,ax=ax[0],fraction=0.046, pad=0.04)
    ax[0].set_title(singleCellTitle)
    
    p2=ax[1].imshow(1/EnsembleMatrix[:,:],cmap='Reds',vmin=0,vmax=vmax)
    fig.colorbar(p2,ax=ax[1],fraction=0.046, pad=0.04)
    ax[1].set_title(ensembleTitle)
 
    # removes nans
    for i in range(PWDmatrix.shape[0]):
        PWDmatrix[i,i]=0
        EnsembleMatrix[i,i]=0
        
    # plots trajectories
        
    coordinatesEnsemble = getsCoordinatesFromPWDmatrix(EnsembleMatrix)
    coordinates = 0.1*getsCoordinatesFromPWDmatrix(PWDmatrix)

    ATACseqMatrix = np.array(HiMdata.ListData[HiMdata.datasetName]["BarcodeColormap"]) / 10

    colors = np.atleast_2d(ATACseqMatrix).flatten()
    
    visualize2D([coordinates,coordinatesEnsemble],colors=colors,cmap='rainbow',titles=[singleCellTitle,ensembleTitle])

    visualize3D(coordinatesEnsemble,colors=colors,cmap='rainbow',title=ensembleTitle)

    visualize3D(coordinates,colors=colors,cmap='rainbow',title=singleCellTitle)

        


    print("\nDone\n\n")
