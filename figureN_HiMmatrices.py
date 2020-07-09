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
from alignBarcodesMasks import plotDistanceHistograms, plotMatrix
# import scaleogram as scg

from HIMmatrixOperations import analysisHiMmatrix,normalizeMatrix,shuffleMatrix,plotScalogram,listsSCtoKeep

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
    parser.add_argument("--shuffle", help="Provide shuffle vector: 0,1,2,3... of the same size or smaller than the original matrix. No spaces! comma-separated!")
    parser.add_argument("--scalogram", help="Use if you want scalogram image to be displayed", action="store_true")
    parser.add_argument("--type", help="Provide one of the following: PWD, contact, iPWD")
    parser.add_argument("--pixelSize", help="Provide pixelSize in um")

    args = parser.parse_args()

    runParameters = {}

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = "."
        # rootFolder='/home/marcnol/data'+os.sep+'Experiment_18'

    if args.outputFolder:
        outputFolder = args.outputFolder
    else:
        outputFolder = 'none'
        
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
        runParameters["plottingFileExtension"] = '.'+args.plottingFileExtension
    else:
        runParameters["plottingFileExtension"] = '.svg'

    if args.shuffle:
        runParameters["shuffle"] = args.shuffle
    else:
        runParameters["shuffle"] = 0

    if args.scalogram:
        runParameters['scalogram']= args.scalogram
    else:
        runParameters['scalogram'] = False

    if args.type:
        runParameters["type"] = args.type
    else:
        runParameters["type"] = 'contact'

    if args.pixelSize:
        runParameters["pixelSize"] = float(args.pixelSize)
    else:
        runParameters["pixelSize"] = 0.1

    return rootFolder, outputFolder,runParameters

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    print(">>> Producing HiM matrix")
    rootFolder, outputFolder, runParameters = parseArguments()

    # loads datasets: parameter files
    fileNameListDataJSON = rootFolder + os.sep + runParameters["parametersFileName"]
    with open(fileNameListDataJSON) as json_file:
        ListData = json.load(json_file)
            
    dataSets= list(ListData.keys())
  
    for idataSet in dataSets:
        
        Samples = ListData[idataSet]['Folders']

        if outputFolder=='none':
            outputFolder = rootFolder

        outputFileName = (
            outputFolder
            + os.sep
            + "Fig_HiMmatrix"
            + "_dataset1:"
            + idataSet
            + "_label:"
            + runParameters["label"]
            + "_action:"
            + runParameters["action"]
            + runParameters["plottingFileExtension"]
        )
        
        Nplots=len(Samples)

    
        fig2 = plt.figure(constrained_layout=False,figsize=(5*Nplots, 5), dpi=80, facecolor='w', edgecolor='k')
        # nCols=np.ceil(len(anchors)/2).astype(int)
        nCols=Nplots
        nRows=1
        spec2 = gridspec.GridSpec(ncols=nCols, nrows=nRows, figure=fig2)
      
        FigList,Yticks, Xticks =[], [], []
        for iRow in range(nRows):
            for iCol in range(nCols):
                FigList.append(fig2.add_subplot(spec2[iRow, iCol]))
                if iRow==nRows-1:
                    Xticks.append(False)
                else:
                    Xticks.append(False)
                if iCol==0:
                    Yticks.append(True)
                else:
                    Yticks.append(False)

        FigLabels = [isample.split(os.sep)[-2] for isample in Samples]
        legendList=[False]*len(Samples)
        colorbar=[False]*len(Samples)
        # colorbar[-1]=True
        
        for isample, ifigure, iFigLabel, yticks, xticks,legend, icolorbar in zip(Samples, FigList, FigLabels, Yticks, Xticks,legendList,colorbar):

            HiMdata = analysisHiMmatrix(runParameters, os.path.dirname(isample))
            HiMdata.loadData()
            
            if runParameters["type"]=='contact':
                matrix=HiMdata.data["ensembleContactProbability"]
                # matrix=normalizeMatrix(matrix)
                cScale = matrix.max() / runParameters["scalingParameter"]
            elif runParameters["type"]=='PWD':
                matrixSC=HiMdata.data["SCmatrixCollated"]
                cells2Plot = listsSCtoKeep(runParameters, HiMdata.data["SClabeledCollated"])
                matrix = runParameters["pixelSize"] * np.nanmedian(matrixSC[:,:,cells2Plot], axis=2)
                cScale = 3*np.nanmedian(matrix) / runParameters["scalingParameter"]
                del matrixSC

            elif runParameters["type"]=='iPWD':
                matrixSC=HiMdata.data["SCmatrixCollated"]
                cells2Plot = listsSCtoKeep(runParameters, HiMdata.data["SClabeledCollated"])
                matrixPWD = runParameters["pixelSize"] * np.nanmedian(matrixSC[:,:,cells2Plot], axis=2)
                matrix = np.reciprocal(matrixPWD)
                cScale = 3*np.reciprocal(np.nanmedian(matrix)) / runParameters["scalingParameter"]
                del matrixPWD, matrixSC

            print("scalingParameters, scale={}, {}".format(runParameters["scalingParameter"],cScale))
            
            nCells = HiMdata.nCellsLoaded()
             
            nDatasets = len(HiMdata.data["runName"])
        
            if runParameters["shuffle"]==0:
                index=range(matrix.shape[0])
            else:
                index=[int(i) for i in runParameters["shuffle"].split(',')]
                matrix=shuffleMatrix(matrix,index)

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
                figTitle = iFigLabel
            )
        
            del HiMdata, matrix
        
        cbar_ax = fig2.add_axes([0.92, 0.20, 0.005, 0.6])
        cbar = fig2.colorbar(f2_ax1_im, cax=cbar_ax, fraction=0.046, pad=0.04)
        ticklabs = cbar.ax.get_yticklabels()
        ticklabs1=["{:04.2f}".format(i*cScale/(len(ticklabs)-1)) for i in range(len(ticklabs))]
        cbar.ax.set_yticklabels(ticklabs1, fontsize=12)
        cbar.set_label(runParameters["type"],fontsize=15)
            
        # HiMdata.update_clims(0, cScale, f1)
        print('Output written to {}'.format(outputFileName))
        plt.savefig(outputFileName)
        titleText="N = {} | n = {}".format(nCells,nDatasets)
        print('Title: {}'.format(titleText))
        print("Output figure: {}".format(outputFileName))
        

    print("\nDone\n\n")
