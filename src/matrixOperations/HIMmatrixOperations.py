#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 21:21:02 2020

@author: marcnol

contains functions and classes needed for the analysis and plotting of HiM matrices

"""

# =============================================================================
# IMPORTS
# =============================================================================


import numpy as np
import os
import glob
import json, csv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import interpolate
from pylab import contourf, colorbar
from numba import jit
from tqdm import trange

from scipy.io import loadmat
from sklearn import manifold
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.neighbors import KernelDensity

from astropy.table import Table, vstack

from fileProcessing.fileManagement import writeString2File

from fileProcessing.fileManagement import isnotebook

# =============================================================================
# CLASSES
# =============================================================================


class analysisHiMmatrix:
    """
    this class is used for loading data processed by processHiMmatrix.py
    Main use is to produce paper quality figures of HiM matrices, 3-way interaction matrices and HiM matrix ratios
    """

    def __init__(self, runParameters, rootFolder="."):
        self.dataFolder = rootFolder + os.sep + "scHiMmatrices"
        self.runParameters = runParameters
        self.rootFolder = rootFolder
        self.data = []
        self.dataFiles = []
        self.folders2Load = []
        self.numberBarcodes = 0

    def loadData(self):
        """
        loads dataset

        Returns
        -------
        self.foldes2Load contains the parameters used for the processing of HiM matrices.
        self.dataFiles dictionary containing the extensions needed to load data files
        self.data dictionary containing the datasets loaded
        """

        # loads datasets: parameter files
        fileNameListDataJSON = self.rootFolder + os.sep + self.runParameters["parametersFileName"]
        with open(fileNameListDataJSON) as json_file:
            ListData = json.load(json_file)

        datasetName = list(ListData.keys())[0]
        print("Dataset: {}".format(datasetName))

        outputFileName = (
            self.dataFolder
            + os.sep
            + datasetName
            + "_label:"
            + self.runParameters["label"]
            + "_action:"
            + self.runParameters["action"]
        )

        fileNameParametersJSON = outputFileName + "_parameters.json"
        with open(fileNameParametersJSON) as json_file:
            folders2Load = json.load(json_file)
        print("Loading parameter file:".format(fileNameParametersJSON))

        # Creates filenames to be loaded
        dataFiles = {}
        dataFiles["ensembleContactProbability"] = "_ensembleContactProbability.npy"
        dataFiles["SCmatrixCollated"] = "_SCmatrixCollated.npy"
        dataFiles["SClabeledCollated"] = "_SClabeledCollated.npy"

        if "3wayContacts_anchors" in ListData[datasetName]:
            for iAnchor in ListData[datasetName]["3wayContacts_anchors"]:
                newKey = "anchor:" + str(iAnchor - 1)
                dataFiles[newKey] = "_" + newKey + "_ensemble3wayContacts.npy"
        else:
            print("No anchors found")

        # loads datasets: numpy matrices
        data = {}
        print("Loading datasets from: {}".format(outputFileName))
        for idataFile in dataFiles.keys():
            print("Loaded: {}: <{}>".format(idataFile, os.path.basename(outputFileName + dataFiles[idataFile])))
            data[idataFile] = np.load(outputFileName + dataFiles[idataFile]).squeeze()

        # loads datasets: lists
        runName = loadList(outputFileName + "_runName.csv")
        data["runName"] = runName
        print("Loaded runNames: {}".format(data["runName"]))

        data["uniqueBarcodes"] = loadList(outputFileName + "_uniqueBarcodes.csv")
        print("Loaded barcodes #: {}".format(data["uniqueBarcodes"]))
        self.numberBarcodes = len(data["uniqueBarcodes"])

        print("Total number of cells loaded: {}".format(data["SCmatrixCollated"].shape[2]))
        print("Number Datasets loaded: {}".format(len(data["runName"])))

        # Exports data
        self.data = data
        self.dataFiles = dataFiles
        self.folders2Load = folders2Load
        self.ListData = ListData
        self.datasetName = datasetName

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
        nCells=0,
        nDatasets=0,
        showTitle=False,
        figTitle="",
    ):

        pos = ifigure.imshow(matrix, cmap=cm)  # colormaps RdBu seismic

        if showTitle:
            titleText = "{} | N = {} | n = {}".format(figTitle, nCells, nDatasets)
            ifigure.title.set_text(titleText)

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
            cbar = plt.colorbar(pos, ax=ifigure, fraction=0.046, pad=0.04)
            cbar.minorticks_on()
            cbar.set_label(cmtitle, fontsize=float(fontsize) * 0.85)
            pos.set_clim(vmin=cMin, vmax=cMax)

        pos.set_clim(vmin=cMin, vmax=cMax)

        return pos

    def update_clims(self, cMin, cMax, axes):
        for ax in axes:
            ax.set_clim(vmin=cMin, vmax=cMax)

    def plot1Dprofile1Dataset(self, ifigure, anchor, iFigLabel, yticks, xticks):

        prop_cycle = plt.rcParams["axes.prop_cycle"]
        colors = prop_cycle.by_key()["color"]
        lwbase = plt.rcParams["lines.linewidth"]
        thin, thick = lwbase / 2, lwbase * 3

        profile = self.data["ensembleContactProbability"][:, anchor - 1]
        x = np.linspace(0, profile.shape[0], num=profile.shape[0], endpoint=True)
        # f = interp1d(x, profile,kind = 'linear') # linear
        tck = interpolate.splrep(x, profile, s=0)
        xnew = np.linspace(0, profile.shape[0], num=100, endpoint=True)
        ynew = interpolate.splev(xnew, tck, der=0)
        if self.runParameters["splines"]:
            ifigure.plot(xnew, ynew, "-")  # x, profile, 'o',
        else:
            ifigure.plot(x, profile, "-")  # x, profile, 'o',

        ifigure.set_xlim([0, profile.shape[0]])
        ifigure.axvline(x=anchor - 0.5, color=colors[4], lw=thick, alpha=0.5)
        ifigure.set_ylim([0, self.runParameters["cAxis"]])

        if xticks:
            ifigure.set_xlabel("barcode #", fontsize=self.runParameters["fontsize"])
            if not self.runParameters["axisTicks"]:
                ifigure.set_xticklabels(())
            else:
                ifigure.set_xticklabels(self.data["uniqueBarcodes"])
        else:
            ifigure.set_xticklabels(())

        if yticks:
            ifigure.set_ylabel("Probability", fontsize=self.runParameters["fontsize"])
            if not self.runParameters["axisTicks"]:
                ifigure.set_yticklabels(())
            else:
                ifigure.set_yticks([0, self.runParameters["cAxis"] / 2, self.runParameters["cAxis"]])
        else:
            ifigure.set_yticklabels(())

    def nCellsLoaded(self):
        if self.runParameters["action"] == "labeled":
            cellswithLabel = [idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x > 0]
            nCells = len(cellswithLabel)
        elif self.runParameters["action"] == "unlabeled":
            cellswithLabel = [idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x == 0]
            nCells = len(cellswithLabel)
        else:
            nCells = self.data["SCmatrixCollated"].shape[2]
        print("nCells selected with label: {}".format(nCells))
        return nCells

    def retrieveSCmatrix(self):
        """
        retrieves single cells that have the label requested

        Returns
        -------
        self.SCmatrixSelected

        """
        nCells = self.nCellsLoaded()
        SCmatrixSelected = np.zeros((self.numberBarcodes, self.numberBarcodes, nCells))

        if self.runParameters["action"] == "labeled":
            cellswithLabel = [idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x > 0]
            newCell = 0
            for iCell in cellswithLabel:
                SCmatrixSelected[:, :, newCell] = self.data["SCmatrixCollated"][:, :, iCell]
                newCell += 1
        elif self.runParameters["action"] == "unlabeled":
            cellswithLabel = [idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x == 0]
            newCell = 0
            for iCell in cellswithLabel:
                SCmatrixSelected[:, :, newCell] = self.data["SCmatrixCollated"][:, :, iCell]
                newCell += 1
        else:
            SCmatrixSelected = self.data["SCmatrixCollated"]
        print("nCells retrieved: {}".format(SCmatrixSelected.shape[2]))
        self.SCmatrixSelected = SCmatrixSelected


# =============================================================================
# FUNCTIONS
# =============================================================================


def normalizeProfile(profile1, profile2, runParameters):

    print("Normalization: {}".format(runParameters["normalize"]))

    mode = runParameters["normalize"]

    if "maximum" in mode:  # normalizes by maximum
        profile1 = profile1 / profile1.max() / 2
        profile2 = profile2 / profile2.max() / 2
    elif "none" in mode:  # no normalization
        m1_norm = 1
        m2_norm = 1
    else:  # normalizes by given factor
        normFactor = float(mode)
        profile1 = profile1 / 1
        profile2 = profile2 / normFactor

    return profile1, profile2


def plot1Dprofile2Datasets(ifigure, HiMdata1, HiMdata2, runParameters, anchor, iFigLabel, yticks, xticks, legend=False):

    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]
    lwbase = plt.rcParams["lines.linewidth"]
    thin, thick = lwbase / 2, lwbase * 3

    profile1 = HiMdata1.data["ensembleContactProbability"][:, anchor - 1]
    profile2 = HiMdata2.data["ensembleContactProbability"][:, anchor - 1]

    profile1, profile2 = normalizeProfile(profile1, profile2, runParameters)

    x = np.linspace(0, profile1.shape[0], num=profile1.shape[0], endpoint=True)
    tck1 = interpolate.splrep(x, profile1, s=0)
    tck2 = interpolate.splrep(x, profile2, s=0)
    xnew = np.linspace(0, profile1.shape[0], num=100, endpoint=True)
    ynew1 = interpolate.splev(xnew, tck1, der=0)
    ynew2 = interpolate.splev(xnew, tck2, der=0)
    if runParameters["splines"]:
        ifigure.plot(xnew, ynew1, "-", xnew, ynew2, "-")  # x, profile, 'o',
    else:
        ifigure.plot(x, profile1, "-", x, profile2, "-")  # x, profile, 'o',

    ifigure.set_xlim([0, profile1.shape[0]])
    ifigure.axvline(x=anchor - 0.5, color=colors[4], lw=thick, alpha=0.5)
    ifigure.set_ylim([0, runParameters["cAxis"]])

    if xticks:
        ifigure.set_xlabel("barcode #", fontsize=runParameters["fontsize"])
        if not runParameters["axisTicks"]:
            ifigure.set_xticklabels(())
        else:
            ifigure.set_xticklabels(HiMdata1.data["uniqueBarcodes"])
    else:
        ifigure.set_xticklabels(())

    if yticks:
        ifigure.set_ylabel("Probability", fontsize=runParameters["fontsize"])
        if not runParameters["axisTicks"]:
            ifigure.set_yticklabels(())
        else:
            ifigure.set_yticks([0, runParameters["cAxis"] / 2, runParameters["cAxis"]])
    else:
        ifigure.set_yticklabels(())

    if legend:
        ifigure.legend([HiMdata1.datasetName, HiMdata2.datasetName], loc="best")


def loadList(fileName):
    with open(fileName, newline="") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ", quotechar="|")
        runName = []
        for row in spamreader:
            # print(', '.join(row))
            if len(runName) > 0:
                runName.append(row)
            else:
                runName = row

    return runName


def attributesLabels2cells(SNDtable, ResultsTable, label="doc"):

    sortedSNDTable = SNDtable.group_by("MaskID #")
    listKeys = list(sortedSNDTable.groups.keys["MaskID #"].data)
    indexKey = [index for i, index in zip(listKeys, range(len(listKeys))) if i == label]

    # checks that there is at least one cell with the label
    if len(indexKey) > 0:

        SNDTablewithLabel = sortedSNDTable.groups[indexKey[0]]
        print("\n>>> Matching labels")
        print("Found {} out of {} cells with {} in dataset".format(len(SNDTablewithLabel), len(sortedSNDTable), label))

        # sorts Results Table by ROI
        PWDTableSortedROI = ResultsTable.group_by("ROI #")
        CUIDsList = []
        CUIDs = Table()
        CUIDs["Cuid"] = []

        print("ROIs to process: {}".format(PWDTableSortedROI.groups.keys))

        for ROI, group in zip(PWDTableSortedROI.groups.keys, PWDTableSortedROI.groups):
            # list of cellIDs in ROI
            cells2Process = group["CellID #"].data.compressed()
            cells2ProcessUID = group["Cuid"]

            # list of ROIs detected in SNDTablewithLabel
            ROIsinSNDTablewithLabel = list(SNDTablewithLabel.group_by("ROI #").groups.keys["ROI #"].data)

            # index of ROI within the keys of SNDTablewithLabel
            indexROIs = [
                index
                for i, index in zip(ROIsinSNDTablewithLabel, range(len(ROIsinSNDTablewithLabel)))
                if i == ROI["ROI #"]
            ]

            # subtable of cells with label and ROI that we are looking for
            SNDTablewithLabelROI = SNDTablewithLabel.group_by("ROI #").groups[indexROIs[0]]
            cellswithLabel = list(SNDTablewithLabelROI["CellID #"].data)

            # finds which cell indeces in Table have label
            listofSelectedCells = [
                index for iCell, index in zip(cells2Process, range(len(cells2Process))) if iCell in cellswithLabel
            ]

            if len(listofSelectedCells) > 0:
                print("Detected {} cells in ROI {} with label".format(len(listofSelectedCells), ROI["ROI #"]))
                if len(CUIDs) > 0:
                    # CUIDs = vstack([CUIDs, cells2ProcessUID[listofSelectedCells]])
                    # print('adding {} more cells'.format(len(cells2ProcessUID[listofSelectedCells])))
                    CUIDsList += list(cells2ProcessUID[listofSelectedCells].data.compressed())
                else:
                    CUIDsList = list(cells2ProcessUID[listofSelectedCells].data.compressed())
                    # CUIDs = cells2ProcessUID[listofSelectedCells]

            print(
                "Processed ROI # {}, found {} out of {} cells with {}".format(
                    ROI["ROI #"], len(listofSelectedCells), len(group), label
                )
            )

        # from list of CUIDs from cells that show label, I construct a binary vector of the same size as SCmatrix. Labeled cells have a 1.
        SClabeled = np.zeros(len(ResultsTable))
        # print('CUID list: {}'.format(CUIDsList2))
        # CUIDsList = CUIDs["Cuid"].data.compressed()
        # CUIDsList = CUIDsList2
        # checks that there are cells found with the label
        if len(CUIDsList) > 0:
            indexCellsWithLabel = [
                iRow for Row, iRow in zip(ResultsTable, range(len(ResultsTable))) if Row["Cuid"] in CUIDsList
            ]
            SClabeled[indexCellsWithLabel] = 1
        else:
            SClabeled, CUIDsList = [], []

        return SClabeled, CUIDsList
    else:
        # otherwise returns an empty list
        print("Warning: No cell with a mask labeled <{}> was found".format(label))
        return [], []


def loadsSCdata(ListData, datasetName, p):
    """
    loads SC datasets from a dict of folders (ListData)

    Parameters
    ----------
    ListData : dict
        dict of folders with data that can be loaded.
    dataset2Load : int, optional
        The item in the dictionary that will be loaded. The default is 3.

    Returns
    -------
    SCmatrixCollated : list of np arrays nBarcodes x nBarcodes x nCells
        Cummulative SC PWD matrix.
    uniqueBarcodes : list of np arrays
        containing the barcode identities for each matrix.
    buildsPWDmatrixCollated : list of Tables
        Tables with all the data for cells and barcodes used to produce SCmatrixCollated.

    """
    # tags2process = list(ListData.keys())
    print("Dataset to load: {}\n\n".format(list(ListData.keys())[0]))

    dimTag=""
    if "d3" in p.keys():
        if p["d3"]:
            dimTag="_3D"

    SCmatrixCollated, uniqueBarcodes = [], []
    buildsPWDmatrixCollated, runName, SClabeledCollated = [], [], []

    for rootFolder in ListData[datasetName]["Folders"]:
        # [finds sample name]
        runName.append(os.path.basename(os.path.dirname(rootFolder)))

        # [makes list of files with Tables to load]
        # tries to load files from newer version of proceesingPipeline.py
        files2Process = glob.glob(rootFolder + "/buildsPWDmatrix" + dimTag + "_order*ROI*.ecsv")
        if len(files2Process) == 0:
            # it resorts to old format
            files2Process = glob.glob(rootFolder + "/buildsPWDmatrix" + dimTag + "_*ROI*.ecsv")
        else:
            print("Found {} ECSV files in {}".format(len(files2Process), rootFolder))

        # checks that something was found
        if len(files2Process) > 0:
            print(">>> Loading {} results tables".format(len(files2Process)))

            # [initializes variables]
            buildsPWDmatrix = Table()
            fileOrder, fileOrderStamp, fileTimeStamp = (
                np.zeros(len(files2Process), dtype=int),
                np.zeros(len(files2Process), dtype=int),
                np.zeros(len(files2Process)),
            )

            # [finds what order Table files should be loaded to agree with order in buildsPWDmatrix_HiMscMatrix.npy]
            for fileName, ifileName in zip(files2Process, range(len(files2Process))):
                if "_order" in fileName:
                    for isplit in fileName.split("_"):
                        if "order" in isplit:
                            fileOrderStamp[ifileName] = int(isplit.split(":")[1])
                            print(
                                "order {}= {}--> {}".format(
                                    ifileName, os.path.basename(fileName), fileOrderStamp[ifileName]
                                )
                            )
                    fileTimeStamp[ifileName] = os.path.getmtime(fileName)
                    choosingTimeStamp = False

                else:
                    fileTimeStamp[ifileName] = os.path.getmtime(fileName)
                    choosingTimeStamp = True

            if choosingTimeStamp:
                fileOrder = np.argsort(fileTimeStamp).astype(int)
            else:
                fileOrder = np.argsort(fileOrderStamp).astype(int)

            # print('FileOrder: {}'.format(fileOrder))
            # [loads buildsPWDmatrix Tables]
            for ifileName in range(len(files2Process)):
                fileName = files2Process[fileOrder[ifileName]]
                newbuildsPWDmatrix = Table.read(fileName, format="ascii.ecsv")  # ascii.ecsv
                buildsPWDmatrix = vstack([buildsPWDmatrix, newbuildsPWDmatrix])
                print(
                    "[{}:{}:{}] From {}, Read: {} cells, Cummulative: {} cells".format(
                        ifileName,
                        fileOrder[ifileName],
                        fileTimeStamp[fileOrder[ifileName]],
                        os.path.basename(fileName),
                        len(newbuildsPWDmatrix),
                        len(buildsPWDmatrix),
                    )
                )

            # [loads SNDassignedCells.ecsv files if available]
            fileNameSNDassignedCells = os.path.dirname(rootFolder) + os.sep + "segmentedObjects/SNDassignedCells.ecsv"
            if os.path.exists(fileNameSNDassignedCells):
                print("Reading and processing: {}".format(fileNameSNDassignedCells))
                SNDassignedCells = Table.read(fileNameSNDassignedCells, format="ascii.ecsv")

                # checks that table is not empty
                if len(SNDassignedCells) > 0:
                    # attributes masks to single cells
                    SClabeled, CUIDsList = attributesLabels2cells(SNDassignedCells, buildsPWDmatrix, label=p["label"])

                    # checks that at least one cell was found to have the label
                    if len(SClabeled) == 0:
                        # if not available it makes a mock SClabeled matrix so that pipeline always works
                        # SClabeled = np.ones(len(buildsPWDmatrix)).astype(int)
                        SClabeled = np.zeros(len(buildsPWDmatrix)).astype(int)

            else:
                # if not available it makes a mock SClabeled matrix so that pipeline always works
                # SClabeled = np.ones(len(buildsPWDmatrix)).astype(int)
                SClabeled = np.zeros(len(buildsPWDmatrix)).astype(int)

            SClabeledCollated.append(SClabeled)

            # [loads and accumulates barcodes and scHiM matrix]
            fileNamMatrix = rootFolder + os.sep + "buildsPWDmatrix" + dimTag + "_HiMscMatrix.npy"
            fileNameBarcodes = rootFolder + os.sep + "buildsPWDmatrix" + dimTag + "_uniqueBarcodes.ecsv"

            if os.path.exists(fileNamMatrix):
                SCmatrix1 = np.load(fileNamMatrix)
                SCmatrixCollated.append(SCmatrix1)
            else:
                print("*** Error: could not find {}".format(fileNamMatrix))

            if os.path.exists(fileNameBarcodes):
                uniqueBarcodes.append(np.loadtxt(fileNameBarcodes).astype(int))
            else:
                print("*** Error: could not find {}".format(fileNameBarcodes))

            buildsPWDmatrixCollated.append(buildsPWDmatrix)

            print("\n>>>Merging rootFolder: {}".format(rootFolder))
            print("Cells added after merge: {}\n".format(SCmatrix1.shape[2]))
        else:
            print("No file detected in the folder you provide: {}".format(rootFolder + "/buildsPWDmatrix_*ROI*.ecsv"))
    print("{} datasets loaded\n".format(len(SCmatrixCollated)))

    return (
        SCmatrixCollated,
        uniqueBarcodes,
        buildsPWDmatrixCollated,
        runName,
        SClabeledCollated,
    )


def loadsSCdataMATLAB(ListData, datasetName, p):

    print("Dataset to load: {}\n\n".format(list(ListData.keys())[0]))

    SCmatrixCollated, uniqueBarcodes = [], []
    runName, SClabeledCollated = [], []

    for rootFolder in ListData[datasetName]["Folders"]:
        # [finds sample name]
        runName.append(os.path.basename(os.path.dirname(rootFolder)))

        # [loads and accumulates barcodes and scHiM matrix]
        fileNameMatrix = rootFolder + os.sep + "HiMscMatrix.mat"
        fileNameBarcodes = rootFolder + os.sep + "buildsPWDmatrix_uniqueBarcodes.ecsv"

        # loads barcodes
        if os.path.exists(fileNameBarcodes):
            uniqueBarcodes.append(np.loadtxt(fileNameBarcodes).astype(int))
            print(">>> Loaded {}".format(fileNameMatrix))
        else:
            print("*** Error: could not find {}".format(fileNameBarcodes))

        # loads SC matrix
        if os.path.exists(fileNameMatrix):
            data = loadmat(fileNameMatrix)
            SCmatrix1 = data["distanceMatrixCumulative"]
            # print(">>> SC matrix 1 shape: {}".format(SCmatrix1.shape))
            # SCmatrix2=SCmatrix1[uniqueBarcodes[0]-1,uniqueBarcodes[0]-1,:]
            SCmatrixCollated.append(SCmatrix1)
            print(">>> Loaded: {}\n SC matrix shape: {}".format(fileNameMatrix, SCmatrix1.shape))
        else:
            print("*** Error: could not find {}".format(fileNameMatrix))

        # loads cell attributes
        cellAttributesMatrix = data["cellAttributesMatrix"]
        ResultsTable = cellAttributesMatrix[0, :]

        SClabeled = np.zeros(len(ResultsTable))
        indexCellsWithLabel = [iRow for iRow, Row in enumerate(ResultsTable) if Row > 0]
        SClabeled[indexCellsWithLabel] = 1
        SClabeledCollated.append(SClabeled)

        print("\n>>>Merging rootFolder: {}".format(rootFolder))
        print("Cells added after merge: {}\n".format(SCmatrix1.shape[2]))

    print("{} datasets loaded\n".format(len(SCmatrixCollated)))

    return (
        SCmatrixCollated,
        uniqueBarcodes,
        runName,
        SClabeledCollated,
    )


def listsSCtoKeep(p, mask):
    #print("{}:{}".format(p["label"], p["action"]))
    if p["action"] == "all":
        try:
            cells2Plot = range(len(mask))
        except TypeError:
            print(mask)
            cells2Plot = range(mask.shape[0])
    elif p["action"] == "labeled":
        a = [i for i in range(len(mask)) if mask[i] == 1]
        cells2Plot = a
    elif p["action"] == "unlabeled":
        a = [i for i in range(len(mask)) if mask[i] == 0]
        cells2Plot = a

    print(">> label: {}\t action:{}\t Ncells2plot:{}\t Ncells in SCmatrix:{}".format(p["label"], p["action"],max(cells2Plot),len(mask)))
    
    return cells2Plot


def normalizeMatrix(SCmatrix_wt):
    SCmatrix_wt_normalized = SCmatrix_wt
    nBins = SCmatrix_wt.shape[0]

    for iRow in range(nBins):
        rowSum = np.sum(SCmatrix_wt[iRow, :])
        for iCol in range(nBins):
            SCmatrix_wt_normalized[iRow, iCol] = SCmatrix_wt_normalized[iRow, iCol] / rowSum
            SCmatrix_wt_normalized[iCol, iRow] = SCmatrix_wt_normalized[iCol, iRow] / rowSum
    return SCmatrix_wt_normalized


def plotsEnsemble3wayContactMatrix(
    SCmatrixCollated, uniqueBarcodes, anchors, sOut, runName, iListData, p, fileNameMD="tmp.md", datasetName=""
):

    # combines matrices from different samples and calculates integrated contact probability matrix
    SCmatrixAllDatasets = []  # np.zeros((nBarcodes,nBarcodes))
    for iSCmatrixCollated, iuniqueBarcodes, mask, iTag in zip(
        SCmatrixCollated, uniqueBarcodes, p["SClabeledCollated"], runName
    ):
        cells2Plot = listsSCtoKeep(p, mask)
        if len(cells2Plot) > 0:
            if max(cells2Plot) > iSCmatrixCollated.shape[2]:
                print(
                    "Error: max in cells2plot {} in dataset {} is larger than the number of available cells {}".format(
                        max(cells2Plot), iTag, iSCmatrixCollated.shape[2]
                    )
                )
            else:
                if len(SCmatrixAllDatasets) > 0:
                    SCmatrixAllDatasets = np.concatenate(
                        (SCmatrixAllDatasets, iSCmatrixCollated[:, :, cells2Plot]), axis=2
                    )
                else:
                    SCmatrixAllDatasets = iSCmatrixCollated[:, :, cells2Plot]
                commonSetUniqueBarcodes = iuniqueBarcodes
        else:
            print("Dataset: {} - {}  did not have any cell to plot".format(datasetName, iTag))

    # print(commonSetUniqueBarcodes)

    # loops over anchors
    for anchor in anchors:
        print("nCells processed: {}".format(SCmatrixAllDatasets.shape[2]))

        # calculates the 3-way matrix for a given anchor
        SCmatrix = calculate3wayContactMatrix(
            SCmatrixAllDatasets,
            uniqueBarcodes,
            p["pixelSize"],
            anchor,
            sOut,
            threshold=iListData["ContactProbability_distanceThreshold"],
            norm="nonNANs",
        )  # norm: nonNANs (default)

        # outputFileName = p['outputFolder'] + os.sep + datasetName + "_Cells:" + p['action'] + "_ensemble3wayContacts"
        outputFileName = (
            p["outputFolder"]
            + os.sep
            + datasetName
            + "_label:"
            + p["label"]
            + "_action:"
            + p["action"]
            + "_ensemble3wayContacts"
        )

        # outputs result
        outputFileName += "_anchor_" + str(anchor)
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")

        cScale = np.max(SCmatrix)
        # print(cScale)
        # print(SCmatrix.shape)

        # plots 3-way matrix
        plotMatrix(
            SCmatrix,
            commonSetUniqueBarcodes,  # before uniqueBarcodes
            p["pixelSize"],
            cm=iListData["ContactProbability_cm"],
            outputFileName=outputFileName,
            clim=cScale,
            cMin=iListData["ContactProbability_cmin"],
            figtitle="3way contacts",
            cmtitle=sOut,
            nCells=0,
            mode="counts",
        )  # twilight_shifted_r

        # saves matrices as individual files for further plotting
        rootOutputFileName = (
            p["outputFolder"]
            + os.sep
            + datasetName
            + "_label:"
            + p["label"]
            + "_action:"
            + p["action"]
            + "_anchor:"
            + str(anchor)
        )
        np.save(rootOutputFileName + "_ensemble3wayContacts.npy", SCmatrix)


def calculate3wayContactMatrix(
    iSCmatrixCollated, iuniqueBarcodes, pixelSize, anchor, sOut, threshold=0.25, norm="nonNANs"
):

    nX = nY = iSCmatrixCollated.shape[0]
    SCmatrix = np.zeros((nX, nY))

    # transform distance matrix from pixel to µm
    mat = pixelSize * iSCmatrixCollated

    # print(nX, nY)
    for bait1 in range(nX):
        for bait2 in range(nY):
            if bait1 == bait2:
                continue

            # print("current bait1", bait1, "bait2", bait2)
            n_contacts, n_nonNaN = getMultiContact(mat, anchor, bait1, bait2, threshold)
            if sOut == "Counts":
                SCmatrix[bait1, bait2] = n_contacts
            elif sOut == "Probability":
                SCmatrix[bait1, bait2] = n_contacts / n_nonNaN
            else:
                print("Unexpected sOut.")
                return -1

            # print(n_contacts / n_nonNaN)
            # print(type(n_contacts), type(n_nonNaN))

    SCmatrix[np.isnan(SCmatrix)] = 0  # set NaN to zero
    return SCmatrix


def getMultiContact(mat, anchor, bait1, bait2, threshold):
    """
    Input:
    mat        : pwd matrix, including only the bins of used RTs
    anchor     : anchor bin
    bait1      : bait bin #1
    bait2      : bait bin #2
    threshold  : contact threshold
    Output:
    n_contacts : number of contacts between bins anchor, bait1, and bait2
    n_nonNaN   : number of cells where the distances anchor-bait1 and anchor-bait2 are present
    """

    A = mat[anchor, bait1, :]
    B = mat[anchor, bait2, :]

    # get fraction of points in quadrant
    n1 = np.sum((A < threshold) & (B < threshold))
    totN = np.sum((~np.isnan(A)) & (~np.isnan(B)))

    return n1, totN


def plotsSinglePWDmatrices(
    SCmatrixCollated, uniqueBarcodes, runName, iListData, p, fileNameMD="tmp.md", datasetName="",
):
    # plots distance matrix for each dataset
    for iSCmatrixCollated, iuniqueBarcodes, iTag, mask in zip(
        SCmatrixCollated, uniqueBarcodes, runName, p["SClabeledCollated"]
    ):
        outputFileName = p["outputFolder"] + os.sep + iTag + "_Cells:" + p["action"] + "_PWDmatrix"

        # selects cels according to label
        cells2Plot = listsSCtoKeep(p, mask)

        plotMatrix(
            iSCmatrixCollated,
            iuniqueBarcodes,
            p["pixelSize"],
            outputFileName=outputFileName,
            figtitle="PWD:" + datasetName + iTag,
            cm=iListData["PWD_cm"],
            clim=iListData["PWD_clim"],
            mode=iListData["PWD_mode"],
            nCells=iSCmatrixCollated.shape[2],
            cells2Plot=cells2Plot,
        )  # twilight_shifted_r 1.4, mode: median KDE coolwarm terrain
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def plotsInversePWDmatrix(SCmatrixCollated, uniqueBarcodes, runName, iListData, p, fileNameMD, datasetName=""):
    # plots inverse distance matrix for each dataset
    for iSCmatrixCollated, iuniqueBarcodes, iTag, mask in zip(
        SCmatrixCollated, uniqueBarcodes, runName, p["SClabeledCollated"]
    ):
        outputFileName = p["outputFolder"] + os.sep + iTag + "_Cells:" + p["action"] + "_invPWDmatrix"

        # selects cels according to label
        cells2Plot = listsSCtoKeep(p, mask)
        # print('Dataset {} cells2plot: {}'.format(iTag,cells2Plot))

        plotMatrix(
            iSCmatrixCollated,
            iuniqueBarcodes,
            p["pixelSize"],
            cm=iListData["iPWD_cm"],
            outputFileName=outputFileName,
            clim=iListData["iPWD_clim"],
            mode=iListData["iPWD_mode"],
            figtitle="inverse PWD:" + datasetName + iTag,
            cmtitle="inverse distance, 1/nm",
            inverseMatrix=True,
            nCells=iSCmatrixCollated.shape[2],
            cells2Plot=cells2Plot,
        )  # twilight_shifted_r, mode: median KDE
        writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def plotsSingleContactProbabilityMatrix(
    SCmatrixCollated, uniqueBarcodes, runName, iListData, p, fileNameMD="tmp.md", datasetName="",
):
    # Plots contact probability matrices for each dataset
    if "minNumberContacts" in iListData.keys():
        minNumberContacts = iListData["minNumberContacts"]
    else:
        minNumberContacts = 0

    print("$ Min number contacts: {}".format(minNumberContacts))

    for iSCmatrixCollated, iuniqueBarcodes, iTag, mask in zip(
        SCmatrixCollated, uniqueBarcodes, runName, p["SClabeledCollated"]
    ):
        # selects cels according to label
        cells2Plot = listsSCtoKeep(p, mask)

        if not cells2Plot:
            break

        if max(cells2Plot) > iSCmatrixCollated.shape[2]:
            print(
                "Error with range in cells2plot {} as it is larger than the number of available cells {}".format(
                    max(cells2Plot), iSCmatrixCollated.shape[2]
                )
            )
        else:
            SCmatrix, nCells = calculateContactProbabilityMatrix(
                iSCmatrixCollated[:, :, cells2Plot],
                iuniqueBarcodes,
                p["pixelSize"],
                threshold=iListData["ContactProbability_distanceThreshold"],
                minNumberContacts = minNumberContacts,
                norm="nonNANs",
            )  # norm: nCells (default), nonNANs
            outputFileName = (
                p["outputFolder"] + os.sep + datasetName + iTag + "_Cells:" + p["action"] + "_contactProbability"
            )

            print("Dataset {} cells2plot: {}".format(iTag, cells2Plot))
            cScale = SCmatrix.max() / iListData["ContactProbability_scale"]

            plotMatrix(
                SCmatrix,
                iuniqueBarcodes,
                p["pixelSize"],
                cm=iListData["ContactProbability_cm"],
                outputFileName=outputFileName,
                cMin=iListData["ContactProbability_cmin"],
                clim=cScale,
                figtitle="HiM:" + datasetName + iTag,
                cmtitle="probability",
                nCells=nCells,
                cells2Plot=cells2Plot,
            )  # twilight_shifted_r terrain coolwarm
            writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")


def fusesSCmatrixCollatedFromDatasets(SCmatrixCollated, uniqueBarcodes, p, runName, iListData):
    # combines matrices from different embryos and calculates integrated contact probability matrix

    SCmatrixAllDatasets = []
    cells2Plot = []
    NcellsTotal = 0

    for iSCmatrixCollated, iuniqueBarcodes, mask, iTag in zip(
        SCmatrixCollated, uniqueBarcodes, p["SClabeledCollated"], runName
    ):
        NcellsTotal += mask.shape[0]
        # selects cels according to label
        cells2Plot = listsSCtoKeep(p, mask)

        if len(cells2Plot) > 0:
            if max(cells2Plot) > iSCmatrixCollated.shape[2]:
                print(
                    "Error: max in cells2plot {} in dataset {} is larger than the number of available cells {}".format(
                        max(cells2Plot), iTag, iSCmatrixCollated.shape[2]
                    )
                )
            else:
                if len(SCmatrixAllDatasets) > 0:
                    SCmatrixAllDatasets = np.concatenate(
                        (SCmatrixAllDatasets, iSCmatrixCollated[:, :, cells2Plot]), axis=2
                    )
                else:
                    SCmatrixAllDatasets = iSCmatrixCollated[:, :, cells2Plot]

                commonSetUniqueBarcodes = iuniqueBarcodes

    if p["saveMatrix"]:
        # write out the ensemble PWD map
        Nbarcodes = SCmatrixAllDatasets.shape[0]
        pixelSize = p["pixelSize"]
        meanSCmatrix = np.zeros((Nbarcodes, Nbarcodes))
        for bin1 in range(Nbarcodes):
            for bin2 in range(Nbarcodes):
                if bin1 != bin2:
                    (maximumKernelDistribution, _, _, _,) = distributionMaximumKernelDensityEstimation(
                        SCmatrixAllDatasets, bin1, bin2, pixelSize, optimizeKernelWidth=False
                    )
                    meanSCmatrix[bin1, bin2] = maximumKernelDistribution

        outputFileName = (
            p["outputFolder"]
            + os.sep
            + "CombinedMatrix_PWD_KDE"
            + ":"
            + list(iListData.keys())[0]
            + "_Cells:"
            + p["action"]
            + ".dat"
        )

        print(">>> Saving fused SCmatrix to {}".format(outputFileName))

        np.savetxt(
            outputFileName,
            meanSCmatrix,
            fmt="%.4f",
            delimiter=" ",
            newline="\n",
            header="Combined pairwise distance map (kernel density estimator)",
            footer="",
            comments="# ",
            encoding=None,
        )

    if p["getStructure"]:
        ## multi-dimensional scaling to get coordinates from PWDs
        # make sure meanSCmatrix is symmetric
        meanSCmatrix = 0.5 * (meanSCmatrix + np.transpose(meanSCmatrix))
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
        XYZ = mds.fit(meanSCmatrix).embedding_
        # print(XYZ.shape)
        print(XYZ)
        outputFileNamePDB = (
            p["outputFolder"]
            + os.sep
            + "CombinedMatrix_PWD_KDE"
            + ":"
            + list(iListData.keys())[0]
            + "_Cells:"
            + p["action"]
            + "_python.pdb"
        )
        write_XYZ_2_pdb(outputFileNamePDB, XYZ)

    return SCmatrixAllDatasets, commonSetUniqueBarcodes, cells2Plot, NcellsTotal


def plotsEnsembleContactProbabilityMatrix(
    SCmatrixCollated, uniqueBarcodes, runName, iListData, p, fileNameMD="tmp.md", datasetName="",
):

    if "minNumberContacts" in iListData.keys():
        minNumberContacts = iListData["minNumberContacts"]
    else:
        minNumberContacts = 0

    # combines matrices from different samples/datasers and calculates integrated contact probability matrix
    SCmatrixAllDatasets, commonSetUniqueBarcodes, cells2Plot, NcellsTotal = fusesSCmatrixCollatedFromDatasets(
        SCmatrixCollated, uniqueBarcodes, p, runName, iListData
    )

    print("nCells selected / processed: {}/{}".format(SCmatrixAllDatasets.shape[2], NcellsTotal))

    # calculates contact probability matrix from merged samples/datasets
    SCmatrix, nCells = calculateContactProbabilityMatrix(
        SCmatrixAllDatasets,
        commonSetUniqueBarcodes,
        p["pixelSize"],
        threshold=iListData["ContactProbability_distanceThreshold"],
        minNumberContacts = minNumberContacts,
        norm=p["HiMnormalization"],
    )  # norm: nCells (default), nonNANs

    # outputs line for MD file and sets output filename
    cScale = SCmatrix.max() / iListData["ContactProbability_scale"]
    outputFileName = p["outputFolder"] + os.sep + datasetName + "_Cells:" + p["action"] + "_ensembleContactProbability"
    writeString2File(fileNameMD, "![]({})\n".format(outputFileName + "_HiMmatrix.png"), "a")

    # plots results
    plotMatrix(
        SCmatrix,
        commonSetUniqueBarcodes,
        p["pixelSize"],
        cm=iListData["ContactProbability_cm"],
        outputFileName=outputFileName,
        clim=cScale,
        cMin=iListData["ContactProbability_cmin"],
        figtitle="HiM counts",
        cmtitle="probability",
        nCells=nCells,
    )  # twilight_shifted_r

    # saves SC matrix in text format
    np.savetxt(
        p["outputFolder"]
        + os.sep
        + "CombinedMatrix"
        + ":"
        + list(iListData.keys())[0]
        + "_Cells:"
        + p["action"]
        + ".dat",
        SCmatrix,
        fmt="%.4f",
        delimiter=" ",
        newline="\n",
        header="Combined contact probability matrix",
        footer="",
        comments="# ",
        encoding=None,
    )

    # saves barcodes in text format
    np.savetxt(
        p["outputFolder"]
        + os.sep
        + "UniqueBarcodes"
        + ":"
        + list(iListData.keys())[0]
        + "_Cells:"
        + p["action"]
        + ".dat",
        commonSetUniqueBarcodes,
        fmt="%.4f",
        delimiter=" ",
        newline="\n",
        header="unique barcodes",
        footer="",
        comments="# ",
        encoding=None,
    )

    return SCmatrix, commonSetUniqueBarcodes


def shuffleMatrix(matrix, index):

    newSize = len(index)
    newMatrix = np.zeros((newSize, newSize))

    if newSize <= matrix.shape[0]:
        for i in range(newSize):
            for j in range(newSize):
                if index[i] < matrix.shape[0] and index[j] < matrix.shape[0]:
                    newMatrix[i, j] = matrix[index[i], index[j]]
    else:
        print("Error: shuffle size {} is larger than matrix dimensions {}".format(newSize, matrix.shape[0]))
        print("Shuffle: {} ".format(index))

    return newMatrix


def comp_func(A, nw):
    n1 = A.shape[0]
    # th= 0;  # a threshold for a minimal distance above with calculate the signal
    # print "Size of the matrix entetered for the Compaction signal :"
    # print n1;

    somme_short = np.zeros((n1, 1))
    # somme_long = np.zeros((n1, 1));
    signal1 = np.zeros((n1, 1))
    n_int = np.zeros((n1, 1))
    # mat_test = np.zeros((n1, n1));

    for i in range(0, n1):
        if i <= nw:
            p1 = n1 + i - nw
            p2 = i + nw
            for k in range(0, n1):
                if k <= p2 or k >= p1:
                    somme_short[i] = somme_short[i] + A[i, k]
                    n_int[i] = n_int[i] + 1

        elif (n1 - nw) <= i:
            p1 = i - nw
            p2 = i + nw - n1
            for k in range(0, n1):
                if k <= p2 or k >= p1:
                    somme_short[i] = somme_short[i] + A[i, k]
                    n_int[i] = n_int[i] + 1

        else:
            p1 = i - nw
            p2 = i + nw
            for k in range(0, n1):
                if p1 <= k and k <= p2:
                    somme_short[i] = somme_short[i] + A[i, k]
                    n_int[i] = n_int[i] + 1

    signal1 = somme_short
    return signal1


def plotScalogram(matrix2plot, outputFileName=""):
    matrixSize = matrix2plot.shape[0]

    scales = range(0, 6)
    comp_scale1 = np.zeros((len(scales), matrixSize))
    for ii, nw in enumerate(scales):
        c = comp_func(matrix2plot, nw)
        comp_scale1[ii,] = c.T

    # definitions for the axes
    left, width = 0.065, 0.65
    bottom, height = 0.1, 0.8
    bottom_h = left_h = left + width + 0.02

    ax1_ = [0, 0.4, width, height]
    ax2_ = [left, bottom, width, 0.2]

    # fig = plt.figure()
    fig1 = plt.figure(constrained_layout=True)
    spec1 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig1)
    f1 = fig1.add_subplot(spec1[1, 0])  # 16
    # ax2 = plt.axes(ax2_)#,sharex=ax1)
    f1.contourf(
        comp_scale1,
        vmin=0.0,
        vmax=1.0,
        levels=[0, 0.15, 0.30, 0.45, 0.6, 0.75, 1.0],
        extent=[0, matrix2plot.shape[0], 0, 400],
        cmap="rainbow",
    )
    imColorBar = contourf(
        comp_scale1,
        vmin=0.0,
        vmax=1.0,
        levels=[0, 0.15, 0.30, 0.45, 0.6, 0.75, 1.0],
        extent=[0, matrix2plot.shape[0], 0, 400],
        cmap="rainbow",
    )
    f1.set_xlim([0, matrixSize])  # list(data['uniqueBarcodes']))
    f1.set_ylabel("Scales")
    bounds = [0.0, 0.5, 1.0]
    cb1 = colorbar(imColorBar, shrink=1, orientation="vertical", ticks=bounds, spacing="proportional", pad=0.04)

    # ax1 = plt.axes(ax1_)
    # pos=ax1.imshow(matrix2plot,cmap='coolwarm')
    # ax1.set_xlabel("barcode #")
    # ax1.set_ylabel("barcode #")
    # ax1.set_xlim([-0.5 , matrixSize-.5])#list(data['uniqueBarcodes']))
    # ax1.set_ylim([matrixSize-.5,-0.5])#list(data['uniqueBarcodes']))
    # cbar = plt.colorbar(pos, fraction=0.046, pad=0.04)

    # ax2 = plt.axes(ax2_)#,sharex=ax1)
    # ax2.contourf(  comp_scale1, vmin=0.0, vmax= 1.0,levels=[0, .15, .30, .45, 0.6,0.75, 1.0],extent=[0,matrix2plot.shape[0],0,400],cmap="rainbow");
    # imColorBar=contourf(  comp_scale1, vmin=0.0, vmax= 1.0,levels=[0, .15, .30, .45, 0.6,0.75, 1.0],extent=[0,matrix2plot.shape[0],0,400],cmap="rainbow");
    # ax2.set_xlim([-.5,matrixSize-.5])#list(data['uniqueBarcodes']))
    # ax2.set_ylabel("Scales");
    # bounds = [0.,0.5, 1.0];
    # cb1 = colorbar(imColorBar,shrink = 1, orientation="vertical",ticks=bounds, spacing='proportional',pad=0.04);

    if outputFileName:
        plt.savefig(outputFileName)
        print("Output scalogram: {}".format(outputFileName))


def write_XYZ_2_pdb(fileName, XYZ):
    # writes XYZ coordinates to a PDB file wth pseudoatoms
    # fileName : string of output file path, e.g. '/foo/bar/test2.pdb'
    # XYZ      : n-by-3 numpy array with atom coordinates

    n_atoms = XYZ.shape[0]
    with open(fileName, "w+") as fid:
        ## atom coordinates
        txt = "HETATM  {: 3d}  C{:02d} PSD P   1      {: 5.3f}  {: 5.3f}  {: 5.3f}  0.00  0.00      PSDO C  \n"
        for i in range(n_atoms):
            fid.write(txt.format(i + 1, i + 1, XYZ[i, 0], XYZ[i, 1], XYZ[i, 2]))

        ## connectivity
        txt1 = "CONECT  {: 3d}  {: 3d}\n"
        txt2 = "CONECT  {: 3d}  {: 3d}  {: 3d}\n"
        # first line of connectivity
        fid.write(txt1.format(1, 2))
        # consecutive lines
        for i in range(2, n_atoms):
            fid.write(txt2.format(i, i - 1, i + 1))
        # last line
        fid.write(txt1.format(i + 1, i))

        print("Done writing {:s} with {:d} atoms.".format(fileName, n_atoms))


def plotDistanceHistograms(
    SCmatrixCollated,
    pixelSize,
    outputFileName="test",
    logNameMD="log.md",
    mode="hist",
    limitNplots=10,
    kernelWidth=0.25,
    optimizeKernelWidth=False,
    maxDistance=4.0,
):

    if not isnotebook():
        NplotsX = NplotsY = SCmatrixCollated.shape[0]
    else:
        if limitNplots == 0:
            NplotsX = NplotsY = SCmatrixCollated.shape[0]
        else:
            NplotsX = NplotsY = min(
                [limitNplots, SCmatrixCollated.shape[0]]
            )  # sets a max of subplots if you are outputing to screen!

    bins = np.arange(0, maxDistance, 0.25)

    sizeX, sizeY = NplotsX * 4, NplotsY * 4

    fig, axs = plt.subplots(figsize=(sizeX, sizeY), ncols=NplotsX, nrows=NplotsY, sharex=True)

    for i in trange(NplotsX):
        for j in range(NplotsY):
            if i != j:
                # print('Printing [{}:{}]'.format(i,j))
                if mode == "hist":
                    axs[i, j].hist(pixelSize * SCmatrixCollated[i, j, :], bins=bins)
                else:
                    (maxKDE, distanceDistribution, KDE, x_d,) = distributionMaximumKernelDensityEstimation(
                        SCmatrixCollated,
                        i,
                        j,
                        pixelSize,
                        optimizeKernelWidth=optimizeKernelWidth,
                        kernelWidth=kernelWidth,
                        maxDistance=maxDistance,
                    )
                    axs[i, j].fill_between(x_d, KDE, alpha=0.5)
                    axs[i, j].plot(
                        distanceDistribution, np.full_like(distanceDistribution, -0.01), "|k", markeredgewidth=1,
                    )
                    axs[i, j].vlines(maxKDE, 0, KDE.max(), colors="r")

            axs[i, j].set_xlim(0, maxDistance)
            axs[i, j].set_yticklabels([])

    plt.xlabel("distance, um")
    plt.ylabel("counts")

    fileExtension = outputFileName.split(".")[-1]

    if len(fileExtension) == 3:
        fileName = outputFileName + "_PWDhistograms." + fileExtension
    else:
        fileName = outputFileName + "_PWDhistograms.png"

    print("Output figure: {}\n".format(fileName))
    plt.savefig(fileName)

    if not isnotebook():
        plt.close()

    writeString2File(logNameMD, "![]({})\n".format(outputFileName + "_PWDhistograms.png"), "a")


def plotMatrix(
    SCmatrixCollated,
    uniqueBarcodes,
    pixelSize,
    numberROIs=1,
    outputFileName="test",
    logNameMD="log.md",
    clim=1.4,
    cm="seismic",
    figtitle="PWD matrix",
    cmtitle="distance, um",
    nCells=0,
    mode="median",
    inverseMatrix=False,
    cMin=0,
    cells2Plot=[],
    fileNameEnding="_HiMmatrix.png",
):

    Nbarcodes = SCmatrixCollated.shape[0]

    ######################################################
    # Calculates ensemble matrix from single cell matrices
    ######################################################

    if len(SCmatrixCollated.shape) == 3:

        # matrix is 3D and needs combining SC matrices into an ensemble matrix
        if len(cells2Plot) == 0:
            cells2Plot = range(SCmatrixCollated.shape[2])

        meanSCmatrix, keepPlotting = calculatesEnsemblePWDmatrix(SCmatrixCollated, pixelSize, cells2Plot, mode=mode)

    else:

        # matrix is 2D and does not need further treatment
        if mode == "counts":
            meanSCmatrix = SCmatrixCollated
            keepPlotting = True
        else:
            meanSCmatrix = pixelSize * SCmatrixCollated
            keepPlotting = True

    if keepPlotting:
        # no errors occurred

        # Calculates the inverse distance matrix if requested in the argument.
        if inverseMatrix:
            meanSCmatrix = np.reciprocal(meanSCmatrix)

        # plots figure
        fig = plt.figure(figsize=(10, 10))
        pos = plt.imshow(meanSCmatrix, cmap=cm)  # colormaps RdBu seismic
        plt.xlabel("barcode #")
        plt.ylabel("barcode #")
        plt.title(
            figtitle
            + " | "
            + str(meanSCmatrix.shape[0])
            + " barcodes | n="
            + str(nCells)
            + " | ROIs="
            + str(numberROIs)
        )
        # print("matrix size: {} | barcodes:{}".format(SCmatrixCollated.shape[0],list(uniqueBarcodes)))
        plt.xticks(np.arange(SCmatrixCollated.shape[0]), uniqueBarcodes)
        plt.yticks(np.arange(SCmatrixCollated.shape[0]), uniqueBarcodes)
        cbar = plt.colorbar(pos, fraction=0.046, pad=0.04)
        cbar.minorticks_on()
        cbar.set_label(cmtitle)
        plt.clim(cMin, clim)

        if len(outputFileName.split(".")) > 1:
            if outputFileName.split(".")[1] != "png":
                if len(outputFileName.split(".")[1]) == 3:
                    # keeps original extension
                    o = outputFileName
                    plt.savefig(outputFileName)
                else:
                    # most likely the full filname contains other '.' in addition to that in the extension
                    o = outputFileName + fileNameEnding
                    plt.savefig(o)
            else:
                o = outputFileName.split(".")[0] + fileNameEnding
                plt.savefig(o)
        else:
            o = outputFileName + fileNameEnding
            plt.savefig(o)

        if not isnotebook():
            plt.close()
        if "png" not in o:
            o += ".png"
        writeString2File(logNameMD, "![]({})\n".format(o), "a")
    else:
        # errors during pre-processing
        print("Error plotting figure. Not executing script to avoid crash.")


def calculateContactProbabilityMatrix(iSCmatrixCollated, iuniqueBarcodes, pixelSize, threshold=0.25, norm="nCells", minNumberContacts = 0):

    nX = nY = iSCmatrixCollated.shape[0]
    nCells = iSCmatrixCollated.shape[2]
    SCmatrix = np.zeros((nX, nY))



    for i in range(nX):
        for j in range(nY):
            if i != j:
                distanceDistribution = pixelSize * iSCmatrixCollated[i, j, :]

                numberContacts = distanceDistribution.squeeze().shape[0]-len(np.nonzero(np.isnan(distanceDistribution))[0])

                if numberContacts < minNumberContacts:
                     print("$ Rejected {}-{} because number contacts: {} < {}".format(i,j,numberContacts,minNumberContacts))
                     probability = 0.0
                else:
                    # normalizes # of contacts by the # of cells
                    if norm == "nCells":
                        probability = len(np.nonzero(distanceDistribution < threshold)[0]) / nCells

                    # normalizes # of contacts by the # of PWD detected in each bin
                    elif norm == "nonNANs":
                        numberNANs = len(np.nonzero(np.isnan(distanceDistribution))[0])
                        if nCells == numberNANs:
                            probability = np.nan
                        else:
                            probability = len(np.nonzero(distanceDistribution < threshold)[0]) / (nCells - numberNANs)

                SCmatrix[i, j] = probability

    return SCmatrix, nCells


# @jit(nopython=True)
def findsOptimalKernelWidth(distanceDistribution):
    bandwidths = 10 ** np.linspace(-1, 1, 100)
    grid = GridSearchCV(KernelDensity(kernel="gaussian"), {"bandwidth": bandwidths}, cv=LeaveOneOut())
    grid.fit(distanceDistribution[:, None])
    return grid.best_params_


# @jit(nopython=True)
def retrieveKernelDensityEstimator(distanceDistribution0, x_d, optimizeKernelWidth=False, kernelWidth=0.25):
    """
    Gets the kernel density function and maximum from a distribution of PWD distances

    Parameters
    ----------
    distanceDistribution0 : nd array
        List of PWD distances.
    x_d : nd array
        x grid.
    optimizeKernelWidth : Boolean, optional
        whether to optimize bandwidth. The default is False.

    Returns
    -------
    np array
        KDE distribution.
    np array
        Original distribution without NaNs

    """

    nan_array = np.isnan(distanceDistribution0)

    not_nan_array = ~nan_array

    distanceDistribution = distanceDistribution0[not_nan_array]

    # instantiate and fit the KDE model
    if optimizeKernelWidth:
        kernelWidth = findsOptimalKernelWidth(distanceDistribution)["bandwidth"]
    else:
        kernelWidth = kernelWidth

    kde = KernelDensity(bandwidth=kernelWidth, kernel="gaussian")

    # makes sure the list is not full of NaNs.
    if distanceDistribution.shape[0] > 0:
        kde.fit(distanceDistribution[:, None])
    else:
        return np.array([0]), np.array([0])

    # score_samples returns the log of the probability density
    logprob = kde.score_samples(x_d[:, None])

    return logprob, distanceDistribution


# @jit(nopython=True)
def distributionMaximumKernelDensityEstimation(
    SCmatrixCollated, bin1, bin2, pixelSize, optimizeKernelWidth=False, kernelWidth=0.25, maxDistance=4.0
):
    """
    calculates the kernel distribution and its maximum from a set of PWD distances

    Parameters
    ----------
    SCmatrixCollated : np array 3 dims
        SC PWD matrix.
    bin1 : int
        first bin.
    bin2 : int
        first bin.
    pixelSize : float
        pixel size in um
    optimizeKernelWidth : Boolean, optional
        does kernel need optimization?. The default is False.

    Returns
    -------
    float
        maximum of kernel.
    np array
        list of PWD distances used.
    np array
        kernel distribution.
    x_d : np array
        x grid.

    """
    distanceDistributionUnlimited = pixelSize * SCmatrixCollated[bin1, bin2, :]  # full distribution
    distanceDistributionUnlimited = distanceDistributionUnlimited[
        ~np.isnan(distanceDistributionUnlimited)
    ]  # removes nans

    if bin1 == bin2:
        # protection agains bins in the diagonal
        distanceDistribution0 = distanceDistributionUnlimited
    else:
        # removes values larger than maxDistance
        distanceDistribution0 = distanceDistributionUnlimited[np.nonzero(distanceDistributionUnlimited < maxDistance)]
    x_d = np.linspace(0, maxDistance, 2000)

    # checks that distribution is not empty
    if distanceDistribution0.shape[0] > 0:
        logprob, distanceDistribution = retrieveKernelDensityEstimator(
            distanceDistribution0, x_d, optimizeKernelWidth, kernelWidth
        )
        if logprob.shape[0] > 1:
            kernelDistribution = 10 * np.exp(logprob)
            maximumKernelDistribution = x_d[np.argmax(kernelDistribution)]
            return maximumKernelDistribution, distanceDistribution, kernelDistribution, x_d
        else:
            return np.nan, np.zeros(x_d.shape[0]), np.zeros(x_d.shape[0]), x_d
    else:
        return np.nan, np.zeros(x_d.shape[0]), np.zeros(x_d.shape[0]), x_d


def getRgFromPWD(PWDmatrix0, minNumberPWD=4, threshold=6):
    """
    Calculates the Rg from a 2D pairwise distance matrix
    while taking into account that some of the PWD might be NaN

    PWDmatrix:       numpy array, NxN
    minFracNotNaN:   require a minimal fraction of PWDs to be not NaN, return NaN otherwise

    for the math, see https://en.wikipedia.org/wiki/Radius_of_gyration#Molecular_applications
    """

    PWDmatrix = PWDmatrix0.copy()

    # check that PWDmatrix is of right shape
    if PWDmatrix.ndim != 2:
        raise SystemExit("getRgFromPWD: Expected 2D input but got {}D.".format(PWDmatrix.ndim))
    if PWDmatrix.shape[0] != PWDmatrix.shape[1]:
        raise SystemExit("getRgFromPWD: Expected square matrix as input.")

    # make sure the diagonal is NaN
    np.fill_diagonal(PWDmatrix, np.NaN)

    # filters out PWD
    PWDmatrix[PWDmatrix > threshold] = np.nan

    # get the number of PWDs that are not NaN
    numPWDs = PWDmatrix.shape[0] * (PWDmatrix.shape[0] - 1) / 2
    numNotNan = np.sum(~np.isnan(PWDmatrix)) / 2  # default is to compute the sum of the flattened array

    if numNotNan < minNumberPWD:
        return np.NaN

    # calculate Rg
    sq = np.square(PWDmatrix)
    sq = np.nansum(sq)  # default is to compute the sum of the flattened array

    Rg_sq = sq / (2 * (2 * numNotNan + PWDmatrix.shape[0]))  # replaces 1/(2*N^2)

    Rg = np.sqrt(Rg_sq)

    return Rg


def getDetectionEffBarcodes(SCmatrixCollated):
    """
    Return the detection efficiency of all barcodes.
    Assumes a barcode is detected as soon as one PWD with this barcode is detected.
    """

    # check that PWDmatrix is of right shape
    if SCmatrixCollated.ndim != 3:
        raise SystemExit("getBarcodeEff: Expected 3D input but got {}D.".format(SCmatrixCollated.ndim))
    if SCmatrixCollated.shape[0] != SCmatrixCollated.shape[1]:
        raise SystemExit("getBarcodeEff: Expected axis 0 and 1 to have the same length.")

    # make sure the diagonal is NaN
    for i in range(SCmatrixCollated.shape[0]):
        SCmatrixCollated[i, i, :] = np.NaN

    # calculate barcode efficiency
    nCells = SCmatrixCollated.shape[2]

    eff = np.sum(~np.isnan(SCmatrixCollated), axis=0)
    eff[eff > 1] = 1

    eff0 = eff.copy()
    nCells2 = np.nonzero(np.sum(eff0, axis=0) > 2)[0].shape[0]

    eff = np.sum(eff, axis=-1)  # sum over all cells

    eff = eff / nCells2

    print("\n\n *** nCells={} | nCells2={}".format(nCells, nCells2))
    return eff


def getBarcodesPerCell(SCmatrixCollated):
    """
    Returns the number of barcodes that were detected in each cell of SCmatrixCollated.
    """

    # make sure the diagonal is NaN
    for i in range(SCmatrixCollated.shape[0]):
        SCmatrixCollated[i, i, :] = np.NaN

    numBarcodes = np.sum(~np.isnan(SCmatrixCollated), axis=0)
    numBarcodes[numBarcodes > 1] = 1
    numBarcodes = np.sum(numBarcodes, axis=0)

    return numBarcodes


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


def sortsCellsbyNumberPWD(HiMdata):

    # SCmatrix = HiMdata.data["SCmatrixCollated"]
    SCmatrix = HiMdata.SCmatrixSelected

    nCells = SCmatrix.shape[2]
    # print("Number of cells loaded: {}".format(nCells))

    # finds the number of barcodes detected per cell.
    nBarcodePerCell = list()
    values = list()
    dtype = [("cellID", int), ("nPWD", int)]

    for iCell in range(nCells):
        SCmatrixCell = SCmatrix[:, :, iCell]
        nPWD = int(np.count_nonzero(~np.isnan(SCmatrixCell)) / 2)
        nBarcodePerCell.append(nPWD)
        values.append((iCell, nPWD))

    valuesArray = np.array(values, dtype=dtype)  # create a structured array
    sortedValues = np.sort(valuesArray, order="nPWD")

    return SCmatrix, sortedValues, nCells


def kdeFit(x, x_d, bandwidth=0.2, kernel="gaussian"):

    kde = KernelDensity(bandwidth=bandwidth, kernel="gaussian")
    kde.fit(x[:, None])

    logprob = kde.score_samples(x_d[:, None])

    return logprob, kde


def calculatesEnsemblePWDmatrix(SCmatrix, pixelSize, cells2Plot, mode="median"):
    """
    performs a KDE or median to calculate the max of the PWD distribution

    Parameters
    ----------
    SCmatrix : TYPE
        DESCRIPTION.
    pixelSize : TYPE
        DESCRIPTION.

    Returns
    -------
    matrix = 2D npy array.

    """

    Nbarcodes = SCmatrix.shape[0]
    # cells2Plot = range(SCmatrix.shape[2])

    meanSCmatrix = np.zeros((Nbarcodes, Nbarcodes))

    if mode == "median":
        # calculates the median of all values #
        #######################################
        if max(cells2Plot) > SCmatrix.shape[2]:
            print(
                "Error with range in cells2plot {} as it is larger than the number of available cells {}".format(
                    max(cells2Plot), SCmatrix.shape[2]
                )
            )
            keepPlotting = False
        else:
            meanSCmatrix = pixelSize * np.nanmedian(SCmatrix[:, :, cells2Plot], axis=2)
            nCells = SCmatrix[:, :, cells2Plot].shape[2]
            keepPlotting = True

    elif mode == "KDE":
        keepPlotting = True

        if max(cells2Plot) > SCmatrix.shape[2]:
            print(
                "Error with range in cells2plot {} as it is larger than the number of available cells {}".format(
                    max(cells2Plot), SCmatrix.shape[2]
                )
            )
            keepPlotting = False
        else:
            for bin1 in trange(Nbarcodes):
                for bin2 in range(Nbarcodes):
                    if bin1 != bin2:
                        # print(f"cells2Plot:{cells2Plot}, ncells:{SCmatrix.shape}")
                        (maximumKernelDistribution, _, _, _,) = distributionMaximumKernelDensityEstimation(
                            SCmatrix[:, :, cells2Plot], bin1, bin2, pixelSize, optimizeKernelWidth=False,
                        )
                        meanSCmatrix[bin1, bin2] = maximumKernelDistribution

    return meanSCmatrix, keepPlotting
