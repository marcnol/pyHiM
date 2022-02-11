#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 14:11:58 2022

@author: marcnol

This script will build chromatin traces using a segmentObjects_barcode table

The methods that will be implemented are:
    1= assigment by mask (either DAPI mask or other)
    2= spatial clusterization using KDtree. This method is mask-free.



Method 1:
    - iterates over ROIs
        - assigns barcode localizations to masks
        - applies local drift correction, if available
        - removes localizations using flux and driftTolerance
        - calculates the pair-wise distances for each single-cell mask
        - outputs are:
            - Table with #cell #PWD #coordinates (e.g. buildsPWDmatrix_3D_order:0_ROI:1.ecsv)
            - NPY array with single cell PWD single cell matrices (e.g. buildsPWDmatrix_3D_HiMscMatrix.npy)
            - NPY array with barcode identities (e.g. buildsPWDmatrix_3D_uniqueBarcodes.ecsv)
            - the files with no "3D" tag contain data analyzed using 2D localizations.

    - Single-cell results are combined together to calculate:
        - Distribution of pairwise distance for each barcode combination
        - Ensemble mean pairwise distance matrix using mean of distribution
        - Ensemble mean pairwise distance matrix using Kernel density estimation
        - Ensemble Hi-M matrix using a predefined threshold
        - For each of these files, there is an image in PNG format saved. Images containing "3D" are for 3D other are for 2D.


"""

# =============================================================================
# IMPORTS
# =============================================================================

import glob, os, sys
import uuid
import re
import numpy as np
from tqdm.contrib import tzip
from tqdm import trange
import matplotlib.pyplot as plt
from skimage.segmentation import expand_labels

from sklearn.metrics import pairwise_distances
from astropy.table import Table
from photutils.segmentation import SegmentationImage

from fileProcessing.fileManagement import (
    folders,
    writeString2File,
    printLog,
    getDictionaryValue,
)
from matrixOperations.HIMmatrixOperations import plotMatrix, plotDistanceHistograms, calculateContactProbabilityMatrix
from imageProcessing.localization_table import localization_table
from matrixOperations.chromatin_trace_table import chromatin_trace_table

from matrixOperations.filter_localizations import get_file_table_new_name

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")

# =============================================================================
# CLASSES
# =============================================================================


class build_traces:
    def __init__(self, param):
        self.param = param

        self.initialize_parameters()

        # initialize with default values
        self.currentFolder = []
        self.maskIdentifier = ['DAPI'] # default mask label

    def initializes_masks(self, Masks):
        self.Masks = Masks
        self.NcellsAssigned = 0
        self.NcellsUnAssigned = 0
        self.NbarcodesinMask = 0
        self.SegmentationMask = SegmentationImage(self.Masks)
        self.numberMasks = self.SegmentationMask.nlabels
        self.barcodesinMask = dict()

        for mask in range(self.numberMasks + 1):
            self.barcodesinMask["maskID_" + str(mask)] = []

    def initialize_parameters(self):
        # initializes parameters from param

        self.tracing_method = getDictionaryValue(self.param.param["buildsPWDmatrix"], "tracing_method",  default="masking")
        self.zBinning = getDictionaryValue(self.param.param["acquisition"], "zBinning", default=1)
        self.pixelSizeXY = getDictionaryValue(self.param.param["acquisition"], "pixelSizeXY", default=0.1)
        self.pixelSizeZ_0 = getDictionaryValue(self.param.param["acquisition"], "pixelSizeZ", default=0.25)
        self.pixelSizeZ = self.zBinning * self.pixelSizeZ_0
        self.pixelSize = [self.pixelSizeXY, self.pixelSizeXY, self.pixelSizeZ]
        self.availableMasks = getDictionaryValue(self.param.param["buildsPWDmatrix"], "masks2process",  default={"nuclei":"DAPI"})
        self.logNameMD = self.param.param["fileNameMD"]
        self.mask_expansion = getDictionaryValue(self.param.param["buildsPWDmatrix"], "mask_expansion", default=8)
        self.availableMasks = self.param.param["buildsPWDmatrix"]["masks2process"]

    def initializeLists(self):
        self.ROIs, self.cellID, self.nBarcodes, self.barcodeIDs, self.cuid, self.buid, self.barcodeCoordinates = (
            [],
            [],
            [],
            [],
            [],
            [],
            [],
        )

    def alignByMasking(self):
        """
        Assigns barcodes to masks and creates <NbarcodesinMask>
        And by filling in the "Cell #" key of barcodeMapROI
        This routine will only select which barcodes go to each cell mask

        Returns
        -------
        self.barcodesinMask # dictionnary with the identities of barcodes contained in each mask.
            Keys: 'maskID_1', 'maskID_2', and so on

        self.NbarcodesinMask # vector containing the number of barcodes for each mask
        self.NcellsAssigned # number of cells assigned
        self.NcellsUnAssigned # number of cells unassigned
        """

        NbarcodesinMask = np.zeros(self.numberMasks + 2)
        NbarcodesROI = 0

        # loops over barcode Table rows in a given ROI
        printLog("> Aligning by masking...")
        for i in trange(len(self.barcodeMapROI.groups[0])):  # i is the index of the barcode in barcodeMapROI
            barcode = self.barcodeMapROI.groups[0]["Barcode #"][i]

            # gets xyz coordinates
            x_corrected = self.barcodeMapROI.groups[0]["ycentroid"][i]
            y_corrected = self.barcodeMapROI.groups[0]["xcentroid"][i]

            if self.ndims == 2:
                z_corrected = self.barcodeMapROI.groups[0]["zcentroid"][i] = 0.0
            else:
                z_corrected = self.barcodeMapROI.groups[0]["zcentroid"][i]

            # binarizes coordinate
            y_int = int(y_corrected)
            x_int = int(x_corrected)

            # finds what mask label this barcode is sitting on
            maskID = self.Masks[x_int][y_int]

            # attributes CellID to a barcode
            self.barcodeMapROI["CellID #"][i] = maskID

            # if it is not background,
            if maskID > 0:
                # increments counter of number of barcodes in the cell mask attributed
                NbarcodesinMask[maskID] += 1

                # stores the identify of the barcode in a mask dictionary
                self.barcodesinMask["maskID_" + str(maskID)].append(i)

                # keeps statistics
                if int(self.barcodeMapROI.groups[0]["ROI #"][i]) == int(self.nROI):
                    NbarcodesROI += 1

        # Total number of masks assigned and not assigned
        self.NcellsAssigned = np.count_nonzero(NbarcodesinMask > 0)
        self.NcellsUnAssigned = self.numberMasks - self.NcellsAssigned

        # this list contains which barcodes are allocated to which masks
        self.NbarcodesinMask = NbarcodesinMask

        printLog("$ Number of cells assigned: {} | discarded: {}".format(self.NcellsAssigned, self.NcellsUnAssigned))

    def buildsVector(self, groupKeys, x, y, z):
        """
        Builds vector from coordinates

        Parameters
        ----------
        groupKeys : list
            list of headers in the barcodes table
        x : float
            x coordinates
        y : float
            y coordinates
        z : float
            z coordinates

        Returns
        -------
        R : np array
            vector with coordinates in nanometers.

        """

        R = np.column_stack((x * self.pixelSize["x"], y * self.pixelSize["y"], z * self.pixelSize["z"]))

        return R

    def buildsSCdistanceTable(self):
        """
        iterates over all masks, calculates PWD for each mask, assigns them to SCdistanceTable

        Returns
        -------
        SCdistanceTable

        """
        # sorts Table by cellID
        barcodeMapROI = self.barcodeMapROI

        # indexes table by cellID
        barcodeMapROI_cellID = barcodeMapROI.group_by("CellID #")

        self.initializeLists()

        # iterates over all traces in an ROI
        printLog("> Building single traces")
        for key, group in tzip(barcodeMapROI_cellID.groups.keys, barcodeMapROI_cellID.groups):
            if key["CellID #"] > 1:  # excludes trace 0 as this is background

                groupKeys, CellID, ROI = group.keys(), key["CellID #"], group["ROI #"].data[0]
                trace_ID = str(uuid.uuid4())

                # gets lists of x, y and z coordinates for barcodes assigned to a cell mask
                x, y, z = (
                    np.array(group["xcentroid"].data),
                    np.array(group["ycentroid"].data),
                    np.array(group["zcentroid"].data),
                )

                # gets vector in nm
                R_mum = self.buildsVector(groupKeys, x, y, z)

                for i in range(x.shape[0]):
                    entry = [group["Buid"].data[i], # spot uid
                            trace_ID,     # trace uid
                            R_mum[i][0],            # x, microns
                            R_mum[i][1],            # y, microns
                            R_mum[i][2],            # z, microns
                            "Chr3R",               # chromosome
                            129999,                # start sequence
                            149999,                # end sequence
                            ROI,                   # ROI number
                            CellID,                # Mask number
                            group["Barcode #"].data[i], # Barcode name
                            ]
                    self.trace_table.data.add_row(entry)

        printLog("$ Coordinates dimensions: {}".format(self.ndims))

    def load_mask(self,TIF_files_in_folder,):

        # finds files with cell masks
        channel = self.param.param["acquisition"][self.maskType+"_channel"]

        fileList2Process = [
            file
            for file in TIF_files_in_folder
            if self.param.decodesFileParts(file)["channel"] == channel # typically "ch00"
            and self.maskIdentifier in os.path.basename(file).split("_")
            and int(self.param.decodesFileParts(file)["roi"]) == self.nROI
        ]

        if len(fileList2Process) > 0:

            # loads file with cell masks
            fileNameROImasks = os.path.basename(fileList2Process[0]).split(".")[0] + "_Masks.npy"
            fullFileNameROImasks = self.dataFolder.outputFolders["segmentedObjects"] + os.sep + fileNameROImasks

            if os.path.exists(fullFileNameROImasks):

                # loads and initializes masks
                segmented_masks = np.load(fullFileNameROImasks)

                # expands mask without overlap by a maximmum of 'distance' pixels
                self.Masks= expand_labels(segmented_masks, distance = self.mask_expansion)

                # initializes masks
                self.initializes_masks(self.Masks)
                return True

            else:
                # Could not find a file with masks to assign. Report and continue with next ROI
                debug_mask_fileName(TIF_files_in_folder,fullFileNameROImasks,self.maskIdentifier,self.nROI,label=self.param.param["acquisition"]["label_channel"])

        else:
            printLog(f"$ Did not identified any filename for mask: {self.maskIdentifier}, channel: {channel}","WARN")
            printLog("-"*80)

        return False

    def assign_masks(self,
                     outputFileName,
                     barcodeMap,
                 ):
        """
        Main function that:
            loads and processes barcode localization files, local alignment file, and masks
            initializes <cellROI> class and assigns barcode localizations to masks
            then constructs the single cell PWD matrix and outputs it toghether with the contact map and the N-map.

        Parameters
        ----------
        outputFileName : string
        self.param : Parameters Class
        self.currentFolder : string
        self.dataFolder : Folder Class
            information to find barcode localizations, local drift corrections and masks

        self.pixelSize : dict, optional
            pixelSize = {'x': pixelSizeXY,
                        'y': pixelSizeXY,
                        'z': pixelSizeZ}
            The default is 0.1 for x and y, 0.0 for z. Pixelsize in um

        self.logNameMD : str, optional
            Filename of Markdown output. The default is "log.md".
        self.ndims : int, optional
            indicates whether barcodes were localized in 2 or 3D. The default is 2.
        self.maskIdentifier:

        Returns
        -------
        None.

        """

        # indexes localization tables by ROI
        barcodeMapROI = barcodeMap.group_by("ROI #")
        numberROIs = len(barcodeMapROI.groups.keys)

        printLog("-"*80)
        printLog(" Loading masks and pre-processing barcodes for Mask <{}> for {} ROIs".format(self.maskIdentifier, numberROIs))

        # finds TIFs in currentFolder
        TIF_files_in_folder = glob.glob(self.currentFolder + os.sep + "*.tif")

        # loops over ROIs
        processingOrder = 0

        for ROI in range(numberROIs):
            self.nROI = barcodeMapROI.groups.keys[ROI][0]  # need to iterate over the first index
            self.barcodeMapROI = barcodeMap.group_by("ROI #").groups[ROI]

            mask_loaded = self.load_mask(TIF_files_in_folder,)
            if mask_loaded:

                printLog("> Processing ROI# {}".format(self.nROI))

                # initializes trace table
                self.trace_table.initialize()

                # finds what barcodes are in each cell mask
                self.alignByMasking()
                printLog(f"$ ROI: {ROI}, N cells assigned: {self.NcellsAssigned - 1} out of {self.numberMasks}\n")

                # builds SCdistanceTable
                self.buildsSCdistanceTable()
                printLog("$ Number of entries in trace table: {}".format(len(self.trace_table.data)))

                # saves trace table with results per ROI
                output_table_fileName = outputFileName + "_mask:" + str(self.maskIdentifier) + "_ROI:" + str(self.nROI) + ".ecsv"
                self.trace_table.save(output_table_fileName, self.trace_table.data)

                # plots results
                self.trace_table.plots_traces([output_table_fileName.split(".")[0], "_traces_XYZ", ".png"], Masks = self.Masks)

                printLog(f"$ Saved output table as {output_table_fileName}")
                processingOrder += 1

    def build_trace_by_masking(self, barcodeMap):

        printLog("> Masks labels: {}".format(self.availableMasks))

        for maskLabel in self.availableMasks.keys():

            self.maskIdentifier = self.availableMasks[maskLabel]
            if "DAPI" in self.maskIdentifier:
                self.maskType = "DAPI"
            else:
                self.maskType = "mask"

            tag = '_' + str(self.ndims) + 'D'

            outputFileName = self.dataFolder.outputFiles["buildsPWDmatrix"] + tag

            # creates and initializes trace table
            self.trace_table = chromatin_trace_table()

            self.assign_masks(
                outputFileName,
                barcodeMap,
            )

            printLog("$ Trace built using mask assignment. Output saved in: {} ".format(self.currentFolder), "info")

    def launch_analysis(self,file):

        # loads barcode coordinate Tables
        table = localization_table()
        barcodeMap, self.uniqueBarcodes = table.load(file)

        if "3D" in file:
            self.ndims = 3
            self.pixelSize = {"x": self.pixelSizeXY, "y": self.pixelSizeXY, "z": self.pixelSizeZ}
        else:
            self.ndims = 2
            self.pixelSize = {"x": self.pixelSizeXY, "y": self.pixelSizeXY, "z": 0}

        if "masking" in self.tracing_method:
            self.build_trace_by_masking(barcodeMap)

    def run(self):
        """
        Function that assigns barcode localizations to masks and constructs single cell cummulative PWD matrix.

        Parameters
        ----------
        param : class
            Parameters
        log1 : class
            logging class.
        session1 : class
            session information

        Returns
        -------
        None.

        """
        # initializes sessionName, dataFolder, currentFolder
        label = "barcode"
        self.dataFolder, self.currentFolder  = initialize_module(self.param, module_name="build_traces",label = label)

        printLog("> Masks labels: {}".format(self.availableMasks))

        # iterates over barcode localization tables in the current folder
        files = [x for x in glob.glob(self.dataFolder.outputFiles["segmentedObjects"] + "_*" + label + ".dat")]

        if len(files) < 1:
            printLog("$ No localization table found to process!","WARN")
            return

        for file in files:
            self.launch_analysis(file)

        printLog(f"$ {len(files)} barcode tables processed in {self.currentFolder}")



def initialize_module(param, module_name="build_traces",label = "barcode"):

    sessionName = module_name

    # processes folders and files
    dataFolder = folders(param.param["rootFolder"])
    printLog("\n"+"="*35+f"{sessionName}"+"="*35+"\n")
    printLog("$ folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(param.param["fileNameMD"], "## {}\n".format(sessionName), "a")

    currentFolder = param.param["rootFolder"]
    dataFolder.createsFolders(currentFolder, param)
    printLog("> Processing Folder: {}".format(currentFolder))

    return dataFolder, currentFolder

def debug_mask_fileName(filesinFolder,fullFileNameROImasks,maskIdentifier,nROI,label=''):

    printLog(f"# Error, no mask file found for ROI: {nROI}\n")
    printLog("# File I was searching for: {}".format(fullFileNameROImasks))
    printLog("# Debug: ")
    for file in filesinFolder:
        if (
            file.split("_")[-1].split(".")[0]
            == label  # typically "ch00"
            and maskIdentifier in file.split("_")
            and int(os.path.basename(file).split("_")[3]) == nROI
        ):
            printLog("$ Hit found!")
        printLog(
            "fileSplit:{}, {} in filename: {}, ROI: {}".format(
                file.split("_")[-1].split(".")[0],
                maskIdentifier,
                maskIdentifier in os.path.basename(file).split("_"),
                int(os.path.basename(file).split("_")[3]),
            )
        )
