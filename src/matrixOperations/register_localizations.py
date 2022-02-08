#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 15:00:35 2022

@author: marcnol

This class will handle correction of barcode positions from a table of local alignments

Remember that global alignments have already been corrected.

"""

# =============================================================================
# IMPORTS
# =============================================================================

import os, glob
import numpy as np
from tqdm import trange
from astropy.table import Table

from fileProcessing.fileManagement import (
    folders,
    writeString2File,
    printLog,
)

from imageProcessing.localization_table import localization_table

from matrixOperations.filter_localizations import get_file_table_new_name

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")

class register_localizations:
    def __init__(self, param):
        """
        Parameters
        ----------
        param : class
            Parameters
        """

        self.param = param
        self.alignmentResultsTableRead = Table()
        self.foundMatch = []

        if "toleranceDrift" in self.param.param["buildsPWDmatrix"]:
            self.toleranceDrift = self.param.param["buildsPWDmatrix"]["toleranceDrift"]
        else:
            self.toleranceDrift = 1
            printLog("# toleranceDrift not found. Set to {}!".format(self.toleranceDrift))

    def searchLocalShift(self, ROI, barcode, zxy_uncorrected):

        if self.alignmentResultsTableRead:
            return self.searchLocalShift_block3D(ROI, barcode, zxy_uncorrected)
        else: # no correction was applied because the localAlignmentTable was not found
            return zxy_uncorrected

    def searchLocalShift_block3D(self, ROI, barcode, zxy_uncorrected):
        """
        Searches for local drift for a specific barcode in a given ROI.
        If it exists then it adds to the uncorrected coordinates

        Parameters
        ----------
        ROI : string
            ROI used
        x_uncorrected : float
            x coordinate.
        y_uncorrected : float
            y coordinate.

        Returns
        -------
        x_corrected : float
            corrected x coordinate.
        y_corrected : float
            corrected y coordinate.

        """
        _foundMatch = False

        # gets blockSize
        blockSizeXY = self.alignmentResultsTable[0]["blockXY"]

        # zxy coord in block reference coord system
        zxyBlock = [np.floor(a/blockSizeXY).astype(int) for a in zxy_uncorrected]

        for row in self.alignmentResultsTable:

            # I need to check that the XY coordinates from localization are the same as the ij indices from the block decomposition!

            if row["ROI #"] == ROI and row["label"] == "RT" + str(barcode) and row["block_i"] == zxyBlock[1] and row["block_j"] == zxyBlock[2]:
                _foundMatch = True
                shifts = [row["shift_z"],row["shift_x"],row["shift_y"]]

                # checks that drifts > self.toleranceDrift are not applied
                if max(shifts) < self.toleranceDrift:
                    zxy_corrected = [a+shift for a,shift in zip(zxy_uncorrected,shifts)]
                else:
                    zxy_corrected = zxy_uncorrected

                # check for quality of shift correction before applying it !!
                #!TODO

        # keeps uncorrected values if no match is found
        if not _foundMatch:
            printLog("# Did not find match for ROI #{} barcode #{}".format(ROI, barcode))
            zxy_corrected = zxy_uncorrected
            self.foundMatch.append(False)
        else:
            self.foundMatch.append(True)

        return zxy_corrected

    def register_barcodes(self, barcodeMap):
        """
        This function will take a barcodeMap and a Table of 3D alignments to register barcode coordinates

        Returns
        -------
        None.

        """

        referenceFiducial = self.param.param["alignImages"]["referenceFiducial"]

        if "blockSize" in self.param.param["alignImages"]:
            blockSize = self.param.param["alignImages"]["blockSize"]
        else:
            blockSize = 256
            printLog("# blockSize not found. Set to {}!".format(blockSize))

        print(f"\n$ Reference barcode = {referenceFiducial }")

        NbarcodesROI = [], [], 0

        # loops over barcode Table rows in a given ROI
        for i in trange(len(barcodeMap.groups[0])): # i is the index of the barcode in barcodeMapROI
            barcode = barcodeMap.groups[0]["Barcode #"][i]
            ROI = barcodeMap.groups[0]["ROI #"][i]

            # keeps the particle if the test passed
            x_uncorrected = barcodeMap.groups[0]["ycentroid"][i]
            y_uncorrected = barcodeMap.groups[0]["xcentroid"][i]
            z_uncorrected = barcodeMap.groups[0]["zcentroid"][i]

            # Corrects XYZ coordinate of barcode if localDriftCorrection is available
            zxy_uncorrected = [z_uncorrected, x_uncorrected , y_uncorrected]
            RTbarcode = "RT" + str(barcode)

            if  RTbarcode not in self.param.param["alignImages"]["referenceFiducial"]:
                zxy_corrected = self.searchLocalShift(ROI, barcode, zxy_uncorrected)
            else:
                # if it is the reference cycle, then it does not correct coordinates
                zxy_corrected = zxy_uncorrected

            # rewrites corrected XYZ values to Table
            barcodeMap.groups[0]["ycentroid"][i] = zxy_corrected[1]
            barcodeMap.groups[0]["xcentroid"][i] = zxy_corrected[2]
            if self.ndims>2:
                barcodeMap.groups[0]["zcentroid"][i] = zxy_corrected[0]

        return barcodeMap

    def loadsLocalAlignment(self):

        if "None" in self.param.param["alignImages"]["localAlignment"]:
            printLog("\n\n$ localAlignment option set to {}".format(self.param.param["alignImages"]["localAlignment"]))
            return False, Table()
        else:
            return self._loadsLocalAlignment()

    def _loadsLocalAlignment(self):
        mode = self.param.param["alignImages"]["localAlignment"]
        localAlignmentFileName = self.dataFolder.outputFiles["alignImages"].split(".")[0] + "_" + mode + ".dat"

        if os.path.exists(localAlignmentFileName):
            alignmentResultsTable = Table.read(localAlignmentFileName, format="ascii.ecsv")
            alignmentResultsTableRead = True
            printLog("$ LocalAlignment file loaded: {}\n$ Will correct coordinates using {} alignment".format(localAlignmentFileName,mode))
            printLog("$ Number of records: {}".format(len(alignmentResultsTable)))
        else:
            printLog("\n\n# Warning: could not find localAlignment: {}\n Proceeding with only global alignments...".format(
                    localAlignmentFileName
                )
            )
            alignmentResultsTableRead = False
            alignmentResultsTable = Table()

        return alignmentResultsTable, alignmentResultsTableRead

    def register_barcodeMap_file(self, file):

        if "3D" in file:
            self.ndims = 3
        else:
            self.ndims = 2

        # loads barcode coordinate Tables
        table = localization_table()
        barcodeMapFull, uniqueBarcodes = table.load(file)

        # checks that everything is OK
        if len(barcodeMapFull) < 1:
            print(f"\nWARNING>{file} contains an empty table!")
            return None

        if 'comments' in barcodeMapFull.meta.keys():
            if "registered" in barcodeMapFull.meta['comments']:
                print(f"\nWARNING>{file} contains a table thas was already registered! \nWill not do anything")
                return None

        # preserves original copy of table for safe keeping
        new_file = get_file_table_new_name(file)
        table.save(new_file, barcodeMapFull)
        barcodeMapFull_unregistered = barcodeMapFull.copy()

        # indexes table by ROI
        barcodeMapROI, numberROIs = table.decode_ROIs(barcodeMapFull)

        for iROI in range(numberROIs):

            # creates sub Table for this ROI
            barcodeMap = barcodeMapROI.groups[iROI]
            nROI = barcodeMap['ROI #'][0]
            print(f"\nProcessing barcode localization Table for ROI: {nROI}")

            # registers localizations
            barcodeMap = self.register_barcodes(barcodeMap)

        # saves and plots registered barcode coordinate Tables
        table.save(file, barcodeMap, comments = 'registered')
        table.plots_distributionFluxes(barcodeMap, [file.split('.')[0], "_registered", "_barcode_stats", ".png"])
        table.plots_localizations(barcodeMap, [file.split('.')[0], "_registered", "_barcode_localizations",".png"])
        table.compares_localizations(barcodeMap,barcodeMapFull_unregistered,[file.split('.')[0], "_registered", "_barcode_comparisons",".png"])


    def register(self):

        """
        Function that registers barcodes using a local drift correction table produced by *alignImages3D*


        Returns
        -------
        None.

        """
        sessionName = "register_localizations"

        # processes folders and files
        self.dataFolder = folders(self.param.param["rootFolder"])
        printLog("\n===================={}====================\n".format(sessionName))
        printLog("$ folders read: {}".format(len(self.dataFolder.listFolders)))
        writeString2File(self.param.param["fileNameMD"], "## {}\n".format(sessionName), "a")
        label = 'barcode'

        currentFolder = self.param.param["rootFolder"]
        self.dataFolder.createsFolders(currentFolder, self.param)
        printLog("> Processing Folder: {}".format(currentFolder))

        # Loads localAlignment if it exists
        alignmentResultsTable, alignmentResultsTable_read_result = self.loadsLocalAlignment()

        if not alignmentResultsTable_read_result:
            print(f"Unable to find aligment table.\nDid you run alignImages3D?\n\n Aborted.")
            return

        # iterates over barcode localization tables in the current folder
        files = [x for x in glob.glob(self.dataFolder.outputFiles["segmentedObjects"] + "_*" + label + ".dat")]

        if len(files)<1:
            print("No file found to process!")
            return

        for file in files:

            self.register_barcodeMap_file(file)

        print(f" {len(files)} barcode tables processed in {currentFolder}")
