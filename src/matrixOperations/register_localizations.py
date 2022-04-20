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

from imageProcessing.localization_table import LocalizationTable

from matrixOperations.filter_localizations import get_file_table_new_name

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")


class RegisterLocalizations:
    def __init__(self, param):
        """
        Parameters
        ----------
        param : class
            Parameters
        """

        self.param = param
        self.alignmentResultsTableRead = False
        self.foundMatch = []

        if "toleranceDrift" in self.param.param["buildsPWDmatrix"]:
            self.toleranceDrift = self.param.param["buildsPWDmatrix"]["toleranceDrift"]
        else:
            self.toleranceDrift = 1
            printLog("# toleranceDrift not found. Set to {}!".format(self.toleranceDrift))

        if "remove_uncorrected_localizations" in self.param.param["buildsPWDmatrix"]:
            self.remove_uncorrected_localizations = self.param.param["buildsPWDmatrix"]["remove_uncorrected_localizations"]
        else:
            self.remove_uncorrected_localizations = True

        if self.remove_uncorrected_localizations:
            printLog("# Uncorrected localizations will be removed!!")
        else:
            printLog("# Uncorrected localizations will be kept!!")                


    def searchLocalShift(self, ROI, barcode, zxy_uncorrected):

        if self.alignmentResultsTableRead:
            return self.searchLocalShift_block3D(ROI, barcode, zxy_uncorrected)
        else:  # no correction was applied because the localAlignmentTable was not found
            return zxy_uncorrected, {'below_tolerance':False}
            print("ERROR> did not found alignmentResultsTable")

    def searchLocalShift_block3D(self, ROI, barcode, zxy_uncorrected):
        """
        Searches for local drift for a specific barcode in a given ROI.
        If it exists then it adds to the uncorrected coordinates

        Parameters
        ----------
        ROI : int
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

        # gets blockSize
        blockSizeXY = self.alignmentResultsTable[0]["blockXY"]

        # zxy coord in block reference coord system
        zxyBlock = [np.floor(a / blockSizeXY).astype(int) for a in zxy_uncorrected]

        # defines ROI, barcode ID, and blocks for the localization
        nROI = "ROI:" + str(ROI)
        nBarcode = "barcode:" + "RT" + str(barcode)
        nBlock_i = "block_i:" + str(zxyBlock[1])
        nBlock_j = "block_j:" + str(zxyBlock[2])

        # finds the corresponding shift int the dictionary

        shifts = [self.dictErrorBlockMasks[nROI][nBarcode][nBlock_i][nBlock_j]["shift_z"],
                  self.dictErrorBlockMasks[nROI][nBarcode][nBlock_i][nBlock_j]["shift_x"],
                  self.dictErrorBlockMasks[nROI][nBarcode][nBlock_i][nBlock_j]["shift_y"]]

        if max(np.abs(shifts)) < self.toleranceDrift:
            zxy_corrected = [a + shift for a, shift in zip(zxy_uncorrected, shifts)]
            quality_correction = {'below_tolerance':True}
        else:
            zxy_corrected = zxy_uncorrected
            quality_correction = {'below_tolerance':False}

        return zxy_corrected, quality_correction 

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

        printLog(
            f"\n$ Parameters:\n Blocksize = {blockSize}\n Tolerance = {self.toleranceDrift}\n Reference barcode = {referenceFiducial}"
        )

        NbarcodesROI = [], [], 0
        list_uncorrected_barcodes = list()

        # loops over barcode Table rows in a given ROI
        for i in trange(len(barcodeMap.groups[0])):  # i is the index of the barcode in barcodeMapROI
            barcode = barcodeMap.groups[0]["Barcode #"][i]
            ROI = barcodeMap.groups[0]["ROI #"][i]

            # keeps the particle if the test passed
            x_uncorrected = barcodeMap.groups[0]["ycentroid"][i]
            y_uncorrected = barcodeMap.groups[0]["xcentroid"][i]
            z_uncorrected = barcodeMap.groups[0]["zcentroid"][i]

            # Corrects XYZ coordinate of barcode if localDriftCorrection is available
            zxy_uncorrected = [z_uncorrected, x_uncorrected, y_uncorrected]
            RTbarcode = "RT" + str(barcode)

            if RTbarcode not in self.param.param["alignImages"]["referenceFiducial"]:
                zxy_corrected, quality_correction = self.searchLocalShift(ROI, barcode, zxy_uncorrected)
            else:
                # if it is the reference cycle, then it does not correct coordinates
                zxy_corrected = zxy_uncorrected

            if not quality_correction['below_tolerance']:
                list_uncorrected_barcodes.append(i) 

                if self.remove_uncorrected_localizations:
                    # will remove localizations that cannot be corrected
                    zxy_corrected = [np.nan, np.nan, np.nan]
                else:
                    # will keep uncorrected localizations
                    pass
                
            # rewrites corrected XYZ values to Table
            barcodeMap.groups[0]["ycentroid"][i] = zxy_corrected[1]
            barcodeMap.groups[0]["xcentroid"][i] = zxy_corrected[2]
            if self.ndims > 2:
                barcodeMap.groups[0]["zcentroid"][i] = zxy_corrected[0]

        if self.remove_uncorrected_localizations:
            printLog(f"$ {len(list_uncorrected_barcodes)} localizations out of {len(barcodeMap.groups[0])} were removed.")
        else:
            printLog(f"$ {len(list_uncorrected_barcodes)} localizations out of {len(barcodeMap.groups[0])} were uncorrected.")            

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
            self.alignmentResultsTable = Table.read(localAlignmentFileName, format="ascii.ecsv")
            self.alignmentResultsTableRead = True

            # builds dict of local alignments
            self.build_local_alignment_dict()

            printLog(
                "$ LocalAlignment file loaded: {}\n$ Will correct coordinates using {} alignment".format(
                    localAlignmentFileName, mode
                )
            )
            printLog("$ Number of records: {}".format(len(self.alignmentResultsTable)))
        else:
            printLog(
                "\n\n# Warning: could not find localAlignment: {}\n Proceeding with only global alignments...".format(
                    localAlignmentFileName
                )
            )
            self.alignmentResultsTableRead = False
            self.alignmentResultsTable = Table()
            self.dictErrorBlockMasks = Table()

    def build_local_alignment_dict(self):
        """
        Builds dictionary of local corrections for each ROI, barcode cycle, and block combination

        Parameters
        ----------
        self.alignmentResultsTable: astropy Table
            alignmentResultsTable table

        self.alignmentResultsTableRead: Boolean
            True when alignmentResultsTable table was read from disk

        Returns
        -------
        exit_code: Boolean

        self.dictErrorBlockMasks: dict

        """
        if not self.alignmentResultsTableRead:
            printLog("Did not find alignmentResultsTable. Cannot continue")
            return False
        else:
            alignmentResultsTable = self.alignmentResultsTable

        # gets blockSize
        blockSizeXY = alignmentResultsTable[0]["blockXY"]

        dictErrorBlockMasks = dict()

        for row in alignmentResultsTable:
            nROI = "ROI:" + str(row["ROI #"])
            nBarcode = "barcode:" + row["label"]
            nBlock_i = "block_i:" + str(row["block_i"])
            nBlock_j = "block_j:" + str(row["block_j"])

            if nROI not in dictErrorBlockMasks.keys():
                dictErrorBlockMasks[nROI] = {}

            if nBarcode not in dictErrorBlockMasks[nROI].keys():
                dictErrorBlockMasks[nROI][nBarcode] = {}

            if nBlock_i not in dictErrorBlockMasks[nROI][nBarcode].keys():
                dictErrorBlockMasks[nROI][nBarcode][nBlock_i] = {}

            if nBlock_j not in dictErrorBlockMasks[nROI][nBarcode][nBlock_i].keys():
                dictErrorBlockMasks[nROI][nBarcode][nBlock_i][nBlock_j] = {}

            dictErrorBlockMasks[nROI][nBarcode][nBlock_i][nBlock_j] = {"shift_z":row["shift_z"],
                                                                       "shift_x":row["shift_x"],
                                                                       "shift_y":row["shift_y"],
                                                                       "quality_xy":row["quality_xy"],
                                                                       "quality_zy":row["quality_zy"],
                                                                       "quality_zx":row["quality_zx"]}


        self.dictErrorBlockMasks = dictErrorBlockMasks
        return True


    def register_barcodeMap_file(self, file):

        if "3D" in file:
            self.ndims = 3
        else:
            self.ndims = 2

        # loads barcode coordinate Tables
        table = LocalizationTable()
        barcodeMapFull, uniqueBarcodes = table.load(file)

        # checks that everything is OK
        if len(barcodeMapFull) < 1:
            printLog(f"\nWARNING>{file} contains an empty table!")
            return None

        if "comments" in barcodeMapFull.meta.keys():
            if "registered" in barcodeMapFull.meta["comments"]:
                printLog(f"\nWARNING>{file} contains a table thas was already registered! \nWill not do anything")
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
            nROI = barcodeMap["ROI #"][0]
            printLog(f"\nProcessing barcode localization Table for ROI: {nROI}")

            # registers localizations
            barcodeMap = self.register_barcodes(barcodeMap)

        # saves and plots registered barcode coordinate Tables
        table.save(file, barcodeMap, comments="registered")
        table.plots_distributionFluxes(barcodeMap, [file.split(".")[0], "_registered", "_barcode_stats", ".png"])
        table.plots_localizations(barcodeMap, [file.split(".")[0], "_registered", "_barcode_localizations", ".png"])
        table.compares_localizations(
            barcodeMap, barcodeMapFull_unregistered, [file.split(".")[0], "_registered", "_barcode_comparisons", ".png"]
        )

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
        label = "barcode"

        currentFolder = self.param.param["rootFolder"]
        self.dataFolder.createsFolders(currentFolder, self.param)
        printLog("> Processing Folder: {}".format(currentFolder))

        # Loads localAlignment if it exists
        self.loadsLocalAlignment()

        if not self.alignmentResultsTableRead:
            printLog(f"Unable to find aligment table.\nDid you run alignImages3D?\n\n Aborted.")
            return

        # iterates over barcode localization tables in the current folder
        files = [x for x in glob.glob(self.dataFolder.outputFiles["segmentedObjects"] + "_*" + label + ".dat")]

        if len(files) < 1:
            printLog("No localization table found to process!")
            return

        for file in files:
            self.register_barcodeMap_file(file)

        printLog(f" {len(files)} barcode tables processed in {currentFolder}")
