#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 16:45:44 2022

@author: marcnol
"""
# =============================================================================
# IMPORTS
# =============================================================================

import glob, os, sys
import re
import numpy as np
from tqdm.contrib import tzip
from tqdm import trange
import matplotlib.pyplot as plt

from astropy.table import Table

from photutils.segmentation import SegmentationImage

from fileProcessing.fileManagement import (
    folders,
    writeString2File,
    printLog,
)

from imageProcessing.localization_table import localization_table

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")

class filter_localizations:
    def __init__(self, param):
        """
        Parameters
        ----------
        param : class
            Parameters
        session1 : class
            session information
        """

        self.param = param
        
    def filterLocalizations_Quality(self, barcodeMap, i):
        """
        [filters barcode localizations either by brigthness or 3D localization accuracy]

        Parameters
        ----------
        i : int
            index in barcodeMap Table


        Returns
        -------
        keep : Boolean
            True if the test is passed.

        """
        if "3DfitKeep" in barcodeMap.keys() and self.ndims == 3:
            # [reading the flag in barcodeMapROI assigned by the 3D localization routine]
            keep = barcodeMap["3DfitKeep"][i] and barcodeMap["flux"][i] > self.flux_min
        else:
            # [or by reading the flux from 2D localization]
            keep = barcodeMap["flux"][i] > self.flux_min

        return keep

    def filterLocalizations_BlockAlignment(self, barcodeMap,i):
        """
        [filters barcode per blockAlignmentMask, if existing]
        runs only if localAligment was not run!

        Parameters
        ----------
        i : int
            index in barcodeMap Table


        Returns
        -------
        keepAlignment : Boolean
            True if the test is passed.

        """
        y_int = int(barcodeMap["xcentroid"][i])
        x_int = int(barcodeMap["ycentroid"][i])
        keepAlignment = True
        if not self.alignmentResultsTableRead: # only proceeds if localAlignment was not performed
            barcodeID = "barcode:" + str(barcodeMap["Barcode #"][i])
            barcodeROI = "ROI:" + str(barcodeMap["ROI #"][i])

            if len(self.dictErrorBlockMasks) > 0:
                if barcodeROI in self.dictErrorBlockMasks.keys():
                    if barcodeID in self.dictErrorBlockMasks[barcodeROI].keys():
                        errorMask = self.dictErrorBlockMasks[barcodeROI][barcodeID]
                        keepAlignment = (
                            errorMask[int(np.floor(x_int / self.blockSize)), int(np.floor(y_int / self.blockSize))]
                            < self.toleranceDrift
                        )

            # keeps it always if barcode is fiducial
            if (
                "RT" + str(barcodeMap["Barcode #"][i])
                in self.param.param["alignImages"]["referenceFiducial"]
            ):
                keepAlignment = True

        return keepAlignment

    def filter_barcode_table(self,barcodeMap):
        """
        iterates over rows of a barcode localization table and filters unwanted rows

        Parameters
        ----------
        barcodeMapROI : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        rows_to_remove=list()
        for i in trange(len(barcodeMap)): # i is the index of the barcode in barcodeMapROI


            # [filters barcode localizations either by]
            keepQuality = self.filterLocalizations_Quality(barcodeMap, i)

            # [filters barcode per blockAlignmentMask, if existing]
            #keepAlignment = self.filterLocalizations_BlockAlignment(barcodeMap, i)
            keepAlignment =True
            
            if keepQuality and keepAlignment:
                rows_to_remove.append(i)

        # removes rows from table
        barcodeMap.remove_rows(rows_to_remove)
        
        print(f"$ Removed {len(rows_to_remove)} barcode localizations from table.")
                
        return barcodeMap
    
    def setup_filter_values(self):
        """
        

        Returns
        -------
        self.toleranceDrift : float
            tolerance to keep barcode localization, in pixel units
        self.blockSize : int
            size of blocks used for blockAlignment.
        self.flux_min : float
            Minimum flux to keep barcode localization

        """
        flux_key = "flux_min_3D" if self.ndims==3 else "flux_min"
        if flux_key in self.param.param["buildsPWDmatrix"]:
            self.flux_min = self.param.param["buildsPWDmatrix"][flux_key]
        else:
            self.flux_min = 0
            printLog("# Flux min not found. Set to {}!".format(self.flux_min))

        if "toleranceDrift" in self.param.param["buildsPWDmatrix"]:
            self.tolerance_drift = self.param.param["buildsPWDmatrix"]["toleranceDrift"]
        else:
            self.tolerance_drift= 1
            printLog("# toleranceDrift not found. Set to {}!".format(self.toleranceDrift))
    
        if "blockSize" in self.param.param["alignImages"]:
            self.blockSize = self.param.param["alignImages"]["blockSize"]
        else:
            self.blockSize = 256
            printLog("# blockSize not found. Set to {}!".format(self.blockSize))

                
    def filter_folder(self):
        """
        Function that filters barcodes using a number of user-provided parameters
    
    
        Returns
        -------
        None.
    
        """
        sessionName = "filter_localizations"
    
        # processes folders and files
        self.dataFolder = folders(self.param.param["rootFolder"])
        printLog("\n===================={}====================\n".format(sessionName))
        printLog("$ folders read: {}".format(len(self.dataFolder.listFolders)))
        writeString2File(self.param.param["fileNameMD"], "## {}\n".format(sessionName), "a")
        label = 'barcode'


        for currentFolder in self.dataFolder.listFolders:
            self.dataFolder.createsFolders(currentFolder, self.param)
            printLog("> Processing Folder: {}".format(currentFolder))
    
            files = [x for x in glob.glob(self.dataFolder.outputFiles["segmentedObjects"] + "_*" + label + ".dat")]
            
                
            if len(files)>0:
                
                for file in files:

                    if "3D" in file:
                        self.ndims = 3
                    else:
                        self.ndims = 2
                    
                    self.setup_filter_values()

                    # Loads barcode coordinate Tables
                    table = localization_table() 
                    barcodeMap, uniqueBarcodes = table.load(file)

                    # processes tables
                    barcodeMapROI = barcodeMap.group_by("ROI #")
                    numberROIs = len(barcodeMapROI.groups.keys)
                    print("\n$ ROIs detected: {}".format(numberROIs))
                    
                    # NEED TO ITERATE OVER ROIs !!
                    
                    # Filters barcode coordinate Tables
                    barcodeMap = self.filter_barcode_table(barcodeMap)
                    
                    # saves filtered barcode coordinate Tables
                    table.save(file, barcodeMapROI, tag = '_filtered',ext = 'dat') 
                    
            else:
                printLog("No barcode tables found!")
                
   
            printLog("Barcode tables {} filtered".format(currentFolder), "info")
            

