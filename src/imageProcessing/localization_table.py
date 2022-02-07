#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 17:06:14 2022

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

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")


class localization_table:
    def __init__(self):
        self.a = 1

    def load(self, fileNameBarcodeCoordinates):
        """
        Loads barcodeMap

        Parameters
        ----------
        fileNameBarcodeCoordinates : string
            filename with barcodeMap

        Returns
        -------
        barcodeMap : Table()
        uniqueBarcodes: list
            lis of unique barcodes read from barcodeMap

        """
        if os.path.exists(fileNameBarcodeCoordinates):
            barcodeMap = Table.read(fileNameBarcodeCoordinates, format="ascii.ecsv")
            printLog("$ Successfully loaded barcode localizations file: {}".format(fileNameBarcodeCoordinates))

            uniqueBarcodes = np.unique(barcodeMap["Barcode #"].data)
            numberUniqueBarcodes = uniqueBarcodes.shape[0]

            print("Number Barcodes read from barcodeMap: {}".format(numberUniqueBarcodes))
            print("Unique Barcodes detected: {}".format(uniqueBarcodes))
        else:
            print("\n\n# ERROR: could not find coordinates file: {}".format(fileNameBarcodeCoordinates))
            sys.exit()

        return barcodeMap, uniqueBarcodes

    def save(self, fileNameBarcodeCoordinates, barcodeMap, tag="_", ext = 'ecsv'):
        
        file = fileNameBarcodeCoordinates.split('.'+ext)[0] + tag + "." + ext
        
        print(f"Saving {tag} results at {file}")
        
        barcodeMap.write(
            file,
            format="ascii.ecsv",
            overwrite=True,
        )