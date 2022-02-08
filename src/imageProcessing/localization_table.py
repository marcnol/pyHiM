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

    def save(self, fileName, barcodeMap, tag="_", ext = 'ecsv', comments=list()):
        """
        Saves output table

        Parameters
        ----------
        fileNameBarcodeCoordinates : string
            filename of table.
        barcodeMap : astropy Table
            Table to be written to file.
        tag : string, optional
            tag to be added to filename. The default is "_".
        ext : string, optional
            file extension. The default is 'ecsv'.
        comments : list of strings, optional
            Will output as comments to the header. The default is [].

        Returns
        -------
        None.

        """
        file = fileName.split('.'+ext)[0] + tag + "." + ext

        print(f"Saving {tag} results at {file}")

        barcodeMap.meta['comments']=comments

        barcodeMap.write(
            file,
            format="ascii.ecsv",
            overwrite=True,
        )

    def plots_distributionFluxes(self, barcodeMap, fileName):


        fig, axes = plt.subplots(1, 2)
        ax = axes.ravel()
        fig.set_size_inches((10, 5))

        fluxes = barcodeMap["flux"]
        sharpness = barcodeMap["sharpness"]
        roundness = barcodeMap["roundness1"]
        peak = barcodeMap["peak"]
        mag = barcodeMap["mag"]

        p1 = ax[0].scatter(fluxes, sharpness, c=peak, cmap="terrain", alpha=0.5)
        ax[0].set_title("color: peak intensity")
        ax[0].set_xlabel("flux")
        ax[0].set_ylabel("sharpness")

        p2 = ax[1].scatter(roundness, mag, c=peak, cmap="terrain", alpha=0.5)
        ax[1].set_title("color: peak intensity")
        ax[1].set_xlabel("roundness")
        ax[1].set_ylabel("magnitude")
        fig.colorbar(p2, ax=ax[1], fraction=0.046, pad=0.04)

        fig.savefig(fileName)

        plt.close(fig)
