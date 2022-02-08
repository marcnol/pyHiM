#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 17:06:14 2022

@author: marcnol

This class will contain methods to load, save, plot barcode localizations and statistics

"""

# =============================================================================
# IMPORTS
# =============================================================================

import os, sys
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table

from fileProcessing.fileManagement import (
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

            print("Number of barcodes read from barcodeMap: {}".format(numberUniqueBarcodes))
            print("Unique Barcodes detected: {}".format(uniqueBarcodes))
        else:
            print("\n\n# ERROR: could not find coordinates file: {}".format(fileNameBarcodeCoordinates))
            sys.exit()

        return barcodeMap, uniqueBarcodes

    def save(self, fileName, barcodeMap, comments=''):
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
        print(f"Saving output table as {fileName} ...")

        try:
            barcodeMap.meta['comments'].append(comments)
        except KeyError:
            barcodeMap.meta['comments']=[comments]

        barcodeMap.write(
            fileName,
            format="ascii.ecsv",
            overwrite=True,
        )

    def plots_distributionFluxes(self, barcodeMap, fileName_list):
        """
        This function plots the distribution of fluxes, sharpness, roundness, magnitude and peak intensity from a Table

        Parameters
        ----------
        barcodeMap : TYPE
            DESCRIPTION.
        fileName_list: list
            filename

        Returns
        -------
        None.

        """

        # initializes figure
        fig, axes = plt.subplots(1, 2)
        ax = axes.ravel()
        fig.set_size_inches((10, 5))

        # initializes variables
        fluxes = barcodeMap["flux"]
        sharpness = barcodeMap["sharpness"]
        roundness = barcodeMap["roundness1"]
        peak = barcodeMap["peak"]
        mag = barcodeMap["mag"]

        # plots data
        p1 = ax[0].scatter(fluxes, sharpness, c=peak, cmap="terrain", alpha=0.5)
        ax[0].set_title("color: peak intensity")
        ax[0].set_xlabel("flux")
        ax[0].set_ylabel("sharpness")

        p2 = ax[1].scatter(roundness, mag, c=peak, cmap="terrain", alpha=0.5)
        ax[1].set_title("color: peak intensity")
        ax[1].set_xlabel("roundness")
        ax[1].set_ylabel("magnitude")
        fig.colorbar(p2, ax=ax[1], fraction=0.046, pad=0.04)

        # saves figure
        fig.savefig("".join(fileName_list))

        plt.close(fig)

    def plots_localization_projection(self, coord1, coord2, axis, colors, title=''*3):
        """
        This function will produce the scatter plot and add title

        Parameters
        ----------
        coord1 : 1D Numpy array, float
            first coordinate (x).
        coord2 : 1D Numpy array, float
            first coordinate (y).
        axis : matplotlib axis
            figure axis handle.
        colors : 1D Numpy array, float
            colorcode used in scatter plot.
        title : string, optional
            title of subpanel. The default is ''*3.

        Returns
        -------
        None.

        """
        axis.scatter(coord1,coord2, s=5, c=colors, alpha=.9, cmap = 'hsv') #nipy_spectral
        axis.set_title(title)

    def plots_localizations(self, barcodeMapFull, fileName_list):

        """
        This function plots 3 subplots (xy, xz, yz) with the localizations.
        One figure is produced per ROI.

        Parameters
        ----------
        image : List of numpy ndarray (N-dimensional array)
            3D raw image of format .tif

        label : List of numpy ndarray (N-dimensional array)
            3D labeled image of format .tif

        fileName_list: list
            filename
        """

        # indexes table by ROI
        barcodeMapROI,numberROIs = self.decode_ROIs(barcodeMapFull)

        for iROI in range(numberROIs):

            # creates sub Table for this ROI
            barcodeMap = barcodeMapROI.groups[iROI]
            nROI = barcodeMap['ROI #'][0]
            print(f"Plotting barcode localization map for ROI: {nROI}")

            # initializes figure
            fig = plt.figure(constrained_layout=False)
            im_size = 60
            fig.set_size_inches((im_size * 2, im_size))
            gs = fig.add_gridspec(2, 2)
            ax = [fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])]

            # defines variables
            x = barcodeMap["xcentroid"]
            y = barcodeMap["ycentroid"]
            z = barcodeMap["zcentroid"]
            colors =  barcodeMap["Barcode #"] # np.random.rand(len(x))
            titles = ["Z-projection", "X-projection", "Y-projection"]

            # makes plot
            self.plots_localization_projection(x,y,ax[0], colors, titles[0])
            self.plots_localization_projection(x,z,ax[1], colors, titles[1])
            self.plots_localization_projection(y,z,ax[2], colors, titles[2])

            fig.tight_layout()

            # saves output figure
            fileName_list_i=fileName_list.copy()
            fileName_list_i.insert(-1,'_ROI' + str(nROI))
            fig.savefig("".join(fileName_list_i))

    def decode_ROIs(self, barcodeMap):

        barcodeMapROI = barcodeMap.group_by("ROI #")

        numberROIs = len(barcodeMapROI.groups.keys)

        print("\n$ ROIs detected: {}".format(numberROIs))

        return barcodeMapROI,numberROIs