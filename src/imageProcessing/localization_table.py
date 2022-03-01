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

from apifish.stack.io import save_table_to_ecsv
from apifish.stack.io import read_table_from_ecsv

from fileProcessing.fileManagement import (
    printLog,
)

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")


class LocalizationTable:

    def __init__(self):
        self.a = 1

    def load(self, file):
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
        if os.path.exists(file):
            # barcodeMap = Table.read(fileNameBarcodeCoordinates, format="ascii.ecsv")
            barcodeMap = read_table_from_ecsv(file)

            printLog("$ Successfully loaded barcode localizations file: {}".format(file))

            uniqueBarcodes = np.unique(barcodeMap["Barcode #"].data)
            numberUniqueBarcodes = uniqueBarcodes.shape[0]

            print("Number of barcodes read from barcodeMap: {}".format(numberUniqueBarcodes))
            print("Unique Barcodes detected: {}".format(uniqueBarcodes))
        else:
            print("\n\n# ERROR: could not find coordinates file: {}".format(file))
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

        save_table_to_ecsv(barcodeMap,fileName)

        '''
        barcodeMap.write(
            fileName,
            format="ascii.ecsv",
            overwrite=True,
        )
        '''

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


    def build_color_dict(self,barcodeMap, key='Barcode #'):

        color_dict = dict()

        unique_barcodes = np.unique(barcodeMap[key])
        output_array = range(unique_barcodes.shape[0])

        for barcode, output in zip(unique_barcodes,output_array):
            color_dict[str(barcode)]=output



        return color_dict


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
            color_dict = self.build_color_dict(barcodeMap, key='Barcode #')
            Nbarcodes = np.unique(barcodeMap['Barcode #']).shape[0]

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
            colors =  [color_dict[str(x)] for x in barcodeMap["Barcode #"]]
            titles = ["Z-projection", "X-projection", "Y-projection"]

            # makes plot
            plots_localization_projection(x,y,ax[0], colors, titles[0])
            plots_localization_projection(x,z,ax[1], colors, titles[1])
            plots_localization_projection(y,z,ax[2], colors, titles[2])

            fig.tight_layout()

            # saves output figure
            fileName_list_i=fileName_list.copy()
            fileName_list_i.insert(-1,'_ROI' + str(nROI))
            fig.savefig("".join(fileName_list_i))

    def decode_ROIs(self, barcodeMap):

        return decode_ROIs(barcodeMap)

        '''
        barcodeMapROI = barcodeMap.group_by("ROI #")

        numberROIs = len(barcodeMapROI.groups.keys)

        print("\n$ ROIs detected: {}".format(numberROIs))

        return barcodeMapROI,numberROIs
        '''

    def compares_localizations(self,barcodeMap1,barcodeMap2,fileName_list, fontsize=20):
        """
        Compares the localizations of two barcode tables

        Parameters
        ----------
        barcodeMap1 : astropy Table
            localization table 1.
        barcodeMap2 : astropy Table
            localization table 2.

        Returns
        -------
        None.

        """

        barcodeMap2.add_index('Buid')
        number_localizations = len(barcodeMap1)

        diffs = dict()
        labels = ['xcentroid','ycentroid','zcentroid']
        for label in labels:
            diffs[label]=[]

        # iterates over rows in barcodeMap1
        for row in range(number_localizations):
            Buid_1 = barcodeMap2[row]['Buid']
            barcode_found = True

            # finds same Buid in barcodeMap2
            try:
                barcodeMap2.loc[Buid_1]
            except KeyError:
                barcode_found = False
                pass

            # collects differences in values between same localization in both tables
            if barcode_found:
                for label in labels:
                    diff = barcodeMap2.loc[Buid_1][label]-barcodeMap1[row][label]
                    if np.isnan(diff):
                        diff = 0
                    diffs[label].append(diff)

        # plots figures
        fig, axes = plt.subplots(2, 2)
        ax = axes.ravel()
        fig.set_size_inches((30, 30))

        for label, axis in zip(labels,ax):
            r = np.array(diffs[label])
            axis.hist(r, bins=20)
            axis.set_xlabel(label+" correction, px", fontsize=fontsize)
            axis.set_ylabel("counts", fontsize=fontsize)

        ax[3].scatter(np.array(diffs['ycentroid']),np.array(diffs['xcentroid']),s=3, alpha=.8)
        ax[3].set_xlabel("dx-position, px", fontsize=fontsize)
        ax[3].set_ylabel("dy-position, px", fontsize=fontsize)

        fig.savefig("".join(fileName_list))



def decode_ROIs(data):

    data_indexed = data.group_by("ROI #")

    numberROIs = len(data_indexed.groups.keys)

    print("\n$ ROIs detected: {}".format(numberROIs))

    return data_indexed, numberROIs


def build_color_dict(data, key='Barcode #'):

    color_dict = dict()

    unique_barcodes = np.unique(data[key])
    output_array = range(unique_barcodes.shape[0])

    for barcode, output in zip(unique_barcodes,output_array):
        color_dict[str(barcode)]=output

    return color_dict

def plots_localization_projection(coord1, coord2, axis, colors, title=''*3):
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