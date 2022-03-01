#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:33:57 2022

@author: marcnol

trace table management class

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

from imageProcessing.localization_table import decode_ROIs, build_color_dict, plots_localization_projection

from stardist import random_label_cmap
lbl_cmap = random_label_cmap()

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")


class ChromatinTraceTable:

    def __init__(self, XYZ_unit='micron', genome_assembly='mm10'):
        self.a = 1
        self.XYZ_unit=XYZ_unit
        self.genome_assembly=genome_assembly

    def initialize(self):
        self.data = Table(
                names=(
                    "Spot_ID",
                    "Trace_ID",
                    "x",
                    "y",
                    "z",
                    "Chrom",
                    "Chrom_Start",
                    "Chrom_End",
                    "ROI #",
                    "Mask_id",
                    "Barcode #",
                ),
                dtype=("S2",
                       "S2",
                       "f4",
                       "f4",
                       "f4",
                       "S2",
                       "int",
                       "int",
                       "int",
                       "int",
                       "int",
                       ),
            )
        self.data.meta['comments']=["XYZ_unit={}".format(self.XYZ_unit),
                                    "genome_assembly={}".format(self.genome_assembly),
                                    ]

    def load(self, file):
        """
        Loads chromatin trace table

        Parameters
        ----------
        fileNameBarcodeCoordinates : string
            filename with chromatin trace table

        Returns
        -------
        chromatin trace table : Table()
        uniqueBarcodes: list
            lis of unique barcodes read from chromatin trace table

        """
        if os.path.exists(file):
            trace_table = Table.read(file, format="ascii.ecsv")
            printLog("$ Successfully loaded chromatin trace table: {}".format(file))
        else:
            print("\n\n# ERROR: could not find chromatin trace table: {}".format(file))
            sys.exit()

        self.data = trace_table

        return trace_table


    def save(self, fileName, table, comments=''):
        """
        Saves output table

        Parameters
        ----------
        fileName: string
            filename of table.
        table: astropy Table
            Table to be written to file.
        comments : list of strings, optional
            Will output as comments to the header. The default is [].

        Returns
        -------
        None.

        """
        print(f"Saving output table as {fileName} ...")

        try:
            table.meta['comments'].append(comments)
        except KeyError:
            table.meta['comments']=[comments]

        table.write(
            fileName,
            format="ascii.ecsv",
            overwrite=True,
        )

    def plots_traces(self, fileName_list, Masks = np.zeros((2048,2048)),pixelSize = [0.1,0.1,0.25] ):

        """
        This function plots 3 subplots (xy, xz, yz) with the localizations.
        One figure is produced per ROI.

        Parameters
        ----------

        fileName_list: list
            filename
        """

        data = self.data

        # indexes table by ROI
        data_indexed, numberROIs = decode_ROIs(data)

        for iROI in range(numberROIs):

            # creates sub Table for this ROI
            data_ROI = data_indexed.groups[iROI]
            nROI = data_ROI['ROI #'][0]
            print(f"Plotting barcode localization map for ROI: {nROI}")
            color_dict = build_color_dict(data_ROI, key='Barcode #')
            Nbarcodes = np.unique(data_ROI['Barcode #']).shape[0]

            # initializes figure
            fig = plt.figure(constrained_layout=False)
            im_size = 60
            fig.set_size_inches((im_size * 2, im_size))
            gs = fig.add_gridspec(2, 2)
            ax = [fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])]

            # defines variables
            x = data_ROI["x"]/pixelSize[0]
            y = data_ROI["y"]/pixelSize[1]
            z = data_ROI["z"]/pixelSize[2]

            colors =  [color_dict[str(x)] for x in data_ROI["Barcode #"]]
            titles = ["Z-projection", "X-projection", "Y-projection"]

            # plots masks if available
            ax[0].imshow(Masks, cmap=lbl_cmap,alpha=.3)

            # makes plot
            plots_localization_projection(x,y,ax[0], colors, titles[0])
            plots_localization_projection(x,z,ax[1], colors, titles[1])
            plots_localization_projection(y,z,ax[2], colors, titles[2])

            fig.tight_layout()

            # calculates mean trace positions and sizes by looping over traces
            data_traces = data_ROI.group_by("Trace_ID")
            number_traces = len(data_traces.groups.keys)
            color_dict_traces = build_color_dict(data_traces, key='Trace_ID')
            colors_traces =  [color_dict_traces[str(x)] for x in data_traces["Trace_ID"]]
            s_traces=[]
            for trace,color,trace_id in zip(data_traces.groups,colors_traces,data_traces.groups.keys):
                x_trace = np.mean(trace['x'].data)/pixelSize[0]
                y_trace = np.mean(trace['y'].data)/pixelSize[1]
                z_trace = np.mean(trace['z'].data)/pixelSize[2]
                s_trace = 300*(np.mean([np.std(trace['x'].data),np.std(trace['y'].data),np.std(trace['z'].data)]))/pixelSize[0]

                # plots circles for each trace
                sc = ax[0].scatter(x_trace,y_trace, s = s_trace, c = color, marker="$\u25EF$", cmap = "nipy_spectral", linewidths=1, alpha = 0.7)

            # saves output figure
            fileName_list_i=fileName_list.copy()
            fileName_list_i.insert(-1,'_ROI' + str(nROI))

            try:
                fig.savefig("".join(fileName_list_i))
            except ValueError:
                print("\nValue error while saving output figure with traces:{}".format("".join(fileName_list_i)))
