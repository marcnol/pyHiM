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

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from apifish.stack.io import read_table_from_ecsv, save_table_to_ecsv
from astropy.table import Table, vstack
from stardist import random_label_cmap

from fileProcessing.fileManagement import print_log
from imageProcessing.localization_table import (
    build_color_dict,
    decode_rois,
    plots_localization_projection,
)

lbl_cmap = random_label_cmap()

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")


class ChromatinTraceTable:
    def __init__(self, xyz_unit="micron", genome_assembly="mm10"):
        self.a = 1
        self.xyz_unit = xyz_unit
        self.genome_assembly = genome_assembly

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
                "label",
            ),
            dtype=("S2", "S2", "f4", "f4", "f4", "S2", "int", "int", "int", "int", "int", "S2",),
        )

        self.data.meta["comments"] = [
            "xyz_unit={}".format(self.xyz_unit),
            "genome_assembly={}".format(self.genome_assembly),
        ]

    def load(self, file):
        """
        Loads chromatin trace table

        Parameters
        ----------
        filename_barcode_coordinates : string
            filename with chromatin trace table

        Returns
        -------
        chromatin trace table : Table()
        unique_barcodes: list
            lis of unique barcodes read from chromatin trace table

        """
        if os.path.exists(file):
            trace_table = read_table_from_ecsv(file)
            print_log("$ Successfully loaded chromatin trace table: {}".format(file))
        else:
            print("\n\n# ERROR: could not find chromatin trace table: {}".format(file))
            sys.exit()

        self.data = trace_table

        return trace_table

    def save(self, file_name, table, comments=""):
        """
        Saves output table

        Parameters
        ----------
        file_name: string
            filename of table.
        table: astropy Table
            Table to be written to file.
        comments : list of strings, optional
            Will output as comments to the header. The default is [].

        Returns
        -------
        None.

        """
        print(f"Saving output table as {file_name} ...")

        try:
            table.meta["comments"].append(comments)
        except KeyError:
            table.meta["comments"] = [comments]

        save_table_to_ecsv(table, file_name)

        """
        table.write(
            file_name,
            format="ascii.ecsv",
            overwrite=True,
        )
        """

    def append(self, table):
        """
        appends <table> to self.data

        Parameters
        ----------
        table : astropy table
            table to append to existing self.data table.

        Returns
        -------
        None.

        """

        self.data = vstack([self.data, table])

    def filter_traces_by_n(self, minimum_number_barcodes=2):
        """
        Removes rows in trace table with less than `minimum_number_barcodes` barcodes

        Parameters
        ----------
        trace_table : ASTROPY Table
            input trace table.
        minimum_number_barcodes : TYPE, optional
            minimum number of barcodes in trace. The default is 1.

        Returns
        -------
        trace_table : ASTROPY Table
            output trace table.

        """

        trace_table = self.data

        # indexes trace file
        trace_table_indexed = trace_table.group_by("Trace_ID")

        # iterates over traces
        print(f"\n$ WIll keep traces with {minimum_number_barcodes } spots")
        print(
            f"$ Number of original spots / traces: {len(trace_table)} / {len(trace_table_indexed.groups)}"
        )

        barcodes_to_remove = []

        for idx, trace in enumerate(trace_table_indexed.groups):

            number_unique_barcodes = len(list(set(trace["Barcode #"].data)))
            
            if number_unique_barcodes < minimum_number_barcodes:
                barcodes_to_remove.append(list(trace["Spot_ID"].data))

        print(f"$ Number of traces to remove: {len(barcodes_to_remove)}")

        list_barcode_to_remove = []
        for barcodes in barcodes_to_remove:
            for x in barcodes:
                list_barcode_to_remove.append(x)

        rows_to_remove = []

        for idx, row in enumerate(trace_table):
            spot_id = row["Spot_ID"]

            if spot_id in list_barcode_to_remove:
                rows_to_remove.append(idx)

        trace_table.remove_rows(rows_to_remove)
        if len(trace_table) > 0:
            trace_table_indexed = trace_table.group_by("Trace_ID")
            number_traces_left = len(trace_table_indexed.groups)
        else:
            number_traces_left = 0
            
        print(f"$ Number of spots / traces left: {len(trace_table)} / {number_traces_left}")

        self.data = trace_table

    def plots_traces(
        self, filename_list, masks=np.zeros((2048, 2048)), pixel_size=[0.1, 0.1, 0.25]
    ):

        """
        This function plots 3 subplots (xy, xz, yz) with the localizations.
        One figure is produced per ROI.

        Parameters
        ----------

        filename_list: list
            filename
        """

        data = self.data

        # indexes table by ROI
        data_indexed, number_rois = decode_rois(data)

        for i_roi in range(number_rois):

            # creates sub Table for this ROI
            data_roi = data_indexed.groups[i_roi]
            n_roi = data_roi["ROI #"][0]
            print(f"Plotting barcode localization map for ROI: {n_roi}")
            color_dict = build_color_dict(data_roi, key="Barcode #")
            n_barcodes = np.unique(data_roi["Barcode #"]).shape[0]

            # initializes figure
            fig = plt.figure(constrained_layout=False)
            im_size = 60
            fig.set_size_inches((im_size * 2, im_size))
            gs = fig.add_gridspec(2, 2)
            ax = [
                fig.add_subplot(gs[:, 0]),
                fig.add_subplot(gs[0, 1]),
                fig.add_subplot(gs[1, 1]),
            ]

            # defines variables
            x = data_roi["x"] / pixel_size[0]
            y = data_roi["y"] / pixel_size[1]
            z = data_roi["z"] / pixel_size[2]

            colors = [color_dict[str(x)] for x in data_roi["Barcode #"]]
            titles = ["Z-projection", "X-projection", "Y-projection"]

            # plots masks if available
            ax[0].imshow(masks, cmap=lbl_cmap, alpha=0.3)

            # makes plot
            plots_localization_projection(x, y, ax[0], colors, titles[0])
            plots_localization_projection(x, z, ax[1], colors, titles[1])
            plots_localization_projection(y, z, ax[2], colors, titles[2])

            fig.tight_layout()

            # calculates mean trace positions and sizes by looping over traces
            data_traces = data_roi.group_by("Trace_ID")
            number_traces = len(data_traces.groups.keys)
            color_dict_traces = build_color_dict(data_traces, key="Trace_ID")
            colors_traces = [color_dict_traces[str(x)] for x in data_traces["Trace_ID"]]
            for trace, color, trace_id in zip(
                data_traces.groups, colors_traces, data_traces.groups.keys
            ):
                x_trace = np.mean(trace["x"].data) / pixel_size[0]
                y_trace = np.mean(trace["y"].data) / pixel_size[1]
                z_trace = np.mean(trace["z"].data) / pixel_size[2]
                s_trace = (
                    300
                    * (
                        np.mean(
                            [
                                np.std(trace["x"].data),
                                np.std(trace["y"].data),
                                np.std(trace["z"].data),
                            ]
                        )
                    )
                    / pixel_size[0]
                )

                # plots circles for each trace
                ax[0].scatter(
                    x_trace,
                    y_trace,
                    s=s_trace,
                    c=color,
                    marker="$\u25EF$",
                    cmap="nipy_spectral",
                    linewidths=1,
                    alpha=0.7,
                )

            # saves output figure
            filename_list_i = filename_list.copy()
            filename_list_i.insert(-1, "_ROI" + str(n_roi))

            try:
                fig.savefig("".join(filename_list_i))
            except ValueError:
                print(
                    "\nValue error while saving output figure with traces:{}".format(
                        "".join(filename_list_i)
                    )
                )
