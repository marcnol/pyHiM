#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 16:45:44 2022

@author: marcnol
"""
# =============================================================================
# IMPORTS
# =============================================================================

import glob
import os

from tqdm import trange

from core.parameters import MatrixParams
from core.pyhim_logging import print_log, print_session_name, write_string_to_file
from imageProcessing.localization_table import LocalizationTable
from imageProcessing.makeProjections import Feature


class FilterLocalizationsTempo(Feature):
    def __init__(self, params: MatrixParams):
        super().__init__(params)
        self.out_folder = self.params.folder
        self.name = "FilterLocalizations"


class FilterLocalizations:
    def __init__(self, param):
        """
        Parameters
        ----------
        param : class
            Parameters
        """

        self.current_param = param

    def filter_localizations_quality(self, barcode_map, i):
        """
        [filters barcode localizations either by brigthness or 3D localization accuracy]

        Parameters
        ----------
        i : int
            index in barcode_map Table


        Returns
        -------
        keep : Boolean
            True if the test is passed.

        """
        return barcode_map["flux"][i] > self.flux_min

    def filter_barcode_table(self, barcode_map):
        """
        iterates over rows of a barcode localization table and filters unwanted rows

        Parameters
        ----------
        barcode_map_roi : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        rows_to_remove = []
        n_barcodes = len(barcode_map)
        print_log(f"$ Minimum flux: {self.flux_min}")
        # [filters barcode per blockAlignmentMask, if existing]
        # keep_alignment = self.filter_localizations_block_alignment(barcode_map, i)
        keep_alignment = True

        for i in trange(n_barcodes):  # i is the index of the barcode in barcode_map_roi
            # [filters barcode localizations either by]
            keep_quality = self.filter_localizations_quality(barcode_map, i)

            if not keep_quality or not keep_alignment:
                rows_to_remove.append(i)

        # removes rows from table
        barcode_map.remove_rows(rows_to_remove)

        print_log(
            f"$ Removed {len(rows_to_remove)} barcode localizations from table out of {n_barcodes}."
        )

        return barcode_map

    def setup_filter_values(self):
        """
        Returns
        -------
        self.block_size : int
            size of blocks used for blockAlignment.
        self.flux_min : float
            Minimum flux to keep barcode localization

        """
        flux_key = "flux_min_3D" if self.ndims == 3 else "flux_min"
        if flux_key in self.current_param.param_dict["buildsPWDmatrix"]:
            self.flux_min = self.current_param.param_dict["buildsPWDmatrix"][flux_key]
        else:
            self.flux_min = 0
            print_log(f"# Flux min not found. Set to {self.flux_min}!")

        if "blockSize" in self.current_param.param_dict["alignImages"]:
            self.block_size = self.current_param.param_dict["alignImages"]["blockSize"]
        else:
            self.block_size = 256
            print_log(f"# blockSize not found. Set to {self.block_size}!")

    def filter_folder(self, data_path, seg_params):
        """
        Function that filters barcodes using a number of user-provided parameters


        Returns
        -------
        None.

        """
        session_name = "filter_localizations"

        # processes folders and files
        print_session_name(session_name)
        write_string_to_file(
            self.current_param.param_dict["fileNameMD"],
            f"## {session_name}\n",
            "a",
        )

        current_folder = self.current_param.param_dict["rootFolder"]
        print_log(f"> Processing Folder: {current_folder}")

        data_file_base = (
            data_path
            + os.sep
            + seg_params.folder
            + os.sep
            + "data"
            + os.sep
            + seg_params.outputFile
        )
        files = list(glob.glob(data_file_base + "_*barcode.dat"))
        if files:
            for file in files:
                self.ndims = 3 if "3D" in os.path.basename(file) else 2
                self.setup_filter_values()

                # Loads barcode coordinate Tables
                table = LocalizationTable()
                barcode_map, unique_barcodes = table.load(file)

                if len(barcode_map) > 0:
                    # plots and saves original barcode coordinate Tables for safe keeping
                    new_file = get_file_table_new_name(file)
                    table.save(new_file, barcode_map)
                    table.plot_distribution_fluxes(
                        barcode_map,
                        [new_file.split(".")[0], "_barcode_stats", ".png"],
                    )
                    table.plots_localizations(
                        barcode_map,
                        [new_file.split(".")[0], "_barcode_localizations", ".png"],
                    )

                    # processes tables
                    barcode_map_roi = barcode_map.group_by("ROI #")
                    number_rois = len(barcode_map_roi.groups.keys)
                    print_log(f"\n$ rois detected: {number_rois}")

                    # Filters barcode coordinate Tables
                    barcode_map = self.filter_barcode_table(barcode_map)

                    # saves and plots filtered barcode coordinate Tables
                    table.save(file, barcode_map, comments="filtered")
                    table.plot_distribution_fluxes(
                        barcode_map, [file.split(".")[0], "_barcode_stats", ".png"]
                    )
                    table.plots_localizations(
                        barcode_map,
                        [file.split(".")[0], "_barcode_localizations", ".png"],
                    )

                else:
                    print_log(f"\nWARNING>{file} contains an empty table!")

        else:
            print_log("No barcode tables found!")

        print_log(f"Barcode tables {current_folder} filtered", "info")


def get_file_table_new_name(file):
    existing_versions = glob.glob(file.split(".")[0] + "_version_*.dat")

    if len(existing_versions) < 1:
        new_version = 0
    else:
        version_numbers = [
            int(x.split("_version_")[1].split("_")[0]) for x in existing_versions
        ]

        new_version = max(version_numbers) + 1 if version_numbers else 0
    return file.split(".dat")[0] + "_version_" + str(new_version) + "_.dat"
