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

# to remove in a future version
import warnings

import numpy as np
from tqdm import trange

from fileProcessing.fileManagement import Folders, print_log, write_string_to_file
from imageProcessing.localization_table import LocalizationTable

warnings.filterwarnings("ignore")


class FilterLocalizations:
    def __init__(self, param):
        """
        Parameters
        ----------
        param : class
            Parameters
        current_session : class
            session information
        """

        self.current_param = param

    def filter_localizations__quality(self, barcode_map, i):
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
        if self.ndims == 3:  # and  "3DfitKeep" in barcode_map.keys()
            # [reading the flag in barcode_map_roi assigned by the 3D localization routine]
            keep = (
                barcode_map["flux"][i] > self.flux_min
            )  # and barcode_map["3DfitKeep"][i]
        else:
            # [or by reading the flux from 2D localization]
            keep = barcode_map["flux"][i] > self.flux_min

        return keep

    # def filter_localizations_block_alignment(self, barcode_map, i):
    #     """
    #     [filters barcode per blockAlignmentMask, if existing]
    #     runs only if localAligment was not run!

    #     Parameters
    #     ----------
    #     i : int
    #         index in barcode_map Table


    #     Returns
    #     -------
    #     keep_alignment : Boolean
    #         True if the test is passed.

    #     """
    #     y_int = int(barcode_map["xcentroid"][i])
    #     x_int = int(barcode_map["ycentroid"][i])
    #     keep_alignment = True
    #     if (
    #         not self.alignment_results_table_read
    #     ):  # only proceeds if localAlignment was not performed
    #         barcode_id = "barcode:" + str(barcode_map["Barcode #"][i])
    #         barcode_roi = "ROI:" + str(barcode_map["ROI #"][i])

    #         if len(self.dict_error_block_masks) > 0:
    #             if barcode_roi in self.dict_error_block_masks.keys():
    #                 if barcode_id in self.dict_error_block_masks[barcode_roi].keys():
    #                     error_mask = self.dict_error_block_masks[barcode_roi][barcode_id]
    #                     keep_alignment = (
    #                         error_mask[
    #                             int(np.floor(x_int / self.block_size)),
    #                             int(np.floor(y_int / self.block_size)),
    #                         ]
    #                         < self.tolerance_drift
    #                     )

    #         # keeps it always if barcode is fiducial
    #         if (
    #             "RT" + str(barcode_map["Barcode #"][i])
    #             in self.current_param.param_dict["alignImages"]["referenceFiducial"]
    #         ):
    #             keep_alignment = True

    #     return keep_alignment

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
        print(f"$ Minimum flux: {self.flux_min}")
        for i in trange(n_barcodes):  # i is the index of the barcode in barcode_map_roi

            # [filters barcode localizations either by]
            keep_quality = self.filter_localizations__quality(barcode_map, i)

            # [filters barcode per blockAlignmentMask, if existing]
            # keep_alignment = self.filter_localizations_block_alignment(barcode_map, i)
            keep_alignment = True

            if not keep_quality or not keep_alignment:
                rows_to_remove.append(i)

        # removes rows from table
        barcode_map.remove_rows(rows_to_remove)

        print(
            f"$ Removed {len(rows_to_remove)} barcode localizations from table out of {n_barcodes}."
        )

        return barcode_map

    def setup_filter_values(self):
        """


        Returns
        -------
        self.tolerance_drift : float
            tolerance to keep barcode localization, in pixel units
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
            print_log("# Flux min not found. Set to {}!".format(self.flux_min))

        if "toleranceDrift" in self.current_param.param_dict["buildsPWDmatrix"]:
            self.tolerance_drift = self.current_param.param_dict["buildsPWDmatrix"][
                "toleranceDrift"
            ]
        else:
            self.tolerance_drift = 1
            print_log(
                "# toleranceDrift not found. Set to {}!".format(self.tolerance_drift)
            )

        if "blockSize" in self.current_param.param_dict["alignImages"]:
            self.block_size = self.current_param.param_dict["alignImages"]["blockSize"]
        else:
            self.block_size = 256
            print_log("# blockSize not found. Set to {}!".format(self.block_size))

    def filter_folder(self):
        """
        Function that filters barcodes using a number of user-provided parameters


        Returns
        -------
        None.

        """
        session_name = "filter_localizations"

        # processes folders and files
        self.data_folder = Folders(self.current_param.param_dict["rootFolder"])
        print_log("\n===================={}====================\n".format(session_name))
        print_log("$ folders read: {}".format(len(self.data_folder.list_folders)))
        write_string_to_file(
            self.current_param.param_dict["fileNameMD"],
            "## {}\n".format(session_name),
            "a",
        )
        label = "barcode"

        for current_folder in self.data_folder.list_folders:
            self.data_folder.create_folders(current_folder, self.current_param)
            print_log("> Processing Folder: {}".format(current_folder))

            files = [
                x
                for x in glob.glob(
                    self.data_folder.output_files["segmentedObjects"]
                    + "_*"
                    + label
                    + ".dat"
                )
            ]

            if len(files) > 0:

                for file in files:

                    if "3D" in file:
                        self.ndims = 3
                    else:
                        self.ndims = 2

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
                        print("\n$ rois detected: {}".format(number_rois))

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
                        print(f"\nWARNING>{file} contains an empty table!")

            else:
                print_log("No barcode tables found!")

            print_log("Barcode tables {} filtered".format(current_folder), "info")


def get_file_table_new_name(file):

    existing_versions = glob.glob(file.split(".")[0] + "_version_*.dat")

    if len(existing_versions) < 1:
        new_version = 0
    else:
        version_numbers = [
            int(x.split("_version_")[1].split("_")[0]) for x in existing_versions
        ]

        if len(version_numbers) > 0:
            new_version = max(version_numbers) + 1
        else:
            new_version = 0

    new_file = file.split(".dat")[0] + "_version_" + str(new_version) + "_.dat"

    return new_file
