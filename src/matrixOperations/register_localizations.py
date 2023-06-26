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

import glob
import os
import sys

# to remove in a future version
import warnings

import numpy as np
from astropy.table import Table
from tqdm import trange

from fileProcessing.fileManagement import Folders, print_log, write_string_to_file
from imageProcessing.localization_table import LocalizationTable
from matrixOperations.filter_localizations import get_file_table_new_name

warnings.filterwarnings("ignore")


class RegisterLocalizations:
    def __init__(self, param):
        """
        Parameters
        ----------
        param : class
            Parameters
        """

        self.current_param = param
        self.alignment_results_table_read = False
        self.found_match = []
        self.local_alignment_filename = None
        self.alignment_results_table = None
        self.dict_error_block_masks = None
        self.ndims = None
        self.data_folder = None

        if "toleranceDrift" in self.current_param.param_dict["buildsPWDmatrix"]:
            self.tolerance_drift = self.current_param.param_dict["buildsPWDmatrix"][
                "toleranceDrift"
            ]
            if isinstance(self.tolerance_drift, int):
                # defines a tuple suitable for anisotropic tolerance_drift (z,x,y)
                self.tolerance_drift = (
                    self.tolerance_drift,
                    self.tolerance_drift,
                    self.tolerance_drift,
                )
            elif len(self.tolerance_drift) != 3:
                self.tolerance_drift = (
                    3,
                    1,
                    1,
                )  # defines default anisotropic tolerance_drift (z,x,y)
            elif isinstance(self.tolerance_drift, list):
                self.tolerance_drift = tuple(self.tolerance_drift)
        else:
            self.tolerance_drift = (
                3,
                1,
                1,
            )  # defines default anisotropic tolerance_drift (z,x,y)
            print_log(
                "# toleranceDrift not found. Set to {}!".format(self.tolerance_drift)
            )

        if (
            "remove_uncorrected_localizations"
            in self.current_param.param_dict["buildsPWDmatrix"]
        ):
            self.remove_uncorrected_localizations = self.current_param.param_dict[
                "buildsPWDmatrix"
            ]["remove_uncorrected_localizations"]
        else:
            self.remove_uncorrected_localizations = True

        if self.remove_uncorrected_localizations:
            print_log("# Uncorrected localizations will be removed!!")
        else:
            print_log("# Uncorrected localizations will be kept!!")

    def search_local_shift(self, roi, barcode, zxy_uncorrected):
        if self.alignment_results_table_read:
            return self.search_local_shift_block_3d(roi, barcode, zxy_uncorrected)
        else:  # no correction was applied because the localAlignmentTable was not found
            print("ERROR> did not found alignment_results_table")
            return zxy_uncorrected, {"below_tolerance": False}

    def search_local_shift_block_3d(self, roi, barcode, zxy_uncorrected):
        """
        Searches for local drift for a specific barcode in a given roi.
        If it exists then it adds to the uncorrected coordinates

        Parameters
        ----------
        roi : int
            roi used
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
        block_size_xy = self.alignment_results_table[0]["blockXY"]

        # zxy coord in block reference coord system
        zxy_block = [np.floor(a / block_size_xy).astype(int) for a in zxy_uncorrected]

        # defines roi, barcode ID, and blocks for the localization
        n_roi = "ROI:" + str(roi)
        n_barcode = "barcode:" + "RT" + str(barcode)
        n_block_i = "block_i:" + str(zxy_block[1])
        n_block_j = "block_j:" + str(zxy_block[2])

        # finds the corresponding shift int the dictionary
        shifts = [
            self.dict_error_block_masks[n_roi][n_barcode][n_block_i][n_block_j][
                "shift_z"
            ],
            self.dict_error_block_masks[n_roi][n_barcode][n_block_i][n_block_j][
                "shift_x"
            ],
            self.dict_error_block_masks[n_roi][n_barcode][n_block_i][n_block_j][
                "shift_y"
            ],
        ]

        accepts_localization = False
        if isinstance(self.tolerance_drift, tuple):
            # makes list with comparisons per axis
            check = [
                np.abs(shift) < tol for shift, tol in zip(shifts, self.tolerance_drift)
            ]
            if all(check):
                accepts_localization = True  # only if tolerance is passed in all axes the localization is kept
        else:  # defaults to previous usage with isotropic tolerance
            if max(np.abs(shifts)) < self.tolerance_drift:
                accepts_localization = True

        if accepts_localization:
            zxy_corrected = [a + shift for a, shift in zip(zxy_uncorrected, shifts)]
            quality_correction = {"below_tolerance": True}
        else:
            zxy_corrected = zxy_uncorrected
            quality_correction = {"below_tolerance": False}

        return zxy_corrected, quality_correction

    def register_barcodes(self, barcode_map):
        """
        This function will take a barcode_map and a Table of 3D alignments to register barcode coordinates

        Returns
        -------
        None.

        """

        reference_fiducial = self.current_param.param_dict["alignImages"][
            "referenceFiducial"
        ]

        if "blockSize" in self.current_param.param_dict["alignImages"]:
            block_size = self.current_param.param_dict["alignImages"]["blockSize"]
        else:
            block_size = 256
            print_log("# blockSize not found. Set to {}!".format(block_size))

        print_log(
            f"\n$ Parameters:\n Blocksize = {block_size}\n Tolerance = {self.tolerance_drift}\n Reference barcode = {reference_fiducial}"
        )

        list_uncorrected_barcodes = []

        # loops over barcode Table rows in a given roi
        for i in trange(
            len(barcode_map.groups[0])
        ):  # i is the index of the barcode in barcode_map_roi
            barcode = barcode_map.groups[0]["Barcode #"][i]
            roi = barcode_map.groups[0]["ROI #"][i]

            # keeps the particle if the test passed
            x_uncorrected = barcode_map.groups[0]["ycentroid"][i]
            y_uncorrected = barcode_map.groups[0]["xcentroid"][i]
            z_uncorrected = barcode_map.groups[0]["zcentroid"][i]

            # Corrects XYZ coordinate of barcode if localDriftCorrection is available
            zxy_uncorrected = [z_uncorrected, x_uncorrected, y_uncorrected]
            rt_barcode = "RT" + str(barcode)

            if (
                rt_barcode
                not in self.current_param.param_dict["alignImages"]["referenceFiducial"]
            ):
                zxy_corrected, quality_correction = self.search_local_shift(
                    roi, barcode, zxy_uncorrected
                )
                if not quality_correction["below_tolerance"]:
                    list_uncorrected_barcodes.append(i)
                    if self.remove_uncorrected_localizations:
                        # will remove localizations that cannot be corrected
                        zxy_corrected = [np.nan, np.nan, np.nan]
                    else:
                        # will keep uncorrected localizations
                        pass

            else:
                # if it is the reference cycle, then it does not correct coordinates
                zxy_corrected = zxy_uncorrected

            # rewrites corrected XYZ values to Table
            # if not quality_correction['below_tolerance'] and self.remove_uncorrected_localizations:
            # print(f" $ Before correction: {barcodeMap.groups[0][i]} ")

            barcode_map.groups[0]["ycentroid"][i] = zxy_corrected[1]
            barcode_map.groups[0]["xcentroid"][i] = zxy_corrected[2]
            if self.ndims > 2:
                barcode_map.groups[0]["zcentroid"][i] = zxy_corrected[0]

            # if not quality_correction['below_tolerance'] and self.remove_uncorrected_localizations:
            # print(f" $ After correction: {barcodeMap.groups[0][i]} ")

        if self.remove_uncorrected_localizations:
            print_log(
                f"$ {len(list_uncorrected_barcodes)} localizations out of {len(barcode_map.groups[0])} were removed."
            )
        else:
            print_log(
                f"$ {len(list_uncorrected_barcodes)} localizations out of {len(barcode_map.groups[0])} were uncorrected."
            )

        return barcode_map

    def load_local_alignment(self):
        if "None" in self.current_param.param_dict["alignImages"]["localAlignment"]:
            print_log(
                "\n\n$ localAlignment option set to {}".format(
                    self.current_param.param_dict["alignImages"]["localAlignment"]
                )
            )
            return False, Table()
        else:
            return self._load_local_alignment()

    def _load_local_alignment(self):
        mode = self.current_param.param_dict["alignImages"]["localAlignment"]
        self.local_alignment_filename = (
            self.data_folder.output_files["alignImages"].split(".")[0]
            + "_"
            + mode
            + ".dat"
        )

        if os.path.exists(self.local_alignment_filename):
            self.alignment_results_table = Table.read(
                self.local_alignment_filename, format="ascii.ecsv"
            )
            self.alignment_results_table_read = True

            # builds dict of local alignments
            self.build_local_alignment_dict()

            print_log(
                "$ LocalAlignment file loaded: {}\n$ Will correct coordinates using {} alignment".format(
                    self.local_alignment_filename, mode
                )
            )
            print_log(
                "$ Number of records: {}".format(len(self.alignment_results_table))
            )
        else:
            print_log(
                "\n\n# Warning: could not find localAlignment: {}\n Proceeding with only global alignments...".format(
                    self.local_alignment_filename
                )
            )
            self.alignment_results_table_read = False
            self.alignment_results_table = Table()
            self.dict_error_block_masks = Table()

    def build_local_alignment_dict(self):
        """
        Builds dictionary of local corrections for each ROI, barcode cycle, and block combination

        Parameters
        ----------
        self.alignment_results_table: astropy Table
            alignment_results_table table

        self.alignment_results_table_read: Boolean
            True when alignment_results_table table was read from disk

        Returns
        -------
        exit_code: Boolean

        self.dict_error_block_masks: dict

        """
        if not self.alignment_results_table_read:
            print_log("Did not find alignment_results_table. Cannot continue")
            return False
        else:
            alignment_results_table = self.alignment_results_table

        dict_error_block_masks = {}

        for row in alignment_results_table:
            n_roi = "ROI:" + str(row["ROI #"])
            n_barcode = "barcode:" + row["label"]
            n_block_i = "block_i:" + str(row["block_i"])
            n_block_j = "block_j:" + str(row["block_j"])

            if n_roi not in dict_error_block_masks:
                dict_error_block_masks[n_roi] = {}

            if n_barcode not in dict_error_block_masks[n_roi].keys():
                dict_error_block_masks[n_roi][n_barcode] = {}

            if n_block_i not in dict_error_block_masks[n_roi][n_barcode].keys():
                dict_error_block_masks[n_roi][n_barcode][n_block_i] = {}

            if (
                n_block_j
                not in dict_error_block_masks[n_roi][n_barcode][n_block_i].keys()
            ):
                dict_error_block_masks[n_roi][n_barcode][n_block_i][n_block_j] = {}

            dict_error_block_masks[n_roi][n_barcode][n_block_i][n_block_j] = {
                "shift_z": row["shift_z"],
                "shift_x": row["shift_x"],
                "shift_y": row["shift_y"],
                "quality_xy": row["quality_xy"],
                "quality_zy": row["quality_zy"],
                "quality_zx": row["quality_zx"],
            }

        self.dict_error_block_masks = dict_error_block_masks
        return True

    def register_barcode_map_file(self, file):
        if "3D" in os.path.basename(file):
            self.ndims = 3
        else:
            self.ndims = 2

        # loads barcode coordinate Tables
        table = LocalizationTable()
        barcode_map_full, _ = table.load(
            file
        )  # barcode_map_full, unique_barcodes = table.load(file)

        # checks that everything is OK
        if len(barcode_map_full) < 1:
            print_log(f"\nWARNING>{file} contains an empty table!")
            return None

        if "comments" in barcode_map_full.meta.keys():
            if "registered" in barcode_map_full.meta["comments"]:
                print_log(
                    f"\nWARNING>{file} contains a table thas was already registered! \nWill not do anything"
                )
                return None

        # preserves original copy of table for safe keeping
        new_file = get_file_table_new_name(file)
        table.save(new_file, barcode_map_full)
        barcode_map_full_unregistered = barcode_map_full.copy()

        # indexes table by ROI
        barcode_map_roi, number_rois = table.decode_rois(barcode_map_full)

        for i_roi in range(number_rois):
            # creates sub Table for this ROI
            barcode_map = barcode_map_roi.groups[i_roi]
            n_roi = barcode_map["ROI #"][0]
            print_log(f"\nProcessing barcode localization Table for ROI: {n_roi}")

            # registers localizations
            barcode_map = self.register_barcodes(barcode_map)

        # saves and plots registered barcode coordinate Tables
        table.save(file, barcode_map, comments="registered")
        table.plot_distribution_fluxes(
            barcode_map, [file.split(".")[0], "_registered", "_barcode_stats", ".png"]
        )
        table.plots_localizations(
            barcode_map,
            [file.split(".")[0], "_registered", "_barcode_localizations", ".png"],
        )
        table.compares_localizations(
            barcode_map,
            barcode_map_full_unregistered,
            [file.split(".")[0], "_registered", "_barcode_comparisons", ".png"],
        )

    def register(self):
        """
        Function that registers barcodes using a local drift correction table produced by *alignImages3D*


        Returns
        -------
        None.

        """
        session_name = "register_localizations"

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

        current_folder = self.current_param.param_dict["rootFolder"]
        self.data_folder.create_folders(current_folder, self.current_param)
        print_log("> Processing Folder: {}".format(current_folder))

        # Loads localAlignment if it exists otherwise it exits with error
        self.load_local_alignment()

        if not self.alignment_results_table_read:
            print_log("Unable to find aligment table.\nDid you run alignImages3D?\n\n ")
            sys.exit(
                "ERROR: Expected to find: {}--> Aborting.".format(
                    self.local_alignment_filename
                )
            )

        # iterates over barcode localization tables in the current folder
        files = list(
            glob.glob(
                self.data_folder.output_files["segmentedObjects"]
                + "_*"
                + label
                + ".dat"
            )
        )

        if len(files) < 1:
            print_log("No localization table found to process!")
            return

        for file in files:
            self.register_barcode_map_file(file)

        print_log(f" {len(files)} barcode tables processed in {current_folder}")
