#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:23:36 2020

@author: marcnol

This script:
    - iterates over rois
        - assigns barcode localizations to masks
        - applies local drift correction, if available
        - removes localizations using flux and driftTolerance
        - calculates the pair-wise distances for each single-cell mask
        - outputs are:
            - Table with #cell #PWD #coordinates (e.g. buildsPWDmatrix_3D_order:0_ROI:1.ecsv)
            - NPY array with single cell PWD single cell matrices (e.g. buildsPWDmatrix_3D_HiMscMatrix.npy)
            - NPY array with barcode identities (e.g. buildsPWDmatrix_3D_uniqueBarcodes.ecsv)
            - the files with no "3D" tag contain data analyzed using 2D localizations.

    - Single-cell results are combined together to calculate:
        - Distribution of pairwise distance for each barcode combination
        - Ensemble mean pairwise distance matrix using mean of distribution
        - Ensemble mean pairwise distance matrix using Kernel density estimation
        - Ensemble Hi-M matrix using a predefined threshold
        - For each of these files, there is an image in PNG format saved. Images containing "3D" are for 3D other are for 2D.


"""

# =============================================================================
# IMPORTS
# =============================================================================

import glob
import os
import re
import sys
import uuid

# to remove in a future version
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from photutils.segmentation import SegmentationImage
from sklearn.metrics import pairwise_distances
from tqdm import trange
from tqdm.contrib import tzip

from fileProcessing.fileManagement import Folders, print_log, write_string_to_file
from matrixOperations.HIMmatrixOperations import (
    calculate_contact_probability_matrix,
    plot_distance_histograms,
    plot_matrix,
)

warnings.filterwarnings("ignore")

# =============================================================================
# CLASSES
# =============================================================================


class CellID:
    def __init__(self, param, data_folder, barcode_map_roi, masks, roi, ndims=2):
        self.current_param = param
        self.data_folder = data_folder
        self.barcode_map_roi = barcode_map_roi
        self.masks = masks
        self.n_cells_assigned = 0
        self.n_cells_unassigned = 0
        self.n_barcodes_in_mask = 0
        self.ndims = ndims
        self.dict_error_block_masks = (
            {}
        )  # contains the results from blockAlignment, if existing

        self.segmentation_mask = SegmentationImage(self.masks)
        self.number_masks = self.segmentation_mask.nlabels
        self.roi = roi
        self.alignment_results_table = Table()
        self.alignment_results_table_read = False
        self.barcodes_in_mask = {}
        self.log_name_md = ""
        self.found_match = []

        for mask in range(self.number_masks + 1):
            self.barcodes_in_mask["maskID_" + str(mask)] = []

    def initialize_lists(self):
        (
            self.rois,
            self.cell_id,
            self.n_barcodes,
            self.barcode_ids,
            self.p,
            self.cuid,
            self.buid,
            self.barcode_coordinates,
        ) = (
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
        )

    def filter_localizations__quality(self, i, flux_min):
        """
        [filters barcode localizations either by brigthness or 3D localization accuracy]

        Parameters
        ----------
        i : int
            index in barcode_map Table
        flux_min : float
            Minimum flux to keep barcode localization

        Returns
        -------
        keep : Boolean
            True if the test is passed.

        """
        if "3DfitKeep" in self.barcode_map_roi.groups[0].keys() and self.ndims == 3:
            # [reading the flag in barcode_map_roi assigned by the 3D localization routine]
            keep = (
                self.barcode_map_roi.groups[0]["3DfitKeep"][i]
                and self.barcode_map_roi.groups[0]["flux"][i] > flux_min
            )
        else:
            # [or by reading the flux from 2D localization]
            keep = self.barcode_map_roi.groups[0]["flux"][i] > flux_min

        return keep

    def filter_localizations_block_alignment(self, i, tolerance_drift, block_size):
        """
        [filters barcode per blockAlignmentMask, if existing]
        runs only if localAligment was not run!

        Parameters
        ----------
        i : int
            index in barcode_map Table
        tolerance_drift : float
            tolerance to keep barcode localization, in pixel units
        block_size : int
            size of blocks used for blockAlignment.

        Returns
        -------
        keep_alignment : Boolean
            True if the test is passed.

        """
        y_int = int(self.barcode_map_roi.groups[0]["xcentroid"][i])
        x_int = int(self.barcode_map_roi.groups[0]["ycentroid"][i])
        keep_alignment = True
        if (
            not self.alignment_results_table_read
        ):  # only proceeds if localAlignment was not performed
            barcode_id = "barcode:" + str(self.barcode_map_roi.groups[0]["Barcode #"][i])
            barcode_roi = "ROI:" + str(self.barcode_map_roi.groups[0]["ROI #"][i])

            if len(self.dict_error_block_masks) > 0:
                if barcode_roi in self.dict_error_block_masks.keys():
                    if barcode_id in self.dict_error_block_masks[barcode_roi].keys():
                        error_mask = self.dict_error_block_masks[barcode_roi][barcode_id]
                        keep_alignment = (
                            error_mask[
                                int(np.floor(x_int / block_size)),
                                int(np.floor(y_int / block_size)),
                            ]
                            < tolerance_drift
                        )

            # keeps it always if barcode is fiducial
            if (
                "RT" + str(self.barcode_map_roi.groups[0]["Barcode #"][i])
                in self.current_param.param_dict["alignImages"]["referenceFiducial"]
            ):
                keep_alignment = True

        return keep_alignment

    def plot_distribution_fluxes(self):
        file_name = (
            self.data_folder.output_folders["buildsPWDmatrix"]
            + os.sep
            + "BarcodeStats_ROI:"
            + str(self.n_roi)
            + "_"
            + str(self.ndims)
            + "D.png"
        )

        fig, axes = plt.subplots(1, 2)
        ax = axes.ravel()
        fig.set_size_inches((10, 5))

        fluxes = self.barcode_map_roi.groups[0]["flux"]
        sharpness = self.barcode_map_roi.groups[0]["sharpness"]
        roundness = self.barcode_map_roi.groups[0]["roundness1"]
        peak = self.barcode_map_roi.groups[0]["peak"]
        mag = self.barcode_map_roi.groups[0]["mag"]

        p_1 = ax[0].scatter(fluxes, sharpness, c=peak, cmap="terrain", alpha=0.5)
        ax[0].set_title("color: peak intensity")
        ax[0].set_xlabel("flux")
        ax[0].set_ylabel("sharpness")

        p_2 = ax[1].scatter(roundness, mag, c=peak, cmap="terrain", alpha=0.5)
        ax[1].set_title("color: peak intensity")
        ax[1].set_xlabel("roundness")
        ax[1].set_ylabel("magnitude")
        fig.colorbar(p_2, ax=ax[1], fraction=0.046, pad=0.04)

        fig.savefig(file_name)

        plt.close(fig)

        write_string_to_file(
            self.log_name_md,
            "Barcode stats for ROI:{}, dims:{} \n![]({})\n".format(
                self.n_roi, self.ndims, file_name
            ),
            "a",
        )

    def plots_barcodes_alignment(self, block_size):
        """
        plots barcode localizations together with the blockAlignment map

        Returns
        -------
        None.

        """
        file_name = (
            self.data_folder.output_folders["buildsPWDmatrix"]
            + os.sep
            + "BarcodeAlignmentAccuracy_ROI:"
            + str(self.n_roi)
            + "_"
            + str(self.ndims)
            + "D.png"
        )

        fig, axes = plt.subplots()
        fig.set_size_inches((20, 20))

        accuracy, x, y = [], [], []
        print_log("> Plotting barcode alignments...")
        for i in trange(len(self.barcode_map_roi.groups[0])):
            barcode_id = "barcode:" + str(self.barcode_map_roi.groups[0]["Barcode #"][i])
            barcode_roi = "ROI:" + str(self.barcode_map_roi.groups[0]["ROI #"][i])
            y_int = int(self.barcode_map_roi.groups[0]["xcentroid"][i])
            x_int = int(self.barcode_map_roi.groups[0]["ycentroid"][i])

            if len(self.dict_error_block_masks) > 0:
                if barcode_roi in self.dict_error_block_masks.keys():
                    if barcode_id in self.dict_error_block_masks[barcode_roi].keys():
                        error_mask = self.dict_error_block_masks[barcode_roi][barcode_id]
                        accuracy.append(
                            error_mask[
                                int(np.floor(x_int / block_size)),
                                int(np.floor(y_int / block_size)),
                            ]
                        )
                        x.append(self.barcode_map_roi.groups[0]["xcentroid"][i])
                        y.append(self.barcode_map_roi.groups[0]["ycentroid"][i])

        p_1 = axes.scatter(
            x, y, s=5, c=accuracy, cmap="terrain", alpha=0.5, vmin=0, vmax=5
        )
        fig.colorbar(p_1, ax=axes, fraction=0.046, pad=0.04)
        axes.set_title("barcode drift correction accuracy, px")

        axes.axis("off")

        fig.savefig(file_name)

        plt.close(fig)

        write_string_to_file(
            self.log_name_md,
            "Barcode stats for ROI:{}, dims:{} \n![]({})\n".format(
                self.n_roi, self.ndims, file_name
            ),
            "a",
        )

    def align_by_masking(self):
        """
        Assigns barcodes to masks and creates <n_barcodes_in_mask>
        This routine will only select which barcodes go to each cell mask

        Returns
        -------
        self.barcodes_in_mask # dictionnary with the identities of barcodes contained in each mask.
            Keys: 'maskID_1', 'maskID_2', and so on

        self.n_barcodes_in_mask # vector containing the number of barcodes for each mask
        self.n_cells_assigned # number of cells assigned
        self.n_cells_unassigned # number of cells unassigned
        """

        n_barcodes_in_mask = np.zeros(self.number_masks + 2)
        flux_key = "flux_min_3D" if self.ndims == 3 else "flux_min"
        if flux_key in self.current_param.param_dict["buildsPWDmatrix"]:
            flux_min = self.current_param.param_dict["buildsPWDmatrix"][flux_key]
        else:
            flux_min = 0
            print_log("# Flux min not found. Set to {}!".format(flux_min))

        if "toleranceDrift" in self.current_param.param_dict["buildsPWDmatrix"]:
            tolerance_drift = self.current_param.param_dict["buildsPWDmatrix"][
                "toleranceDrift"
            ]
        else:
            tolerance_drift = 1
            print_log("# toleranceDrift not found. Set to {}!".format(tolerance_drift))

        if "blockSize" in self.current_param.param_dict["alignImages"]:
            block_size = self.current_param.param_dict["alignImages"]["blockSize"]
        else:
            block_size = 256
            print_log("# blockSize not found. Set to {}!".format(block_size))

        print_log(
            "\n$ ndims = {}\n$ Flux min = {} \n$ ToleranceDrift = {} px\n$ Reference barcode = {}".format(
                self.ndims,
                flux_min,
                tolerance_drift,
                self.current_param.param_dict["alignImages"]["referenceFiducial"],
            )
        )

        # Produces images of distribution of fluxes.
        self.plot_distribution_fluxes()
        self.plots_barcodes_alignment(block_size)

        keep_quality_all, keep_alignment_all, n_barcodes_roi = [], [], 0
        # loops over barcode Table rows in a given ROI
        print_log("> Aligning by masking...")
        for i in trange(
            len(self.barcode_map_roi.groups[0])
        ):  # i is the index of the barcode in barcode_map_roi
            barcode = self.barcode_map_roi.groups[0]["Barcode #"][i]
            roi = self.barcode_map_roi.groups[0]["ROI #"][i]

            # [filters barcode localizations either by]
            keep_quality = self.filter_localizations__quality(i, flux_min)

            # [filters barcode per blockAlignmentMask, if existing]
            keep_alignment = self.filter_localizations_block_alignment(
                i, tolerance_drift, block_size
            )

            # applies all filters
            if keep_quality and keep_alignment:

                # keeps the particle if the test passed
                x_uncorrected = self.barcode_map_roi.groups[0]["ycentroid"][
                    i
                ]  # control inversion between x-y
                y_uncorrected = self.barcode_map_roi.groups[0]["xcentroid"][i]

                if self.ndims == 2:
                    z_uncorrected = self.barcode_map_roi.groups[0]["zcentroid"][i] = 0.0
                else:
                    z_uncorrected = self.barcode_map_roi.groups[0]["zcentroid"][i]

                y_int = int(y_uncorrected)
                x_int = int(x_uncorrected)

                # finds what mask label this barcode is sitting on
                mask_id = self.masks[x_int][y_int]

                # Corrects XYZ coordinate of barcode if localDriftCorrection is available
                zxy_uncorrected = [z_uncorrected, x_uncorrected, y_uncorrected]
                rt_barcode = "RT" + str(barcode)
                if (
                    rt_barcode
                    not in self.current_param.param_dict["alignImages"][
                        "referenceFiducial"
                    ]
                ):
                    zxy_corrected = self.search_local_shift(
                        roi, mask_id, barcode, zxy_uncorrected, tolerance_drift
                    )
                else:
                    # if it is the reference cycle, then it does not correct coordinates
                    zxy_corrected = zxy_uncorrected

                # rewrites corrected XYZ values to Table
                self.barcode_map_roi.groups[0]["ycentroid"][i] = zxy_corrected[1]
                self.barcode_map_roi.groups[0]["xcentroid"][i] = zxy_corrected[2]
                if self.ndims > 2:
                    self.barcode_map_roi.groups[0]["zcentroid"][i] = zxy_corrected[0]

                # attributes CellID to a barcode
                self.barcode_map_roi["CellID #"][i] = mask_id

                # if it is not background,
                if mask_id > 0:
                    # increments counter of number of barcodes in the cell mask attributed
                    n_barcodes_in_mask[mask_id] += 1

                    # stores the identify of the barcode to the mask
                    self.barcodes_in_mask["maskID_" + str(mask_id)].append(i)

            # keeps statistics
            if int(self.barcode_map_roi.groups[0]["ROI #"][i]) == int(self.n_roi):
                keep_quality_all.append(keep_quality)
                keep_alignment_all.append(keep_alignment)
                n_barcodes_roi += 1

        # Total number of masks assigned and not assigned
        self.n_cells_assigned = np.count_nonzero(n_barcodes_in_mask > 0)
        self.n_cells_unassigned = self.number_masks - self.n_cells_assigned

        # this list contains which barcodes are allocated to which masks
        self.n_barcodes_in_mask = n_barcodes_in_mask

        print_log(
            "$ Number of localizations passing quality test: {} / {}".format(
                sum(keep_quality_all), n_barcodes_roi
            )
        )

        print_log(
            "$ Number of localizations passing alignment test: {} / {}".format(
                sum(keep_alignment_all), n_barcodes_roi
            )
        )

        print_log(
            "$ Number of cells assigned: {} | discarded: {}".format(
                self.n_cells_assigned, self.n_cells_unassigned
            )
        )

    def search_local_shift(self, roi, cell_id, barcode, zxy_uncorrected, tolerance_drift=1):

        if "mask2D" in self.current_param.param_dict["alignImages"]["localAlignment"]:
            return self.search_local_shift_mask_2d(roi, cell_id, zxy_uncorrected)
        elif (
            "block3D" in self.current_param.param_dict["alignImages"]["localAlignment"]
            and self.alignment_results_table_read
        ):
            return self.search_local_shift_block_3d(
                roi, barcode, zxy_uncorrected, tolerance_drift
            )
        else:  # no correction was applied because the localAlignmentTable was not found
            return zxy_uncorrected

    def search_local_shift_block_3d(self, roi, barcode, zxy_uncorrected, tolerance_drift=1):
        """
        Searches for local drift for a specific barcode in a given ROI.
        If it exists then it adds to the uncorrected coordinates

        Parameters
        ----------
        roi : string
            roi used
        cell_id: string
            ID of the cell
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
        _found_match = False

        # gets blockSize
        block_size_xy = self.alignment_results_table[0]["blockXY"]

        # zxy coord in block reference coord system
        zxy_block = [np.floor(a / block_size_xy).astype(int) for a in zxy_uncorrected]

        for row in self.alignment_results_table:

            # I need to check that the XY coordinates from localization are the same as the ij indices from the block decomposition!

            if (
                row["ROI #"] == roi
                and row["label"] == "RT" + str(barcode)
                and row["block_i"] == zxy_block[1]
                and row["block_j"] == zxy_block[2]
            ):
                _found_match = True
                shifts = [row["shift_z"], row["shift_x"], row["shift_y"]]

                # checks that drifts > tolerance_drift are not applied
                if max(shifts) < tolerance_drift:
                    zxy_corrected = [
                        a + shift for a, shift in zip(zxy_uncorrected, shifts)
                    ]
                else:
                    zxy_corrected = zxy_uncorrected

                # check for quality of shift correction before applying it !!
                #!TODO

        # keeps uncorrected values if no match is found
        if not _found_match:
            print_log(
                "# Did not find match for ROI #{} barcode #{}".format(roi, barcode)
            )
            zxy_corrected = zxy_uncorrected
            self.found_match.append(False)
        else:
            self.found_match.append(True)

        return zxy_corrected

    def search_local_shift_mask_2d(self, roi, cell_id, zxy_uncorrected):
        """
        Searches for local drift for current mask. If it exists then id adds it to the uncorrected coordinates

        Parameters
        ----------
        roi : string
            roi used
        cell_id: string
            ID of the cell
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

        _found_match = False
        for row in self.alignment_results_table:
            if row["ROI #"] == roi and row["CellID #"] == cell_id:
                _found_match = True
                shifts = [0, row["shift_x"], row["shift_y"]]
                zxy_corrected = [a + shift for a, shift in zip(zxy_uncorrected, shifts)]

        # keeps uncorrected values if no match is found
        if not _found_match:
            print_log(
                "# Did not find match for CellID #{} in ROI #{}".format(cell_id, roi)
            )
            zxy_corrected = zxy_uncorrected
            self.found_match.append(False)
        else:
            self.found_match.append(True)

        return zxy_corrected

    def build_vector(self, group_keys, x, y, z):
        """
        Builds vector from coordinates

        Parameters
        ----------
        group_keys : list
            list of headers in the barcodes table
        x : float
            x coordinates
        y : float
            y coordinates
        z : float
            z coordinates

        Returns
        -------
        coords : np array
            vector with coordinates in nanometers.

        """

        coords = np.column_stack(
            (x * self.pixel_size["x"], y * self.pixel_size["y"], z * self.pixel_size["z"])
        )

        return coords

    def calculate_pwd_single_mask(self, roi, cell_id, group_keys, x, y, z):
        """
        Calculates PWD between barcodes detected in a given mask. For this:
            - converts xyz pixel coordinates into nm using self.pixel_size dictionary
            - calculates pair-wise distance matrix in nm
            - converts it into pixel units using self.pixel_size['x'] as an isotropic pixelsize.

        Parameters
        ----------
        roi : string
            roi used
        cell_id: string
            ID of the cell
        x_uncorrected: float
            x coordinates uncorrected
        y_uncorrected: float
            y coordinates uncorrected
        z_uncorrected: float
            z coordinates uncorrected

        Returns
        -------
        Returns pairwise distance matrix between corrected barcodes in isotropic pixel units

        """
        r_nm = self.build_vector(group_keys, x, y, z)

        pairwd = pairwise_distances(r_nm)

        pairwd = pairwd / self.pixel_size["x"]

        return pairwd

    def builds_sc_distance_table(self):
        """
        iterates over all masks, calculates PWD for each mask, assigns them to sc_distance_table

        Returns
        -------
        sc_distance_table

        """
        # sorts Table by cellID
        barcode_map_roi = self.barcode_map_roi
        barcode_map_roi_cell_id = barcode_map_roi.group_by(
            "CellID #"
        )  # ROI data sorted by cellID

        self.initialize_lists()

        # iterates over all cell masks in an ROI
        print_log("> Building sc_ distance Tables")
        for key, group in tzip(
            barcode_map_roi_cell_id.groups.keys, barcode_map_roi_cell_id.groups
        ):
            if key["CellID #"] > 1:  # excludes cellID 0 as this is background

                group_keys, cell_id, roi = (
                    group.keys(),
                    key["CellID #"],
                    group["ROI #"].data[0],
                )

                # gets lists of x, y and z coordinates for barcodes assigned to a cell mask
                x, y, z = (
                    np.array(group["xcentroid"].data),
                    np.array(group["ycentroid"].data),
                    np.array(group["zcentroid"].data),
                )

                # calculates the pwd between barcodes in CellID
                pwd = self.calculate_pwd_single_mask(roi, cell_id, group_keys, x, y, z)
                r_nm = self.build_vector(group_keys, x, y, z)

                self.rois.append(group["ROI #"].data[0])
                self.cell_id.append(key["CellID #"])
                self.n_barcodes.append(len(group))
                self.barcode_ids.append(group["Barcode #"].data)
                self.buid.append(group["Buid"].data)
                self.p.append(pwd)
                self.barcode_coordinates.append(r_nm)
                self.cuid.append(str(uuid.uuid4()))  # creates cell unique identifier

        print_log(
            "$ Local correction applied to {}/{} barcodes in ROI {}".format(
                np.nonzero(self.found_match)[0].shape[0],
                len(self.found_match),
                group["ROI #"].data[0],
            )
        )

        print_log("$ Coordinates dimensions: {}".format(self.ndims))

        sc_distance_table = Table()
        sc_distance_table["Cuid"] = self.cuid
        sc_distance_table["ROI #"] = self.rois
        sc_distance_table["CellID #"] = self.cell_id
        sc_distance_table["nBarcodes"] = self.n_barcodes
        sc_distance_table["Barcode #"] = self.barcode_ids
        sc_distance_table["Buid"] = self.buid
        sc_distance_table["PWDmatrix"] = self.p
        sc_distance_table["barcode xyz, nm"] = self.barcode_coordinates

        self.sc_distance_table = sc_distance_table

    def build_distance_matrix(self, mode="mean"):
        """
        Builds pairwise distance matrix from a coordinates table

        Parameters
        ----------
        mode : string, optional
            The default is "mean": calculates the mean distance if there are several combinations possible.
            "min": calculates the minimum distance if there are several combinations possible.
            "last": keeps the last distance calculated

        Returns
        -------
        self.sc_matrix the single-cell PWD matrix
        self.mean_sc_matrix the ensamble PWD matrix (mean of sc_matrix without nans)
        self.unique_barcodes list of unique barcodes

        """
        # [ builds sc_distance_table ]
        self.builds_sc_distance_table()
        print_log("$ Cells with barcodes found: {}".format(len(self.sc_distance_table)))

        # [ builds sc_matrix ]
        number_matrices = len(self.sc_distance_table)  # z dimensions of sc_matrix

        # unique_barcodes = np.unique(self.barcode_map_roi["Barcode #"].data)
        unique_barcodes = self.unique_barcodes

        # number of unique Barcodes for xy dimensions of sc_matrix
        number_unique_barcodes = unique_barcodes.shape[0]
        sc_matrix = np.zeros(
            (number_unique_barcodes, number_unique_barcodes, number_matrices)
        )
        sc_matrix[:] = np.NaN

        # loops over cell masks
        for i_cell, sc_pwd_item in zip(range(number_matrices), self.sc_distance_table):
            barcodes_to_process = sc_pwd_item["Barcode #"]

            # loops over barcodes detected in cell mask: barcode1
            for barcode1, ibarcode1 in zip(
                barcodes_to_process, range(len(barcodes_to_process))
            ):
                index_barcode_1 = np.nonzero(unique_barcodes == barcode1)[0][0]

                # loops over barcodes detected in cell mask: barcode2
                for barcode2, ibarcode2 in zip(
                    barcodes_to_process, range(len(barcodes_to_process))
                ):
                    index_barcode_2 = np.nonzero(unique_barcodes == barcode2)[0][0]

                    if barcode1 != barcode2:

                        # attributes distance from the PWDmatrix field in the sc_pwd_item table
                        newdistance = sc_pwd_item["PWDmatrix"][ibarcode1][ibarcode2]

                        # inserts value into sc_matrix
                        if mode == "last":
                            sc_matrix[index_barcode_1][index_barcode_2][i_cell] = newdistance
                        elif mode == "mean":
                            sc_matrix[index_barcode_1][index_barcode_2][i_cell] = np.nanmean(
                                [
                                    newdistance,
                                    sc_matrix[index_barcode_1][index_barcode_2][i_cell],
                                ]
                            )
                        elif mode == "min":
                            sc_matrix[index_barcode_1][index_barcode_2][i_cell] = np.nanmin(
                                [
                                    newdistance,
                                    sc_matrix[index_barcode_1][index_barcode_2][i_cell],
                                ]
                            )

        self.sc_matrix = sc_matrix
        self.mean_sc_matrix = np.nanmean(sc_matrix, axis=2)
        # self.unique_barcodes = unique_barcodes


# =============================================================================
# FUNCTIONS
# =============================================================================


def calculate_n_matrix(sc_matrix):

    number_cells = sc_matrix.shape[2]

    if number_cells > 0:
        n_matrix = np.sum(~np.isnan(sc_matrix), axis=2)
    else:
        number_barcodes = sc_matrix.shape[0]
        n_matrix = np.zeros((number_barcodes, number_barcodes))

    return n_matrix


def load_local_alignment(current_param, data_folder):

    if "None" in current_param.param_dict["alignImages"]["localAlignment"]:
        print_log(
            "\n\n$ localAlignment option set to {}".format(
                current_param.param_dict["alignImages"]["localAlignment"]
            )
        )
        return False, Table()
    else:
        return _load_local_alignment(
            data_folder, current_param.param_dict["alignImages"]["localAlignment"]
        )


def _load_local_alignment(data_folder, mode):

    local_alignment_filename = (
        data_folder.output_files["alignImages"].split(".")[0] + "_" + mode + ".dat"
    )
    if os.path.exists(local_alignment_filename):
        alignment_results_table = Table.read(
            local_alignment_filename, format="ascii.ecsv"
        )
        alignment_results_table_read = True
        print_log(
            "$ LocalAlignment file loaded: {}\n$ Will correct coordinates using {} alignment".format(
                local_alignment_filename, mode
            )
        )
        print_log("$ Number of records: {}".format(len(alignment_results_table)))
    else:
        print_log(
            "\n\n# Warning: could not find localAlignment: {}\n Proceeding with only global alignments...".format(
                local_alignment_filename
            )
        )
        alignment_results_table_read = False
        alignment_results_table = Table()

    return alignment_results_table, alignment_results_table_read


# def load_local_alignment(data_folder):
#     """
#     reads and returns localAlignmentTable, if it exists

#     Parameters
#     ----------
#     data_folder : folder()
#         DESCRIPTION.

#     Returns
#     -------
#     alignment_results_table : Table()
#         DESCRIPTION.
#     alignment_results_table_read : Boolean
#         DESCRIPTION.

#     """
#     local_alignment_filename = data_folder.output_files["alignImages"].split(".")[0] + "_localAlignment.dat"
#     if os.path.exists(local_alignment_filename):
#         alignment_results_table = Table.read(local_alignment_filename, format="ascii.ecsv")
#         alignment_results_table_read = True
#         print_log("LocalAlignment file loaded !\nWill correct coordinates in XY")
#     else:
#         print_log(
#             "\n\n*** Warning: could not find localAlignment: {}\n Proceeding with only global alignments...".format(
#                 local_alignment_filename
#             )
#         )
#         alignment_results_table_read = False
#         alignment_results_table = Table()

#     return alignment_results_table, alignment_results_table_read


def load_barcode_map(filename_barcode_coordinates, ndims):
    """
    Loads barcode_map

    Parameters
    ----------
    filename_barcode_coordinates : string
        filename with barcode_map
    ndims : int
        either 2 or 3.

    Returns
    -------
    barcode_map : Table()
    localization_dimension : int
        either 2 or 3.
    unique_barcodes: list
        lis of unique barcodes read from barcode_map

    """
    if os.path.exists(filename_barcode_coordinates):
        barcode_map = Table.read(filename_barcode_coordinates, format="ascii.ecsv")
        print_log(
            "$ Successfully loaded barcode localizations file: {}".format(
                filename_barcode_coordinates
            )
        )

        unique_barcodes = np.unique(barcode_map["Barcode #"].data)
        number_unique_barcodes = unique_barcodes.shape[0]

        print_log(
            "Number Barcodes read from barcode_map: {}".format(number_unique_barcodes)
        )
        print_log("Unique Barcodes detected: {}".format(unique_barcodes))
    else:
        print_log(
            "\n\n# ERROR: could not find coordinates file: {}".format(
                filename_barcode_coordinates
            )
        )
        sys.exit()

    return barcode_map, ndims, unique_barcodes


def build_dictionary_error_alignment_masks(current_param, data_folder):
    """
    Builds and returns dictionary with error alignment block masks produced during the alignment process if
    the 'blockAlignment' option was used

    Parameters
    ----------
    current_param : Parameters()
    data_folder : folder()

    Returns
    -------
    dict_error_block_masks : dict

    """
    folder = data_folder.output_folders["alignImages"]
    file_list = glob.glob(folder + os.sep + "*_errorAlignmentBlockMap.npy")

    # decodes files and builds dictionnary
    filename_regexp = current_param.param_dict["acquisition"]["fileNameRegExp"]
    filename_regexp = filename_regexp.split(".")[0]
    list_re = [
        re.search(
            filename_regexp, os.path.basename(x).split("_errorAlignmentBlockMap.npy")[0]
        )
        for x in file_list
    ]

    dict_error_block_masks = {}

    for file, reg_exp in zip(file_list, list_re):
        if "ROI:" + str(int(reg_exp["roi"])) not in dict_error_block_masks.keys():
            dict_error_block_masks["ROI:" + str(int(reg_exp["roi"]))] = {}
        if (
            "barcode:" + reg_exp["cycle"].split("RT")[-1]
            not in dict_error_block_masks.keys()
        ):
            new_mask = np.load(file)
            dict_error_block_masks["ROI:" + str(int(reg_exp["roi"]))][
                "barcode:" + reg_exp["cycle"].split("RT")[-1]
            ] = new_mask

    return dict_error_block_masks


def plots_all_matrices(
    sc_matrix_collated,
    n_matrix,
    unique_barcodes,
    pixel_size,
    number_rois,
    output_filename,
    log_name_md,
    localization_dimension,
):
    """
    Plots all matrices after analysis

    Parameters
    ----------
    sc_matrix_collated : npy array
        pwd matrix for single cells.
    n_matrix : npy array
        2d matrix with number of measurements per barcode combination.
    unique_barcodes : npy array
        barcode identities.
    pixel_size : npy array
        pixelsize in um.
    number_rois : int
        self explanatory.
    output_filename : str
        self explanatory.
    log_name_md : str
        Markdown filename.
    localization_dimension : int
        indicates dimension of barcode localization.

    Returns
    -------
    None.

    """
    # adapts clim depending on whether 2 or 3 dimensions are used for barcode localizations
    if localization_dimension == 2:
        clim = 1.6
    else:
        clim = 2.2

    # plots PWD matrix
    # uses KDE
    plot_matrix(
        sc_matrix_collated,
        unique_barcodes,
        pixel_size,
        number_rois,
        output_filename,
        log_name_md,
        figtitle="PWD matrix - KDE",
        mode="KDE",  # median or KDE
        clim=clim,
        c_m="terrain",
        filename_ending="_PWDmatrixKDE.png",
    )  # need to validate use of KDE. For the moment it does not handle well null distributions

    # uses median
    plot_matrix(
        sc_matrix_collated,
        unique_barcodes,
        pixel_size,
        number_rois,
        output_filename,
        log_name_md,
        figtitle="PWD matrix - median",
        mode="median",  # median or KDE
        clim=clim,
        c_m="coolwarm",
        filename_ending="_PWDmatrixMedian.png",
    )  # need to validate use of KDE. For the moment it does not handle well null distributions

    # calculates and plots contact probability matrix from merged samples/datasets
    him_matrix, n_cells = calculate_contact_probability_matrix(
        sc_matrix_collated, unique_barcodes, pixel_size, norm="nonNANs",
    )  # norm: n_cells (default), nonNANs

    c_scale = him_matrix.max()
    plot_matrix(
        him_matrix,
        unique_barcodes,
        pixel_size,
        number_rois,
        output_filename,
        log_name_md,
        figtitle="Hi-M matrix",
        mode="counts",
        clim=c_scale,
        c_m="coolwarm",
        filename_ending="_HiMmatrix.png",
    )

    # plots n_matrix
    plot_matrix(
        n_matrix,
        unique_barcodes,
        pixel_size,
        number_rois,
        output_filename,
        log_name_md,
        figtitle="N-matrix",
        mode="counts",
        clim=np.max(n_matrix),
        c_m="Blues",
        filename_ending="_Nmatrix.png",
    )

    plot_distance_histograms(
        sc_matrix_collated,
        pixel_size,
        output_filename,
        log_name_md,
        mode="KDE",
        kernel_width=0.25,
        optimize_kernel_width=False,
    )


def build_pwd_matrix(
    current_param,
    current_folder,
    filename_barcode_coordinates,
    output_filename,
    data_folder,
    pixel_size={"x": 0.1, "y": 0.1, "z": 0.0},
    log_name_md="log.md",
    ndims=2,
    mask_identifier="DAPI",
):
    """
    Main function that:
        loads and processes barcode localization files, local alignment file, and masks
        initializes <cell_roi> class and assigns barcode localizations to masks
        then constructs the single cell PWD matrix and outputs it toghether with the contact map and the N-map.

    Parameters
    ----------
    current_param : Parameters Class
    current_folder : string
    filename_barcode_coordinates : string
    output_filename : string
    data_folder : Folder Class
        information to find barcode localizations, local drift corrections and masks

    pixel_size : dict, optional
        pixel_size = {'x': pixelSizeXY,
                    'y': pixelSizeXY,
                    'z': pixel_size_z}
        The default is 0.1 for x and y, 0.0 for z. Pixelsize in um

    log_name_md : str, optional
        Filename of Markdown output. The default is "log.md".
    ndims : int, optional
        indicates whether barcodes were localized in 2 or 3D. The default is 2.

    Returns
    -------
    None.

    """
    # Loads localAlignment if it exists
    alignment_results_table, alignment_results_table_read = load_local_alignment(
        current_param, data_folder
    )

    # Loads coordinate Tables
    barcode_map, localization_dimension, unique_barcodes = load_barcode_map(
        filename_barcode_coordinates, ndims
    )

    # Builds dictionnary with filenames of errorAlignmentBlockMasks for each ROI and each barcode
    dict_error_block_masks = build_dictionary_error_alignment_masks(
        current_param, data_folder
    )

    # processes tables
    barcode_map_roi = barcode_map.group_by("ROI #")
    number_rois = len(barcode_map_roi.groups.keys)
    print_log("\n$ rois detected: {}".format(number_rois))

    # loops over rois
    files_in_folder = glob.glob(current_folder + os.sep + "*.tif")
    # sc_matrix_collated, unique_barcodes, processing_order = [], [], 0
    sc_matrix_collated, processing_order = [], 0

    for roi in range(number_rois):
        n_roi = barcode_map_roi.groups.keys[roi][
            0
        ]  # need to iterate over the first index

        print_log(
            "----------------------------------------------------------------------"
        )
        print_log(
            "> Loading masks and pre-processing barcodes for Mask <{}> ROI# {}".format(
                mask_identifier, n_roi
            )
        )
        print_log(
            "----------------------------------------------------------------------"
        )

        barcode_map_single_roi = barcode_map.group_by("ROI #").groups[roi]

        # finds file with cell masks
        files_to_process = [
            file
            for file in files_in_folder
            if file.split("_")[-1].split(".")[0]
            == current_param.param_dict["acquisition"][
                "label_channel"
            ]  # typically "ch00"
            and mask_identifier in os.path.basename(file).split("_")
            and int(os.path.basename(file).split("_")[3]) == n_roi
        ]

        if len(files_to_process) > 0:

            # loads file with cell masks
            filename_roi_masks = (
                os.path.basename(files_to_process[0]).split(".")[0] + "_Masks.npy"
            )
            full_filename_roi_masks = (
                os.path.dirname(filename_barcode_coordinates) + os.sep + filename_roi_masks
            )
            if os.path.exists(full_filename_roi_masks):
                masks = np.load(full_filename_roi_masks)

                # Assigns barcodes to masks for a given ROI
                cell_roi = CellID(
                    current_param,
                    data_folder,
                    barcode_map_single_roi,
                    masks,
                    roi,
                    ndims=localization_dimension,
                )
                cell_roi.ndims, cell_roi.n_roi, cell_roi.log_name_md, cell_roi.pixel_size = (
                    ndims,
                    n_roi,
                    log_name_md,
                    pixel_size,
                )
                cell_roi.unique_barcodes = unique_barcodes

                if alignment_results_table_read:
                    cell_roi.alignment_results_table = alignment_results_table

                cell_roi.dict_error_block_masks = dict_error_block_masks
                cell_roi.alignment_results_table_read = alignment_results_table_read

                # finds what barcodes are in each cell mask
                cell_roi.align_by_masking()

                # builds the single cell distance Matrix
                cell_roi.build_distance_matrix("min")  # mean min last

                print_log(
                    "$ ROI: {}, N cells assigned: {} out of {}\n".format(
                        roi, cell_roi.n_cells_assigned - 1, cell_roi.number_masks
                    )
                )

                # saves Table with results per roi
                cell_roi.sc_distance_table.write(
                    output_filename
                    + "_order:"
                    + str(processing_order)
                    + "_ROI:"
                    + str(n_roi)
                    + ".ecsv",
                    format="ascii.ecsv",
                    overwrite=True,
                )

                if len(sc_matrix_collated) > 0:
                    sc_matrix_collated = np.concatenate(
                        (sc_matrix_collated, cell_roi.sc_matrix), axis=2
                    )
                else:
                    sc_matrix_collated = cell_roi.sc_matrix
                del cell_roi

                processing_order += 1

            # Could not find a file with masks to assign. Report and continue with next roi
            ###############################################################################
            else:
                print_log(
                    "# Error, no mask file found for ROI: {}, segmentedMasks: {}\n".format(
                        n_roi, filename_barcode_coordinates
                    )
                )
                print_log("# File I was searching for: {}".format(full_filename_roi_masks))
                print_log("# Debug: ")
                for file in files_in_folder:
                    if (
                        file.split("_")[-1].split(".")[0]
                        == current_param.param_dict["acquisition"][
                            "label_channel"
                        ]  # typically "ch00"
                        and mask_identifier in file.split("_")
                        and int(os.path.basename(file).split("_")[3]) == n_roi
                    ):
                        print_log("$ Hit found!")
                    print_log(
                        "fileSplit:{}, {} in filename: {}, ROI: {}".format(
                            file.split("_")[-1].split(".")[0],
                            mask_identifier,
                            mask_identifier in os.path.basename(file).split("_"),
                            int(os.path.basename(file).split("_")[3]),
                        )
                    )

    if processing_order > 0:
        # calculates N-matrix: number of PWD distances for each barcode combination
        n_matrix = calculate_n_matrix(sc_matrix_collated)

        # saves output
        np.save(
            output_filename + "_" + mask_identifier + "_HiMscMatrix.npy",
            sc_matrix_collated,
        )
        np.savetxt(
            output_filename + "_" + mask_identifier + "_uniqueBarcodes.ecsv",
            unique_barcodes,
            delimiter=" ",
            fmt="%d",
        )
        np.save(output_filename + "_" + mask_identifier + "_Nmatrix.npy", n_matrix)
        pixel_size_xy = pixel_size["x"]

        if sc_matrix_collated.shape[2] > 0:
            #################################
            # makes and saves outputs plots #
            #################################
            plots_all_matrices(
                sc_matrix_collated,
                n_matrix,
                unique_barcodes,
                pixel_size_xy,
                number_rois,
                output_filename + "_" + mask_identifier,
                log_name_md,
                localization_dimension,
            )
        else:
            print_log(
                "# Nothing to plot. Single cell matrix is empty. Number of cells: {}".format(
                    sc_matrix_collated.shape[2]
                )
            )


def process_pwd_matrices(current_param, current_session):
    """
    Function that assigns barcode localizations to masks and constructs single cell cummulative PWD matrix.

    Parameters
    ----------
    current_param : class
        Parameters
    current_log : class
        logging class.
    current_session : class
        session information

    Returns
    -------
    None.

    """
    session_name = "buildsPWDmatrix"

    # processes folders and files
    data_folder = Folders(current_param.param_dict["rootFolder"])
    print_log("\n===================={}====================\n".format(session_name))
    print_log("$ folders read: {}".format(len(data_folder.list_folders)))
    write_string_to_file(
        current_param.param_dict["fileNameMD"],
        "## {}\n".format(session_name),
        "a",
    )
    label = "barcode"

    for current_folder in data_folder.list_folders:
        # files_folder=glob.glob(current_folder+os.sep+'*.tif')
        data_folder.create_folders(current_folder, current_param)
        print_log("> Processing Folder: {}".format(current_folder))

        available_masks = current_param.param_dict["buildsPWDmatrix"]["masks2process"]
        print_log("> Masks labels: {}".format(available_masks))

        for mask_label in available_masks.keys():

            mask_identifier = available_masks[mask_label]

            filename_barcode_coordinates = (
                data_folder.output_files["segmentedObjects"] + "_" + label + ".dat"
            )
            if os.path.exists(filename_barcode_coordinates):
                # 2D
                output_filename = data_folder.output_files["buildsPWDmatrix"]
                print_log("> 2D processing: {}".format(output_filename))

                if "pixelSizeXY" in current_param.param_dict["acquisition"].keys():
                    pixel_size_xy = current_param.param_dict["acquisition"]["pixelSizeXY"]
                    pixel_size = {"x": pixel_size_xy, "y": pixel_size_xy, "z": 0.0}
                else:
                    pixel_size = {"x": 0.1, "y": 0.1, "z": 0.0}

                build_pwd_matrix(
                    current_param,
                    current_folder,
                    filename_barcode_coordinates,
                    output_filename,
                    data_folder,
                    pixel_size,
                    current_param.param_dict["fileNameMD"],
                    mask_identifier=mask_identifier,
                )

            # 3D
            filename_barcode_coordinates = (
                data_folder.output_files["segmentedObjects"] + "_3D_" + label + ".dat"
            )
            if os.path.exists(filename_barcode_coordinates):
                output_filename = data_folder.output_files["buildsPWDmatrix"] + "_3D"
                print_log("> 3D processing: {}".format(output_filename))

                if (
                    "pixelSizeZ" in current_param.param_dict["acquisition"].keys()
                ) and ("pixelSizeXY" in current_param.param_dict["acquisition"].keys()):
                    pixel_size_xy = current_param.param_dict["acquisition"]["pixelSizeXY"]

                    if "zBinning" in current_param.param_dict["acquisition"]:
                        z_binning = current_param.param_dict["acquisition"]["zBinning"]
                    else:
                        z_binning = 1

                    pixel_size_z = (
                        z_binning * current_param.param_dict["acquisition"]["pixelSizeZ"]
                    )

                    pixel_size = {
                        "x": pixel_size_xy,
                        "y": pixel_size_xy,
                        "z": pixel_size_z * z_binning,
                    }
                else:
                    pixel_size = {"x": 0.1, "y": 0.1, "z": 0.25}

                build_pwd_matrix(
                    current_param,
                    current_folder,
                    filename_barcode_coordinates,
                    output_filename,
                    data_folder,
                    pixel_size,
                    current_param.param_dict["fileNameMD"],
                    ndims=3,
                    mask_identifier=mask_identifier,
                )

            # tights loose ends
            current_session.add(current_folder, session_name)

            print_log("HiM matrix in {} processed".format(current_folder), "info")
