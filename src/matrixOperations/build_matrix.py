#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 08:40:30 2022

@author: marcnol

This script:
    - iterates over chromatin traces
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

import ast
import glob
import os
import re
import sys
import uuid

# to remove in a future version
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, unique
from photutils.segmentation import SegmentationImage
from sklearn.metrics import pairwise_distances
from tqdm import trange
from tqdm.contrib import tzip

from fileProcessing.fileManagement import (
    get_dictionary_value,
    print_log,
)
from matrixOperations.build_traces import initialize_module
from matrixOperations.chromatin_trace_table import ChromatinTraceTable
from matrixOperations.HIMmatrixOperations import (
    calculate_contact_probability_matrix,
    plot_distance_histograms,
    plot_matrix,
)

warnings.filterwarnings("ignore")

# =============================================================================
# CLASSES
# =============================================================================


class BuildMatrix:
    def __init__(self, param):

        self.current_param = param

        self.initialize_parameters()

        # initialize with default values
        self.current_folder = []

    def initialize_parameters(self):
        # initializes parameters from current_param

        self.tracing_method = get_dictionary_value(
            self.current_param.param_dict["buildsPWDmatrix"],
            "tracing_method",
            default="masking",
        )
        self.z_binning = get_dictionary_value(
            self.current_param.param_dict["acquisition"], "zBinning", default=1
        )
        self.pixel_size_xy = get_dictionary_value(
            self.current_param.param_dict["acquisition"], "pixelSizeXY", default=0.1
        )
        self.pixel_size_z_0 = get_dictionary_value(
            self.current_param.param_dict["acquisition"], "pixelSizeZ", default=0.25
        )
        self.pixel_size_z = self.z_binning * self.pixel_size_z_0
        self.pixel_size = [self.pixel_size_xy, self.pixel_size_xy, self.pixel_size_z]
        self.available_masks = get_dictionary_value(
            self.current_param.param_dict["buildsPWDmatrix"],
            "masks2process",
            default={"nuclei": "DAPI"},
        )
        self.log_name_md = self.current_param.param_dict["fileNameMD"]
        self.mask_expansion = get_dictionary_value(
            self.current_param.param_dict["buildsPWDmatrix"],
            "mask_expansion",
            default=8,
        )
        self.available_masks = self.current_param.param_dict["buildsPWDmatrix"][
            "masks2process"
        ]
        self.colormaps = self.current_param.param_dict["buildsPWDmatrix"]["colormaps"]

    def calculate_pwd_single_mask(self, x, y, z):
        """
        Calculates PWD between barcodes detected in a given mask. For this:
            - converts xyz pixel coordinates into nm using self.pixel_size dictionary
            - calculates pair-wise distance matrix in nm
            - converts it into pixel units using self.pixel_size['x'] as an isotropic pixelsize.

        Parameters
        ----------
        r1: list of floats with xyz coordinates for spot 1 in microns
        r2: list of floats with xyz coordinates for spot 2 in microns

        Returns
        -------
        Returns pairwise distance matrix between barcodes in microns

        """

        # x = np.array([r1[0], r2[0]])
        # y = np.array([r1[1], r2[1]])
        # z = np.array([r1[2], r2[2]])

        # r_mum = np.column_stack((x, y, z))
        r_mum = np.column_stack((x, y, z))

        pairwd = pairwise_distances(r_mum)

        return pairwd

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
        self.unique_barcodes list of unique barcodes

        """
        # detects number of unique traces from trace table
        number_matrices = len(unique(self.trace_table.data, keys="Trace_ID"))

        # finds unique barcodes from trace table
        unique_barcodes = unique(self.trace_table.data, keys="Barcode #")[
            "Barcode #"
        ].data
        number_unique_barcodes = unique_barcodes.shape[0]

        print_log(
            f"$ Found {number_unique_barcodes} barcodes and {number_matrices} traces.",
            "INFO",
        )

        # Initializes sc_matrix
        sc_matrix = np.zeros(
            (number_unique_barcodes, number_unique_barcodes, number_matrices)
        )
        sc_matrix[:] = np.NaN

        # loops over traces
        print_log("> Processing traces...", "INFO")
        data_traces = self.trace_table.data.group_by("Trace_ID")
        for trace, trace_id, itrace in tzip(
            data_traces.groups, data_traces.groups.keys, range(number_matrices)
        ):

            barcodes_to_process = trace["Barcode #"].data

            # gets lists of x, y and z coordinates for barcodes assigned to a cell mask
            x, y, z = (
                np.array(trace["x"].data),
                np.array(trace["y"].data),
                np.array(trace["z"].data),
            )
            pwd_matrix = self.calculate_pwd_single_mask(x, y, z)

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
                        newdistance = pwd_matrix[ibarcode1, ibarcode2]

                        # inserts newdistance into sc_matrix using desired method
                        if mode == "last":
                            sc_matrix[index_barcode_1][index_barcode_2][itrace] = newdistance
                        elif mode == "mean":
                            sc_matrix[index_barcode_1][index_barcode_2][itrace] = np.nanmean(
                                [
                                    newdistance,
                                    sc_matrix[index_barcode_1][index_barcode_2][itrace],
                                ]
                            )
                        elif mode == "min":
                            sc_matrix[index_barcode_1][index_barcode_2][itrace] = np.nanmin(
                                [
                                    newdistance,
                                    sc_matrix[index_barcode_1][index_barcode_2][itrace],
                                ]
                            )

        self.sc_matrix = sc_matrix
        self.unique_barcodes = unique_barcodes

    def calculate_n_matrix(self):

        number_cells = self.sc_matrix.shape[2]

        if number_cells > 0:
            n_matrix = np.sum(~np.isnan(self.sc_matrix), axis=2)
        else:
            number_barcodes = self.sc_matrix.shape[0]
            n_matrix = np.zeros((number_barcodes, number_barcodes))

        self.n_matrix = n_matrix

    def plots_all_matrices(self, file):
        """
        Plots all matrices after analysis

        Parameters
        ----------
        file : str
            trace file name used for get output filenames.

        Returns
        -------
        None.

        """
        number_rois = 1  # by default we plot one ROI at a time.
        output_filename = file.split(".")[0] + "_Matrix"
        clim_scale = 1.0  # factor to multiply the clim by. If 1, the clim will be the mean of the PWD distribution of the whole map
        pixel_size = 1  # this is 1 as coordinates are in microns.
        n_cells = self.sc_matrix.shape[2]

        # plots PWD matrix
        # uses KDE
        plot_matrix(
            self.sc_matrix,
            self.unique_barcodes,
            pixel_size,
            number_rois,
            output_filename,
            self.log_name_md,
            figtitle="PWD matrix - KDE",
            mode="KDE",  # median or KDE
            clim=clim_scale * np.nanmean(self.sc_matrix),
            n_cells=n_cells,
            c_m=self.colormaps["PWD_KDE"],
            cmtitle="distance, um",
            filename_ending="_PWDmatrixKDE.png",
        )

        # uses median
        plot_matrix(
            self.sc_matrix,
            self.unique_barcodes,
            pixel_size,
            number_rois,
            output_filename,
            self.log_name_md,
            figtitle="PWD matrix - median",
            mode="median",  # median or KDE
            clim=clim_scale * np.nanmean(self.sc_matrix),
            cmtitle="distance, um",
            n_cells=n_cells,
            c_m=self.colormaps["PWD_median"],
            filename_ending="_PWDmatrixMedian.png",
        )

        # calculates and plots contact probability matrix from merged samples/datasets
        him_matrix, n_cells = calculate_contact_probability_matrix(
            self.sc_matrix, self.unique_barcodes, pixel_size, norm="nonNANs",
        )  # norm: n_cells (default), nonNANs

        c_scale = him_matrix.max()
        plot_matrix(
            him_matrix,
            self.unique_barcodes,
            pixel_size,
            number_rois,
            output_filename,
            self.log_name_md,
            figtitle="Hi-M matrix",
            mode="counts",
            clim=c_scale,
            n_cells=n_cells,
            c_m=self.colormaps["contact"],
            cmtitle="proximity frequency",
            filename_ending="_HiMmatrix.png",
        )

        # plots n_matrix
        plot_matrix(
            self.n_matrix,
            self.unique_barcodes,
            pixel_size,
            number_rois,
            output_filename,
            self.log_name_md,
            figtitle="N-matrix",
            mode="counts",
            n_cells=n_cells,
            clim=np.max(self.n_matrix),
            c_m=self.colormaps["Nmatrix"],
            cmtitle="number of measurements",
            filename_ending="_Nmatrix.png",
        )

        plot_distance_histograms(
            self.sc_matrix,
            pixel_size,
            output_filename,
            self.log_name_md,
            mode="KDE",
            kernel_width=0.25,
            optimize_kernel_width=False,
        )

    def save_matrices(self, file):

        output_filename = file.split(".")[0] + "_Matrix"

        # saves output
        np.save(output_filename + "_PWDscMatrix.npy", self.sc_matrix)
        np.savetxt(
            output_filename + "_uniqueBarcodes.ecsv",
            self.unique_barcodes,
            delimiter=" ",
            fmt="%d",
        )
        np.save(output_filename + "_Nmatrix.npy", self.n_matrix)

    def launch_analysis(self, file):
        """
        run analysis for a chromatin trace table.

        Returns
        -------
        None.

        """

        # creates and loads trace table
        self.trace_table = ChromatinTraceTable()
        self.trace_table.load(file)

        # runs calculation of PWD matrix
        self.build_distance_matrix("min")  # mean min last

        # calculates N-matrix: number of PWD distances for each barcode combination
        self.calculate_n_matrix()

        # runs plotting operations
        self.plots_all_matrices(file)

        # saves matrix
        self.save_matrices(file)

    def run(self):

        # initializes session_name, data_folder, current_folder
        self.label = "barcode"
        self.data_folder, self.current_folder = initialize_module(
            self.current_param, module_name="build_matrix", label=self.label
        )

        # reads chromatin traces
        files = [
            x
            for x in glob.glob(
                self.data_folder.output_folders["buildsPWDmatrix"]
                + os.sep
                + "Trace_*.ecsv"
            )
            if "uniqueBarcodes" not in x
        ]

        if len(files) < 1:
            print_log("$ No chromatin trace table found !", "WARN")
            return

        print_log(f"> Will process {len(files)} trace tables with names:")
        for file in files:
            print_log(f"{os.path.basename(file)}")

        for file in files:
            self.launch_analysis(file)

        print_log(
            f"$ {len(files)} chromatin trace tables processed in {self.current_folder}"
        )
