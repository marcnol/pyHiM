#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 14:11:58 2022

@author: marcnol

This script will build chromatin traces using a segmentObjects_barcode table

The methods that will be implemented are:
    1= assigment by mask (either DAPI mask or other)
    2= spatial clusterization using KDtree. This method is mask-free.



Method 1:
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
from apifish.stack.io import read_array
from astropy.table import Table
from scipy.spatial import KDTree
from skimage.segmentation import expand_labels
from sklearn.metrics import pairwise_distances
from tqdm import trange
from tqdm.contrib import tzip

from fileProcessing.fileManagement import (
    Folders,
    get_dictionary_value,
    print_log,
    write_string_to_file,
)
from imageProcessing.localization_table import LocalizationTable
from matrixOperations.chromatin_trace_table import ChromatinTraceTable
from matrixOperations.filter_localizations import get_file_table_new_name
from matrixOperations.HIMmatrixOperations import (
    calculate_contact_probability_matrix,
    plot_distance_histograms,
    plot_matrix,
)

warnings.filterwarnings("ignore")

# =============================================================================
# CLASSES
# =============================================================================


class BuildTraces:
    def __init__(self, param):
        self.current_param = param

        self.initialize_parameters()

        # initialize with default values
        self.current_folder = []
        self.mask_identifier = ["DAPI"]  # default mask label

    def initializes_masks(self, masks):
        self.masks = masks
        self.n_cells_assigned = 0
        self.n_cells_unassigned = 0
        self.n_barcodes_in_mask = 0
        self.number_masks = np.max(self.masks).astype(int)
        self.barcodes_in_mask = {}

        for mask in range(self.number_masks + 1):
            self.barcodes_in_mask["maskID_" + str(mask)] = []

    def initialize_parameters(self):
        # initializes parameters from current_param

        self.tracing_method = get_dictionary_value(
            self.current_param.param_dict["buildsPWDmatrix"],
            "tracing_method",
            default=["masking"],
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
        self.kd_tree_distance_threshold_mum = get_dictionary_value(
            self.current_param.param_dict["buildsPWDmatrix"],
            "kd_tree_distance_threshold_mum",
            default=1,
        )

    def initialize_lists(self):
        (
            self.rois,
            self.cell_id,
            self.n_barcodes,
            self.barcode_ids,
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
        )

    def align_by_masking(self):
        """
        Assigns barcodes to masks and creates <n_barcodes_in_mask>
        And by filling in the "Cell #" key of barcode_map_roi
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
        n_barcodes_roi = 0

        image_size = self.masks.shape

        # loops over barcode Table rows in a given ROI
        print_log(f"> Aligning localizations to {self.number_masks} masks...")
        for i in trange(
            len(self.barcode_map_roi.groups[0])
        ):  # i is the index of the barcode in barcode_map_roi
            barcode = self.barcode_map_roi.groups[0]["Barcode #"][i]

            # gets xyz coordinates
            x_corrected = self.barcode_map_roi.groups[0]["ycentroid"][i]
            y_corrected = self.barcode_map_roi.groups[0]["xcentroid"][i]

            if self.ndims == 2:
                z_corrected = self.barcode_map_roi.groups[0]["zcentroid"][i] = 0.0
            else:
                z_corrected = self.barcode_map_roi.groups[0]["zcentroid"][i]

            # binarizes coordinate
            y_int = binarize_coordinate(y_corrected)
            x_int = binarize_coordinate(x_corrected)

            # finds what mask label this barcode is sitting on
            if np.isnan(x_int) or np.isnan(y_int):
                # if a barcode has coordinates with NaNs, it is assigned to background
                mask_id = 0
            else:
                if (
                    x_int < image_size[0]
                    and y_int < image_size[1]
                    and x_int > 0
                    and y_int > 0
                ):
                    mask_id = self.masks[x_int][y_int]
                else:
                    # if a barcode has coordinates outside the image, it is assigned to background
                    mask_id = 0

            # attributes CellID to a barcode
            self.barcode_map_roi["CellID #"][i] = mask_id

            # if it is not background,
            if mask_id > 0:
                # increments counter of number of barcodes in the cell mask attributed
                n_barcodes_in_mask[mask_id] += 1

                # stores the identify of the barcode in a mask dictionary
                self.barcodes_in_mask["maskID_" + str(mask_id)].append(i)

                # keeps statistics
                if int(self.barcode_map_roi.groups[0]["ROI #"][i]) == int(self.n_roi):
                    n_barcodes_roi += 1

        # Total number of masks assigned and not assigned
        self.n_cells_assigned = np.count_nonzero(n_barcodes_in_mask > 0)
        self.n_cells_unassigned = self.number_masks - self.n_cells_assigned

        # this list contains which barcodes are allocated to which masks
        self.n_barcodes_in_mask = n_barcodes_in_mask

        print_log(
            "$ Number of cells assigned: {} | discarded: {}".format(
                self.n_cells_assigned, self.n_cells_unassigned
            )
        )

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

    def builds_sc_distance_table(self):
        """
        iterates over all masks, calculates PWD for each mask, assigns them to sc_distance_table

        Returns
        -------
        sc_distance_table

        """
        # sorts Table by cellID
        barcode_map_roi = self.barcode_map_roi

        # indexes table by cellID
        barcode_map_roi_cell_id = barcode_map_roi.group_by("CellID #")

        self.initialize_lists()

        # iterates over all traces in an ROI
        print_log("> Building single traces")
        for key, group in tzip(
            barcode_map_roi_cell_id.groups.keys, barcode_map_roi_cell_id.groups
        ):
            if key["CellID #"] > 1:  # excludes trace 0 as this is background

                group_keys, cell_id, roi = (
                    group.keys(),
                    key["CellID #"],
                    group["ROI #"].data[0],
                )
                trace_id = str(uuid.uuid4())

                # gets lists of x, y and z coordinates for barcodes assigned to a cell mask
                x, y, z = (
                    np.array(group["xcentroid"].data),
                    np.array(group["ycentroid"].data),
                    np.array(group["zcentroid"].data),
                )

                # gets vector in nm
                r_mum = self.build_vector(group_keys, x, y, z)

                for i in range(x.shape[0]):
                    entry = [
                        group["Buid"].data[i],  # spot uid
                        trace_id,  # trace uid
                        r_mum[i][0],  # x, microns
                        r_mum[i][1],  # y, microns
                        r_mum[i][2],  # z, microns
                        "xxxxx",  # chromosome
                        0,  # start sequence
                        999999999,  # end sequence
                        roi,  # ROI number
                        cell_id,  # Mask number
                        group["Barcode #"].data[i],  # Barcode name
                        "x" * 20,  # label
                    ]
                    self.trace_table.data.add_row(entry)

        print_log("$ Coordinates dimensions: {}".format(self.ndims))

    def load_mask(self,files_in_folder,):
        """
        searches and loads mask files for building chromatin trace

        Parameters
        ----------
        files_in_folder : list of str
            list of TIF files to be explored.

        Returns
        -------
        bool
            True: mask found and loaded
            False: failed to find mask file

        """
        
        # finds files with cell masks
        channel = self.current_param.param_dict["acquisition"][
            self.mask_type + "_channel"
        ]

        files_to_process = [
            file
            for file in files_in_folder
            if self.current_param.decode_file_parts(file)["channel"]
            == channel  # typically "ch00"
            and self.mask_identifier in os.path.basename(file).split("_")
            and int(self.current_param.decode_file_parts(file)["roi"]) == self.n_roi
        ]

        if len(files_to_process) > 0:

            # loads file with cell masks
            filename_roi_masks = (
                os.path.basename(files_to_process[0]).split(".")[0] + "_Masks.npy"
            )
            full_filename_roi_masks = (
                self.data_folder.output_folders["segmentedObjects"]
                + os.sep
                + filename_roi_masks
            )

            if os.path.exists(full_filename_roi_masks):

                # loads and initializes masks
                segmented_masks = read_array(full_filename_roi_masks)
                print(f"$ loaded mask file: {full_filename_roi_masks}")
                
                # expands mask without overlap by a maximmum of 'distance' pixels
                self.masks = expand_labels(
                    segmented_masks, distance=self.mask_expansion
                )

                # initializes masks
                self.initializes_masks(self.masks)
                return True

            else:
                # Could not find a file with masks to assign. Report and continue with next ROI
                debug_mask_filename(
                    files_in_folder,
                    full_filename_roi_masks,
                    self.mask_identifier,
                    self.n_roi,
                    label=self.current_param.param_dict["acquisition"]["label_channel"],
                )

        else:
            print_log(f"$ Did not find any filename for mask: {self.mask_identifier}, channel: {channel}","WARN")
            print_log("-"*80)
            # Could not find a file with masks to assign. Report and continue with next ROI
            debug_mask_filename(files_in_folder,"None",self.mask_identifier,self.n_roi,label=self.current_param.param_dict["acquisition"]["label_channel"])

        return False

    def assign_masks(
        self, output_filename, barcode_map,
    ):
        """
        Main function that:
            loads and processes barcode localization files, local alignment file, and masks
            initializes <cell_roi> class and assigns barcode localizations to masks
            then constructs the single cell PWD matrix and outputs it toghether with the contact map and the N-map.

        Parameters
        ----------
        output_filename : string
        self.current_param : Parameters Class
        self.current_folder : string
        self.data_folder : Folder Class
            information to find barcode localizations, local drift corrections and masks

        self.pixel_size : dict, optional
            pixel_size = {'x': pixelSizeXY,
                        'y': pixelSizeXY,
                        'z': pixel_size_z}
            The default is 0.1 for x and y, 0.0 for z. Pixelsize in um

        self.log_name_md : str, optional
            Filename of Markdown output. The default is "log.md".
        self.ndims : int, optional
            indicates whether barcodes were localized in 2 or 3D. The default is 2.
        self.mask_identifier:

        Returns
        -------
        None.

        """

        # indexes localization tables by ROI
        barcode_map_roi = barcode_map.group_by("ROI #")
        number_rois = len(barcode_map_roi.groups.keys)

        print_log("-" * 80)
        print_log(
            " Loading masks and pre-processing barcodes for Mask <{}> for {} rois".format(
                self.mask_identifier, number_rois
            )
        )

        # finds TIFs in current_folder
        tif_files_in_folder = glob.glob(self.current_folder + os.sep + "*.tif")

        # loops over rois
        processing_order = 0

        for roi in range(number_rois):
            self.n_roi = barcode_map_roi.groups.keys[roi][
                0
            ]  # need to iterate over the first index
            self.barcode_map_roi = barcode_map.group_by("ROI #").groups[roi]

            mask_loaded = self.load_mask(tif_files_in_folder,)
            if mask_loaded:

                print_log("> Processing ROI# {}".format(self.n_roi))

                # initializes trace table
                self.trace_table.initialize()

                # finds what barcodes are in each cell mask
                self.align_by_masking()
                print_log(f"$ ROI: {roi}, N cells assigned: {self.n_cells_assigned - 1} out of {self.number_masks}\n")

                # builds sc_distance_table
                self.builds_sc_distance_table()
                print_log(
                    "$ Number of entries in trace table: {}".format(
                        len(self.trace_table.data)
                    )
                )

                if len(self.trace_table.data) > 0:
                    # saves trace table with results per ROI
                    output_table_filename = (
                        output_filename
                        + "_"
                        + self.label
                        + "_mask:"
                        + str(self.mask_identifier)
                        + "_ROI:"
                        + str(self.n_roi)
                        + ".ecsv"
                    )
                    self.trace_table.save(output_table_filename, self.trace_table.data)

                    # plots results
                    self.trace_table.plots_traces(
                        [output_table_filename.split(".")[0], "_traces_XYZ", ".png"],
                        masks=self.masks,
                    )

                    print_log(f"$ Saved output table as {output_table_filename}")
                else:
                    print_log(f"! Warning: table was empty therefore not saved!")

                processing_order += 1

    def build_trace_by_masking(self, barcode_map):

        print_log("> Masks labels: {}".format(self.available_masks))

        for mask_label in self.available_masks.keys():

            self.mask_identifier = self.available_masks[mask_label]
            if "DAPI" in self.mask_identifier:
                self.mask_type = "DAPI"
            else:
                self.mask_type = "mask"

            tag = str(self.ndims) + "D"

            output_filename = (
                self.data_folder.output_folders["buildsPWDmatrix"]
                + os.sep
                + "Trace_"
                + tag
            )

            # creates and initializes trace table
            self.trace_table = ChromatinTraceTable()

            self.assign_masks(
                output_filename, barcode_map,
            )

            print_log(
                "$ Trace built using mask assignment. Output saved in: {} ".format(
                    self.current_folder
                ),
                "info",
            )

    def group_localizations_by_coordinate(self,):
        """
        Uses a KDTree to group detections by it's coordinates, given a certain distance threshold
        Returns a list of lists. Each list contains the lines if the input data (segmentedObjects_3D_barcode.dat)
        where the detections are less than a pixel away from each other

        Parameters
        ----------
        coordinates : numpy array, float
            Matrix containing the xyz coordinates of barcodes.
        distance_threshold : float, defaul 1.0
            Distance threshold in pixels used to detect neighboring barcodes.

        Returns
        -------
        group_list : list
            list of lists containing the coordinates of barcodes associated together.
        """

        # gets coordinates from trace table
        data_table = self.barcode_map_roi
        pixel_size = self.pixel_size
        len_data_table = len(data_table)

        if self.ndims == 3:
            coordinates = np.concatenate(
                [
                    pixel_size["x"] * data_table["xcentroid"].data.reshape(len_data_table, 1),
                    pixel_size["y"] * data_table["ycentroid"].data.reshape(len_data_table, 1),
                    pixel_size["z"] * data_table["zcentroid"].data.reshape(len_data_table, 1),
                ],
                axis=1,
            )
        elif self.ndims == 2:
            coordinates = np.concatenate(
                                    [
                                        pixel_size['x']*data_table['xcentroid'].data.reshape(len_data_table,1),
                                        pixel_size['y']*data_table['ycentroid'].data.reshape(len_data_table,1),
                                        0.0*data_table['zcentroid'].data.reshape(len_data_table,1),
                                    ],
                                    axis = 1,
                                )
        """ if this code above works and does not introduce bugs, we will remove the commented lines in future
        elif self.ndims == 2:
            coordinates = np.concatenate([pixel_size['x']*data_table['xcentroid'].data.reshape(len_data_table,1),
                                    pixel_size['y']*data_table['ycentroid'].data.reshape(len_data_table,1)], axis = 1)
        """
            
        # gets tree of coordinates
        print_log(f'> Creating KDTree for {self.ndims} dimensions')
        x_tree = KDTree(coordinates)

        ## set distance thresold
        r = self.kd_tree_distance_threshold_mum

        # Groups points when they're less than r away
        points = []
        for element in coordinates:
            points.append(x_tree.query_ball_point(element, r, p=2.0))

        # Get unique groups
        groups = list(set(tuple(x) for x in points))
        group_list = [list(x) for x in groups]
        self.n_cells_assigned = len(group_list)

        # Fills in trace information in localization table

        # iterates over traces
        for trace_id, trace in enumerate(group_list):

            # iterates over localizations in trace
            for i in range(len(trace)):
                # gets index of localization in data_table
                index_localization = trace[i]

                # records trace to which this index belongs
                data_table["CellID #"][index_localization] = trace_id

        self.barcode_map_roi = data_table

    def build_trace_by_clustering(
        self, barcode_map,
    ):

        # decompose by ROI!
        # indexes localization tables by ROI
        barcode_map_roi = barcode_map.group_by("ROI #")
        number_rois = len(barcode_map_roi.groups.keys)

        print_log("-" * 80)
        print_log("> Starting spatial clustering for {} ROI in {} dimensions".format(number_rois,self.ndims))

        tag = str(self.ndims) + "D"

        output_filename = (
            self.data_folder.output_folders["buildsPWDmatrix"] + os.sep + "Trace_" + tag
        )

        # creates and initializes trace table
        self.trace_table = ChromatinTraceTable()

        # loops over rois
        for roi in range(number_rois):
            self.n_roi = barcode_map_roi.groups.keys[roi][
                0
            ]  # need to iterate over the first index
            self.barcode_map_roi = barcode_map.group_by("ROI #").groups[roi]

            print_log("$ Processing ROI# {}".format(self.n_roi))

            # initializes trace table
            self.trace_table.initialize()

            # build traces by spatial clustering
            self.group_localizations_by_coordinate()
            print_log(f"$ ROI: {roi}, N cells assigned: {self.n_cells_assigned - 1}\n")

            # builds sc_distance_table
            self.builds_sc_distance_table()
            if len(self.trace_table.data) > 0:
                print_log(
                    "$ Number of entries in trace table: {}".format(
                        len(self.trace_table.data)
                    )
                )

                # saves trace table with results per ROI
                output_table_filename = (
                    output_filename
                    + "_"
                    + self.label
                    + "_KDtree"
                    + "_ROI:"
                    + str(self.n_roi)
                    + ".ecsv"
                )
                self.trace_table.save(output_table_filename, self.trace_table.data)

                # plots results
                self.trace_table.plots_traces(
                    [output_table_filename.split(".")[0], "_traces_XYZ", ".png"]
                )

                print_log(f"$ Traces built. Saved output table as {output_table_filename}")
            else:
                print_log(f"! Warning: table was empty therefore not saved!")
                

    def launch_analysis(self, file):

        # loads barcode coordinate Tables
        table = LocalizationTable()
        barcode_map, self.unique_barcodes = table.load(file)

        if "3D" in os.path.basename(file):
            self.ndims = 3
            self.pixel_size = {
                "x": self.pixel_size_xy,
                "y": self.pixel_size_xy,
                "z": self.pixel_size_z,
            }
        else:
            self.ndims = 2
            self.pixel_size = {"x": self.pixel_size_xy, "y": self.pixel_size_xy, "z": 0}

        if (
            "clustering" in self.tracing_method and self.ndims == 3
        ):  # for now it only runs for 3D data
            self.build_trace_by_clustering(barcode_map)
        elif self.ndims == 2:
            print_log(f"! Warning: localization files in 2D will not be processed using clustering.\n")
            
        if "masking" in self.tracing_method:
            self.build_trace_by_masking(barcode_map)

    def run(self):
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
        # initializes session_name, data_folder, current_folder
        self.label = "barcode"
        self.data_folder, self.current_folder = initialize_module(
            self.current_param, module_name="build_traces", label=self.label
        )

        print_log("> Masks labels: {}".format(self.available_masks))

        # iterates over barcode localization tables in the current folder
        files = [
            x
            for x in glob.glob(
                self.data_folder.output_files["segmentedObjects"]
                + "_*"
                + self.label
                + ".dat"
            )
        ]

        if len(files) < 1:
            print_log("$ No localization table found to process!", "WARN")
            return

        print_log(f"> Will process {len(files)} localization tables with names:")
        for file in files:
            print_log(f"{os.path.basename(file)}")

        for file in files:
            self.launch_analysis(file)

        print_log(f"$ {len(files)} barcode tables processed in {self.current_folder}")


def initialize_module(current_param, module_name="build_traces", label="barcode"):

    session_name = module_name

    # processes folders and files
    data_folder = Folders(current_param.param_dict["rootFolder"])
    print_log("\n" + "=" * 35 + f"{session_name}" + "=" * 35 + "\n")
    print_log("$ folders read: {}".format(len(data_folder.list_folders)))
    write_string_to_file(
        current_param.param_dict["fileNameMD"],
        "## {}\n".format(session_name),
        "a",
    )

    current_folder = current_param.param_dict["rootFolder"]
    data_folder.create_folders(current_folder, current_param)
    print_log("> Processing Folder: {}".format(current_folder))

    return data_folder, current_folder


def debug_mask_filename(
    files_in_folder, full_filename_roi_masks, mask_identifier, n_roi, label=""
):

    print_log(f"# Error, no mask file found for ROI: {n_roi}\n")
    print_log("# File I was searching for: {}".format(full_filename_roi_masks))
    print_log("# Debug: ")
    for file in files_in_folder:
        if (
            file.split("_")[-1].split(".")[0] == label  # typically "ch00"
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


def binarize_coordinate(x):
    if not np.isnan(x):
        return int(x)
    else:
        return np.nan
