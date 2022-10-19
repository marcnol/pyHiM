#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 21:21:02 2020

@author: marcnol

contains functions and classes needed for the analysis and plotting of HiM matrices

"""

# =============================================================================
# IMPORTS
# =============================================================================


import csv
import glob
import json
import os

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack
from numba import jit
from pylab import colorbar, contourf
from scipy import interpolate
from scipy.io import loadmat
from sklearn import manifold
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.neighbors import KernelDensity
from tqdm import trange

from fileProcessing.fileManagement import is_notebook, write_string_to_file

# =============================================================================
# CLASSES
# =============================================================================


class AnalysisHiMMatrix:
    """
    this class is used for loading data processed by processHiMmatrix.py
    Main use is to produce paper quality figures of HiM matrices, 3-way interaction matrices and HiM matrix ratios
    """

    def __init__(self, run_parameters, root_folder="."):
        self.data_folder = root_folder + os.sep + "scHiMmatrices"
        self.run_parameters = run_parameters
        self.root_folder = root_folder
        self.data = []
        self.data_files = []
        self.folders_to_load = []
        self.number_barcodes = 0

    def load_data(self):
        """
        loads dataset

        Returns
        -------
        self.foldes2Load contains the parameters used for the processing of HiM matrices.
        self.data_files dictionary containing the extensions needed to load data files
        self.data dictionary containing the datasets loaded
        """

        # loads datasets: parameter files
        filename_list_data_json = (
            self.root_folder + os.sep + self.run_parameters["parametersFileName"]
        )
        with open(filename_list_data_json, encoding="utf-8") as json_file:
            list_data = json.load(json_file)

        dataset_name = list(list_data.keys())[0]
        print("Dataset: {}".format(dataset_name))

        output_filename = (
            self.data_folder
            + os.sep
            + dataset_name
            + "_label:"
            + self.run_parameters["label"]
            + "_action:"
            + self.run_parameters["action"]
        )

        filename_parameters_json = output_filename + "_parameters.json"
        with open(filename_parameters_json, encoding="utf-8") as json_file:
            folders_to_load = json.load(json_file)
        print("Loading parameter file:".format(filename_parameters_json))

        # Creates filenames to be loaded
        data_files = {}
        data_files["ensembleContactProbability"] = "_ensembleContactProbability.npy"
        data_files["SCmatrixCollated"] = "_SCmatrixCollated.npy"
        data_files["SClabeledCollated"] = "_SClabeledCollated.npy"

        if "3wayContacts_anchors" in list_data[dataset_name]:
            for i_anchor in list_data[dataset_name]["3wayContacts_anchors"]:
                new_key = "anchor:" + str(i_anchor - 1)
                data_files[new_key] = "_" + new_key + "_ensemble3wayContacts.npy"
        else:
            print("No anchors found")

        # loads datasets: numpy matrices
        data = {}
        print("Loading datasets from: {}".format(output_filename))
        for i_data_file in data_files.keys():
            print(
                "Loaded: {}: <{}>".format(
                    i_data_file, os.path.basename(output_filename + data_files[i_data_file])
                )
            )
            data[i_data_file] = np.load(output_filename + data_files[i_data_file]).squeeze()

        # loads datasets: lists
        run_name = load_list(output_filename + "_runName.csv")
        data["runName"] = run_name
        print("Loaded runNames: {}".format(data["runName"]))

        data["uniqueBarcodes"] = load_list(output_filename + "_uniqueBarcodes.csv")
        print("Loaded barcodes #: {}".format(data["uniqueBarcodes"]))
        self.number_barcodes = len(data["uniqueBarcodes"])

        print(
            "Total number of cells loaded: {}".format(data["SCmatrixCollated"].shape[2])
        )
        print("Number Datasets loaded: {}".format(len(data["runName"])))

        # Exports data
        self.data = data
        self.data_files = data_files
        self.folders_to_load = folders_to_load
        self.list_data = list_data
        self.dataset_name = dataset_name

    # functions

    def plot_2d_matrix_simple(
        self,
        ifigure,
        matrix,
        unique_barcodes,
        yticks,
        xticks,
        cmtitle="probability",
        c_min=0,
        c_max=1,
        c_m="coolwarm",
        fontsize=12,
        colorbar=False,
        axis_ticks=False,
        n_cells=0,
        n_datasets=0,
        show_title=False,
        fig_title="",
    ):

        pos = ifigure.imshow(matrix, cmap=c_m)  # colormaps RdBu seismic

        if show_title:
            title_text = "{} | N = {} | n = {}".format(fig_title, n_cells, n_datasets)
            ifigure.title.set_text(title_text)

        # plots figure
        if xticks:
            ifigure.set_xlabel("barcode #", fontsize=fontsize)
            if not axis_ticks:
                ifigure.set_xticklabels(())
            else:
                print("barcodes:{}".format(unique_barcodes))
                # ifigure.set_xticks(np.arange(matrix.shape[0]),unique_barcodes)
                ifigure.set_xticklabels(unique_barcodes)

        else:
            ifigure.set_xticklabels(())
        if yticks:
            ifigure.set_ylabel("barcode #", fontsize=fontsize)
            if not axis_ticks:
                ifigure.set_yticklabels(())
            else:
                # ifigure.set_yticks(np.arange(matrix.shape[0]), unique_barcodes)
                ifigure.set_yticklabels(unique_barcodes)
        else:
            ifigure.set_yticklabels(())

        for xtick, ytick in zip(
            ifigure.xaxis.get_majorticklabels(), ifigure.yaxis.get_majorticklabels()
        ):
            xtick.set_fontsize(fontsize)
            ytick.set_fontsize(fontsize)

        if colorbar:
            cbar = plt.colorbar(pos, ax=ifigure, fraction=0.046, pad=0.04)
            cbar.minorticks_on()
            cbar.set_label(cmtitle, fontsize=float(fontsize) * 0.85)
            pos.set_clim(vmin=c_min, vmax=c_max)

        pos.set_clim(vmin=c_min, vmax=c_max)

        return pos

    def update_clims(self, c_min, c_max, axes):
        for ax in axes:
            ax.set_clim(vmin=c_min, vmax=c_max)

    def plot_1d_profile1dataset(self, ifigure, anchor, i_fig_label, yticks, xticks):

        prop_cycle = plt.rcParams["axes.prop_cycle"]
        colors = prop_cycle.by_key()["color"]
        lwbase = plt.rcParams["lines.linewidth"]
        thin, thick = lwbase / 2, lwbase * 3

        profile = self.data["ensembleContactProbability"][:, anchor - 1]
        x = np.linspace(0, profile.shape[0], num=profile.shape[0], endpoint=True)
        # f = interp1d(x, profile,kind = 'linear') # linear
        tck = interpolate.splrep(x, profile, s=0)
        xnew = np.linspace(0, profile.shape[0], num=100, endpoint=True)
        ynew = interpolate.splev(xnew, tck, der=0)
        if self.run_parameters["splines"]:
            ifigure.plot(xnew, ynew, "-")  # x, profile, 'o',
        else:
            ifigure.plot(x, profile, "-")  # x, profile, 'o',

        ifigure.set_xlim([0, profile.shape[0]])
        ifigure.axvline(x=anchor - 0.5, color=colors[4], lw=thick, alpha=0.5)
        ifigure.set_ylim([0, self.run_parameters["cAxis"]])

        if xticks:
            ifigure.set_xlabel("barcode #", fontsize=self.run_parameters["fontsize"])
            if not self.run_parameters["axisTicks"]:
                ifigure.set_xticklabels(())
            else:
                ifigure.set_xticklabels(self.data["uniqueBarcodes"])
        else:
            ifigure.set_xticklabels(())

        if yticks:
            ifigure.set_ylabel("Probability", fontsize=self.run_parameters["fontsize"])
            if not self.run_parameters["axisTicks"]:
                ifigure.set_yticklabels(())
            else:
                ifigure.set_yticks(
                    [0, self.run_parameters["cAxis"] / 2, self.run_parameters["cAxis"]]
                )
        else:
            ifigure.set_yticklabels(())

    def n_cells_loaded(self):
        if self.run_parameters["action"] == "labeled":
            cells_with_label = [
                idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x > 0
            ]
            n_cells = len(cells_with_label)
        elif self.run_parameters["action"] == "unlabeled":
            cells_with_label = [
                idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x == 0
            ]
            n_cells = len(cells_with_label)
        else:
            n_cells = self.data["SCmatrixCollated"].shape[2]
        print("n_cells selected with label: {}".format(n_cells))
        return n_cells

    def retrieve_sc_matrix(self):
        """
        retrieves single cells that have the label requested

        Returns
        -------
        self.sc_matrix_selected

        """
        n_cells = self.n_cells_loaded()
        sc_matrix_selected = np.zeros((self.number_barcodes, self.number_barcodes, n_cells))

        if self.run_parameters["action"] == "labeled":
            cells_with_label = [
                idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x > 0
            ]
            new_cell = 0
            for i_cell in cells_with_label:
                sc_matrix_selected[:, :, new_cell] = self.data["SCmatrixCollated"][
                    :, :, i_cell
                ]
                new_cell += 1
        elif self.run_parameters["action"] == "unlabeled":
            cells_with_label = [
                idx for idx, x in enumerate(self.data["SClabeledCollated"]) if x == 0
            ]
            new_cell = 0
            for i_cell in cells_with_label:
                sc_matrix_selected[:, :, new_cell] = self.data["SCmatrixCollated"][
                    :, :, i_cell
                ]
                new_cell += 1
        else:
            sc_matrix_selected = self.data["SCmatrixCollated"]
        print("n_cells retrieved: {}".format(sc_matrix_selected.shape[2]))
        self.sc_matrix_selected = sc_matrix_selected


# =============================================================================
# FUNCTIONS
# =============================================================================


def normalize_profile(profile1, profile2, run_parameters):

    print("Normalization: {}".format(run_parameters["normalize"]))

    mode = run_parameters["normalize"]

    if "maximum" in mode:  # normalizes by maximum
        profile1 = profile1 / profile1.max() / 2
        profile2 = profile2 / profile2.max() / 2
    elif "none" in mode:  # no normalization
        m1_norm = 1
        m2_norm = 1
    else:  # normalizes by given factor
        norm_factor = float(mode)
        profile1 = profile1 / 1
        profile2 = profile2 / norm_factor

    return profile1, profile2


def plot_1d_profile2datasets(
    ifigure,
    him_data_1,
    him_data_2,
    run_parameters,
    anchor,
    i_fig_label,
    yticks,
    xticks,
    legend=False,
):

    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]
    lwbase = plt.rcParams["lines.linewidth"]
    thin, thick = lwbase / 2, lwbase * 3

    profile1 = him_data_1.data["ensembleContactProbability"][:, anchor - 1]
    profile2 = him_data_2.data["ensembleContactProbability"][:, anchor - 1]

    profile1, profile2 = normalize_profile(profile1, profile2, run_parameters)

    x = np.linspace(0, profile1.shape[0], num=profile1.shape[0], endpoint=True)
    tck1 = interpolate.splrep(x, profile1, s=0)
    tck2 = interpolate.splrep(x, profile2, s=0)
    xnew = np.linspace(0, profile1.shape[0], num=100, endpoint=True)
    ynew1 = interpolate.splev(xnew, tck1, der=0)
    ynew2 = interpolate.splev(xnew, tck2, der=0)
    if run_parameters["splines"]:
        ifigure.plot(xnew, ynew1, "-", xnew, ynew2, "-")  # x, profile, 'o',
    else:
        ifigure.plot(x, profile1, "-", x, profile2, "-")  # x, profile, 'o',

    ifigure.set_xlim([0, profile1.shape[0]])
    ifigure.axvline(x=anchor - 0.5, color=colors[4], lw=thick, alpha=0.5)
    ifigure.set_ylim([0, run_parameters["cAxis"]])

    if xticks:
        ifigure.set_xlabel("barcode #", fontsize=run_parameters["fontsize"])
        if not run_parameters["axisTicks"]:
            ifigure.set_xticklabels(())
        else:
            ifigure.set_xticklabels(him_data_1.data["uniqueBarcodes"])
    else:
        ifigure.set_xticklabels(())

    if yticks:
        ifigure.set_ylabel("Probability", fontsize=run_parameters["fontsize"])
        if not run_parameters["axisTicks"]:
            ifigure.set_yticklabels(())
        else:
            ifigure.set_yticks(
                [0, run_parameters["cAxis"] / 2, run_parameters["cAxis"]]
            )
    else:
        ifigure.set_yticklabels(())

    if legend:
        ifigure.legend([him_data_1.dataset_name, him_data_2.dataset_name], loc="best")


def load_list(file_name):
    with open(file_name, newline="", encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ", quotechar="|")
        run_name = []
        for row in spamreader:
            # print(', '.join(row))
            if len(run_name) > 0:
                run_name.append(row)
            else:
                run_name = row

    return run_name


def attributes_labels2cells(snd_table, results_table, label="doc"):

    sorted_snd_table = snd_table.group_by("MaskID #")
    list_keys = list(sorted_snd_table.groups.keys["MaskID #"].data)
    index_key = [index for i, index in zip(list_keys, range(len(list_keys))) if i == label]

    # checks that there is at least one cell with the label
    if len(index_key) > 0:

        snd_table_with_label = sorted_snd_table.groups[index_key[0]]
        print("\n>>> Matching labels")
        print(
            "Found {} out of {} cells with {} in dataset".format(
                len(snd_table_with_label), len(sorted_snd_table), label
            )
        )

        # sorts Results Table by ROI
        pwd_table_sorted_roi = results_table.group_by("ROI #")
        cuids_list = []
        cuids = Table()
        cuids["Cuid"] = []

        print("rois to process: {}".format(pwd_table_sorted_roi.groups.keys))

        for roi, group in zip(pwd_table_sorted_roi.groups.keys, pwd_table_sorted_roi.groups):
            # list of cellIDs in ROI
            cells_to_process = group["CellID #"].data.compressed()
            cells_to_process_uid = group["Cuid"]

            # list of rois detected in snd_table_with_label
            rois_in_snd_table_with_label = list(
                snd_table_with_label.group_by("ROI #").groups.keys["ROI #"].data
            )

            # index of ROI within the keys of snd_table_with_label
            index_rois = [
                index
                for i, index in zip(
                    rois_in_snd_table_with_label, range(len(rois_in_snd_table_with_label))
                )
                if i == roi["ROI #"]
            ]

            # subtable of cells with label and ROI that we are looking for
            snd_table_with_label_roi = snd_table_with_label.group_by("ROI #").groups[
                index_rois[0]
            ]
            cells_with_label = list(snd_table_with_label_roi["CellID #"].data)

            # finds which cell indeces in Table have label
            list_of_selected_cells = [
                index
                for i_cell, index in zip(cells_to_process, range(len(cells_to_process)))
                if i_cell in cells_with_label
            ]

            if len(list_of_selected_cells) > 0:
                print(
                    "Detected {} cells in ROI {} with label".format(
                        len(list_of_selected_cells), roi["ROI #"]
                    )
                )
                if len(cuids) > 0:
                    # cuids = vstack([cuids, cells_to_process_uid[list_of_selected_cells]])
                    # print('adding {} more cells'.format(len(cells_to_process_uid[list_of_selected_cells])))
                    cuids_list += list(
                        cells_to_process_uid[list_of_selected_cells].data.compressed()
                    )
                else:
                    cuids_list = list(
                        cells_to_process_uid[list_of_selected_cells].data.compressed()
                    )
                    # cuids = cells_to_process_uid[list_of_selected_cells]

            print(
                "Processed ROI # {}, found {} out of {} cells with {}".format(
                    roi["ROI #"], len(list_of_selected_cells), len(group), label
                )
            )

        # from list of cuids from cells that show label, I construct a binary vector of the same size as sc_matrix. Labeled cells have a 1.
        sc_labeled = np.zeros(len(results_table))
        # print('CUID list: {}'.format(CUIDsList2))
        # cuids_list = cuids["Cuid"].data.compressed()
        # cuids_list = CUIDsList2
        # checks that there are cells found with the label
        if len(cuids_list) > 0:
            index_cells_with_label = [
                i_row
                for row, i_row in zip(results_table, range(len(results_table)))
                if row["Cuid"] in cuids_list
            ]
            sc_labeled[index_cells_with_label] = 1
        else:
            sc_labeled, cuids_list = [], []

        return sc_labeled, cuids_list
    else:
        # otherwise returns an empty list
        print("Warning: No cell with a mask labeled <{}> was found".format(label))
        return [], []


def load_sc_data(list_data, dataset_name, p):
    """
    loads SC datasets from a dict of folders (list_data)

    Parameters
    ----------
    list_data : dict
        dict of folders with data that can be loaded.
    dataset2Load : int, optional
        The item in the dictionary that will be loaded. The default is 3.

    Returns
    -------
    sc_matrix_collated : list of np arrays n_barcodes x n_barcodes x n_cells
        Cummulative SC PWD matrix.
    unique_barcodes : list of np arrays
        containing the barcode identities for each matrix.
    build_pwd_matrix_collated : list of Tables
        Tables with all the data for cells and barcodes used to produce sc_matrix_collated.

    """
    # tags2process = list(list_data.keys())
    print("Dataset to load: {}\n\n".format(list(list_data.keys())[0]))

    dim_tag = ""
    if "d3" in p.keys():
        if p["d3"]:
            dim_tag = "_3D"
        else:
            dim_tag = "_2D"            

    sc_matrix_collated, unique_barcodes = [], []
    build_pwd_matrix_collated, run_name, sc_labeled_collated = [], [], []

    for root_folder in list_data[dataset_name]["Folders"]:
        # [finds sample name]
        run_name.append(os.path.basename(os.path.dirname(root_folder)))

        # [makes list of files with Tables to load]
        # tries to load files from newer version of proceesingPipeline.py
        files_to_process_compatibility = glob.glob(root_folder + "/buildsPWDmatrix" + dim_tag + "_order*ROI*.ecsv")
        files_to_process = files_to_process_compatibility + glob.glob(root_folder + "/Trace" + dim_tag + "_barcode_*ROI*.ecsv")
        
        print("files_to_process: {}".format(root_folder +"/Trace" + dim_tag + "_barcode_ROI.ecsv"))
        
        if len(files_to_process) == 0:
            # it resorts to old format
            files_to_process = glob.glob(
                root_folder + "/buildsPWDmatrix" + dim_tag + "_*ROI*.ecsv"
            )
        else:
            print(
                "Found {} ECSV files in {}".format(len(files_to_process), root_folder)
            )

        # checks that something was found
        if len(files_to_process) > 0:
            print(">>> Loading {} results tables".format(len(files_to_process)))

            # [initializes variables]
            build_pwd_matrix = Table()
            file_order, file_order_stamp, file_time_stamp = (
                np.zeros(len(files_to_process), dtype=int),
                np.zeros(len(files_to_process), dtype=int),
                np.zeros(len(files_to_process)),
            )

            # [finds what order Table files should be loaded to agree with order in buildsPWDmatrix_HiMscMatrix.npy]
            for file_name, i_filename in zip(
                files_to_process, range(len(files_to_process))
            ):
                if "_order" in file_name:
                    for isplit in file_name.split("_"):
                        if "order" in isplit:
                            file_order_stamp[i_filename] = int(isplit.split(":")[1])
                            print(
                                "order {}= {}--> {}".format(
                                    i_filename,
                                    os.path.basename(file_name),
                                    file_order_stamp[i_filename],
                                )
                            )
                    file_time_stamp[i_filename] = os.path.getmtime(file_name)
                    choosing_time_stamp = False

                else:
                    file_time_stamp[i_filename] = os.path.getmtime(file_name)
                    choosing_time_stamp = True

            if choosing_time_stamp:
                file_order = np.argsort(file_time_stamp).astype(int)
            else:
                file_order = np.argsort(file_order_stamp).astype(int)

            # print('FileOrder: {}'.format(file_order))
            # [loads buildsPWDmatrix Tables]
            for i_filename in range(len(files_to_process)):
                file_name = files_to_process[file_order[i_filename]]
                new_build_pwd_matrix = Table.read(
                    file_name, format="ascii.ecsv"
                )  # ascii.ecsv
                build_pwd_matrix = vstack([build_pwd_matrix, new_build_pwd_matrix])
                print(
                    "[{}:{}:{}] From {}, Read: {} cells, Cummulative: {} cells".format(
                        i_filename,
                        file_order[i_filename],
                        file_time_stamp[file_order[i_filename]],
                        os.path.basename(file_name),
                        len(new_build_pwd_matrix),
                        len(build_pwd_matrix),
                    )
                )

            # [loads snd_assigned_cells.ecsv files if available]
            filename_snd_assigned_cells = (
                os.path.dirname(root_folder)
                + os.sep
                + "segmentedObjects/snd_assigned_cells.ecsv"
            )
            if os.path.exists(filename_snd_assigned_cells):
                print("Reading and processing: {}".format(filename_snd_assigned_cells))
                snd_assigned_cells = Table.read(
                    filename_snd_assigned_cells, format="ascii.ecsv"
                )

                # checks that table is not empty
                if len(snd_assigned_cells) > 0:
                    # attributes masks to single cells
                    sc_labeled, cuids_list = attributes_labels2cells(
                        snd_assigned_cells, build_pwd_matrix, label=p["label"]
                    )

                    # checks that at least one cell was found to have the label
                    if len(sc_labeled) == 0:
                        # if not available it makes a mock sc_labeled matrix so that pipeline always works
                        # sc_labeled = np.ones(len(build_pwd_matrix)).astype(int)
                        sc_labeled = np.zeros(len(build_pwd_matrix)).astype(int)

            else:
                # if not available it makes a mock sc_labeled matrix so that pipeline always works
                # sc_labeled = np.ones(len(build_pwd_matrix)).astype(int)
                sc_labeled = np.zeros(len(build_pwd_matrix)).astype(int)

            sc_labeled_collated.append(sc_labeled)

            # [loads and accumulates barcodes and scHiM matrix]
            filename_matrix = (
                root_folder + os.sep + "buildsPWDmatrix" + dim_tag + "_HiMscMatrix.npy"
            )
            filename_barcodes = (
                root_folder
                + os.sep
                + "buildsPWDmatrix"
                + dim_tag
                + "_uniqueBarcodes.ecsv"
            )

            if os.path.exists(filename_matrix):
                sc_matrix1 = np.load(filename_matrix)
                sc_matrix_collated.append(sc_matrix1)
            else:
                print("*** Error: could not find {}".format(filename_matrix))

            if os.path.exists(filename_barcodes):
                unique_barcodes.append(np.loadtxt(filename_barcodes).astype(int))
            else:
                print("*** Error: could not find {}".format(filename_barcodes))

            build_pwd_matrix_collated.append(build_pwd_matrix)

            print("\n>>>Merging root_folder: {}".format(root_folder))
            print("Cells added after merge: {}\n".format(sc_matrix1.shape[2]))
        else:
            print(
                "No file detected in the folder you provide: {}".format(
                    root_folder + "/buildsPWDmatrix_*ROI*.ecsv"
                )
            )
    print("{} datasets loaded\n".format(len(sc_matrix_collated)))

    return (
        sc_matrix_collated,
        unique_barcodes,
        build_pwd_matrix_collated,
        run_name,
        sc_labeled_collated,
    )


def load_sc_data_matlab(list_data, dataset_name, p):

    print("Dataset to load: {}\n\n".format(list(list_data.keys())[0]))

    sc_matrix_collated, unique_barcodes = [], []
    run_name, sc_labeled_collated = [], []

    for root_folder in list_data[dataset_name]["Folders"]:
        # [finds sample name]
        run_name.append(os.path.basename(os.path.dirname(root_folder)))

        # [loads and accumulates barcodes and scHiM matrix]
        filename_matrix = root_folder + os.sep + "HiMscMatrix.mat"
        filename_barcodes = root_folder + os.sep + "buildsPWDmatrix_uniqueBarcodes.ecsv"

        # loads barcodes
        if os.path.exists(filename_barcodes):
            unique_barcodes.append(np.loadtxt(filename_barcodes).astype(int))
            print(">>> Loaded {}".format(filename_matrix))
        else:
            print("*** Error: could not find {}".format(filename_barcodes))

        # loads SC matrix
        if os.path.exists(filename_matrix):
            data = loadmat(filename_matrix)
            sc_matrix1 = data["distanceMatrixCumulative"]
            sc_matrix_collated.append(sc_matrix1)
            print(
                ">>> Loaded: {}\n SC matrix shape: {}".format(
                    filename_matrix, sc_matrix1.shape
                )
            )
        else:
            print("*** Error: could not find {}".format(filename_matrix))

        # loads cell attributes
        cell_attributes_matrix = data["cellAttributesMatrix"]
        results_table = cell_attributes_matrix[0, :]

        sc_labeled = np.zeros(len(results_table))
        index_cells_with_label = [i_row for i_row, row in enumerate(results_table) if row > 0]
        sc_labeled[index_cells_with_label] = 1
        sc_labeled_collated.append(sc_labeled)

        print("\n>>>Merging root_folder: {}".format(root_folder))
        print("Cells added after merge: {}\n".format(sc_matrix1.shape[2]))

    print("{} datasets loaded\n".format(len(sc_matrix_collated)))

    return (
        sc_matrix_collated,
        unique_barcodes,
        run_name,
        sc_labeled_collated,
    )


def list_sc_to_keep(p, mask):
    # print("{}:{}".format(p["label"], p["action"]))
    if p["action"] == "all":
        try:
            cells_to_plot = range(len(mask))
        except TypeError:
            print(mask)
            cells_to_plot = range(mask.shape[0])
    elif p["action"] == "labeled":
        a = [i for i in range(len(mask)) if mask[i] == 1]
        cells_to_plot = a
    elif p["action"] == "unlabeled":
        a = [i for i in range(len(mask)) if mask[i] == 0]
        cells_to_plot = a

    print(
        ">> label: {}\t action:{}\t Ncells2plot:{}\t Ncells in sc_matrix:{}".format(
            p["label"], p["action"], max(cells_to_plot), len(mask)
        )
    )

    return cells_to_plot


def normalize_matrix(sc_matrix_wt):
    sc_matrix_wt_normalized = sc_matrix_wt
    n_bins = sc_matrix_wt.shape[0]

    for i_row in range(n_bins):
        row_sum = np.sum(sc_matrix_wt[i_row, :])
        for i_col in range(n_bins):
            sc_matrix_wt_normalized[i_row, i_col] = (
                sc_matrix_wt_normalized[i_row, i_col] / row_sum
            )
            sc_matrix_wt_normalized[i_col, i_row] = (
                sc_matrix_wt_normalized[i_col, i_row] / row_sum
            )
    return sc_matrix_wt_normalized


def plot_ensemble_3_way_contact_matrix(
    sc_matrix_collated,
    unique_barcodes,
    anchors,
    s_out,
    run_name,
    i_list_data,
    p,
    markdown_filename="tmp.md",
    dataset_name="",
):

    # combines matrices from different samples and calculates integrated contact probability matrix
    sc_matrix_all_datasets = []  # np.zeros((n_barcodes,n_barcodes))
    for i_sc_matrix_collated, i_unique_barcodes, mask, i_tag in zip(
        sc_matrix_collated, unique_barcodes, p["SClabeledCollated"], run_name
    ):
        cells_to_plot = list_sc_to_keep(p, mask)
        if len(cells_to_plot) > 0:
            if max(cells_to_plot) > i_sc_matrix_collated.shape[2]:
                print(
                    "Error: max in cells2plot {} in dataset {} is larger than the number of available cells {}".format(
                        max(cells_to_plot), i_tag, i_sc_matrix_collated.shape[2]
                    )
                )
            else:
                if len(sc_matrix_all_datasets) > 0:
                    sc_matrix_all_datasets = np.concatenate(
                        (sc_matrix_all_datasets, i_sc_matrix_collated[:, :, cells_to_plot]),
                        axis=2,
                    )
                else:
                    sc_matrix_all_datasets = i_sc_matrix_collated[:, :, cells_to_plot]
                common_set_unique_barcodes = i_unique_barcodes
        else:
            print(
                "Dataset: {} - {}  did not have any cell to plot".format(
                    dataset_name, i_tag
                )
            )

    # print(common_set_unique_barcodes)

    # loops over anchors
    for anchor in anchors:
        print("n_cells processed: {}".format(sc_matrix_all_datasets.shape[2]))

        # calculates the 3-way matrix for a given anchor
        sc_matrix = calculate_3_way_contact_matrix(
            sc_matrix_all_datasets,
            unique_barcodes,
            p["pixelSize"],
            anchor,
            s_out,
            threshold=i_list_data["ContactProbability_distanceThreshold"],
            norm="nonNANs",
        )  # norm: nonNANs (default)

        # output_filename = p['output_folder'] + os.sep + dataset_name + "_Cells:" + p['action'] + "_ensemble3wayContacts"
        output_filename = (
            p["outputFolder"]
            + os.sep
            + dataset_name
            + "_label:"
            + p["label"]
            + "_action:"
            + p["action"]
            + "_ensemble3wayContacts"
        )

        # outputs result
        output_filename += "_anchor_" + str(anchor)
        write_string_to_file(
            markdown_filename,
            "![]({})\n".format(output_filename + "_HiMmatrix.png"),
            "a",
        )

        c_scale = np.max(sc_matrix)
        # print(c_scale)
        # print(sc_matrix.shape)

        # plots 3-way matrix
        plot_matrix(
            sc_matrix,
            common_set_unique_barcodes,  # before unique_barcodes
            p["pixelSize"],
            c_m=i_list_data["ContactProbability_cm"],
            output_filename=output_filename,
            clim=c_scale,
            c_min=i_list_data["ContactProbability_cmin"],
            figtitle="3way contacts",
            cmtitle=s_out,
            n_cells=0,
            mode="counts",
        )  # twilight_shifted_r

        # saves matrices as individual files for further plotting
        root_output_filename = (
            p["outputFolder"]
            + os.sep
            + dataset_name
            + "_label:"
            + p["label"]
            + "_action:"
            + p["action"]
            + "_anchor:"
            + str(anchor)
        )
        np.save(root_output_filename + "_ensemble3wayContacts.npy", sc_matrix)


def calculate_3_way_contact_matrix(
    i_sc_matrix_collated,
    i_unique_barcodes,
    pixel_size,
    anchor,
    s_out,
    threshold=0.25,
    norm="nonNANs",
):

    n_x = n_y = i_sc_matrix_collated.shape[0]
    sc_matrix = np.zeros((n_x, n_y))

    # transform distance matrix from pixel to µm
    mat = pixel_size * i_sc_matrix_collated

    # print(n_x, n_y)
    for bait1 in range(n_x):
        for bait2 in range(n_y):
            if bait1 == bait2:
                continue

            # print("current bait1", bait1, "bait2", bait2)
            n_contacts, n_non_nan = get_multi_contact(mat, anchor, bait1, bait2, threshold)
            if s_out == "Counts":
                sc_matrix[bait1, bait2] = n_contacts
            elif s_out == "Probability":
                sc_matrix[bait1, bait2] = n_contacts / n_non_nan
            else:
                print("Unexpected s_out.")
                return -1

            # print(n_contacts / n_non_nan)
            # print(type(n_contacts), type(n_non_nan))

    sc_matrix[np.isnan(sc_matrix)] = 0  # set NaN to zero
    return sc_matrix


def get_multi_contact(mat, anchor, bait1, bait2, threshold):
    """
    Input:
    mat        : pwd matrix, including only the bins of used rts
    anchor     : anchor bin
    bait1      : bait bin #1
    bait2      : bait bin #2
    threshold  : contact threshold
    Output:
    n_contacts : number of contacts between bins anchor, bait1, and bait2
    n_non_nan   : number of cells where the distances anchor-bait1 and anchor-bait2 are present
    """

    mat_a = mat[anchor, bait1, :]
    mat_b = mat[anchor, bait2, :]

    # get fraction of points in quadrant
    n_1 = np.sum((mat_a < threshold) & (mat_b < threshold))
    tot_n = np.sum((~np.isnan(mat_a)) & (~np.isnan(mat_b)))

    return n_1, tot_n


def plot_single_pwd_matrice(
    sc_matrix_collated,
    unique_barcodes,
    run_name,
    i_list_data,
    p,
    markdown_filename="tmp.md",
    dataset_name="",
):
    # plots distance matrix for each dataset
    for i_sc_matrix_collated, i_unique_barcodes, i_tag, mask in zip(
        sc_matrix_collated, unique_barcodes, run_name, p["SClabeledCollated"]
    ):
        output_filename = (
            p["outputFolder"] + os.sep + i_tag + "_Cells:" + p["action"] + "_PWDmatrix"
        )

        # selects cels according to label
        cells_to_plot = list_sc_to_keep(p, mask)

        plot_matrix(
            i_sc_matrix_collated,
            i_unique_barcodes,
            p["pixelSize"],
            output_filename=output_filename,
            figtitle="PWD:" + dataset_name + i_tag,
            c_m=i_list_data["PWD_cm"],
            clim=i_list_data["PWD_clim"],
            mode=i_list_data["PWD_mode"],
            n_cells=i_sc_matrix_collated.shape[2],
            cells_to_plot=cells_to_plot,
        )  # twilight_shifted_r 1.4, mode: median KDE coolwarm terrain
        write_string_to_file(
            markdown_filename,
            "![]({})\n".format(output_filename + "_HiMmatrix.png"),
            "a",
        )


def plot_inverse_pwd_matrix(
    sc_matrix_collated,
    unique_barcodes,
    run_name,
    i_list_data,
    p,
    markdown_filename,
    dataset_name="",
):
    # plots inverse distance matrix for each dataset
    for i_sc_matrix_collated, i_unique_barcodes, i_tag, mask in zip(
        sc_matrix_collated, unique_barcodes, run_name, p["SClabeledCollated"]
    ):
        output_filename = (
            p["outputFolder"]
            + os.sep
            + i_tag
            + "_Cells:"
            + p["action"]
            + "_invPWDmatrix"
        )

        # selects cels according to label
        cells_to_plot = list_sc_to_keep(p, mask)
        # print('Dataset {} cells2plot: {}'.format(i_tag,cells_to_plot))

        plot_matrix(
            i_sc_matrix_collated,
            i_unique_barcodes,
            p["pixelSize"],
            c_m=i_list_data["iPWD_cm"],
            output_filename=output_filename,
            clim=i_list_data["iPWD_clim"],
            mode=i_list_data["iPWD_mode"],
            figtitle="inverse PWD:" + dataset_name + i_tag,
            cmtitle="inverse distance, 1/nm",
            inverse_matrix=True,
            n_cells=i_sc_matrix_collated.shape[2],
            cells_to_plot=cells_to_plot,
        )  # twilight_shifted_r, mode: median KDE
        write_string_to_file(
            markdown_filename,
            "![]({})\n".format(output_filename + "_HiMmatrix.png"),
            "a",
        )


def plot_single_contact_probability_matrix(
    sc_matrix_collated,
    unique_barcodes,
    run_name,
    i_list_data,
    p,
    markdown_filename="tmp.md",
    dataset_name="",
):
    # Plots contact probability matrices for each dataset
    if "minNumberContacts" in i_list_data.keys():
        min_number_contacts = i_list_data["minNumberContacts"]
    else:
        min_number_contacts = 0

    print("$ Min number contacts: {}".format(min_number_contacts))

    for i_sc_matrix_collated, i_unique_barcodes, i_tag, mask in zip(
        sc_matrix_collated, unique_barcodes, run_name, p["SClabeledCollated"]
    ):
        # selects cels according to label
        cells_to_plot = list_sc_to_keep(p, mask)

        if not cells_to_plot:
            break

        if max(cells_to_plot) > i_sc_matrix_collated.shape[2]:
            print(
                "Error with range in cells2plot {} as it is larger than the number of available cells {}".format(
                    max(cells_to_plot), i_sc_matrix_collated.shape[2]
                )
            )
        else:
            sc_matrix, n_cells = calculate_contact_probability_matrix(
                i_sc_matrix_collated[:, :, cells_to_plot],
                i_unique_barcodes,
                p["pixelSize"],
                threshold=i_list_data["ContactProbability_distanceThreshold"],
                min_number_contacts=min_number_contacts,
                norm="nonNANs",
            )  # norm: n_cells (default), nonNANs
            output_filename = (
                p["outputFolder"]
                + os.sep
                + dataset_name
                + i_tag
                + "_Cells:"
                + p["action"]
                + "_contactProbability"
            )

            print("Dataset {} cells2plot: {}".format(i_tag, cells_to_plot))
            c_scale = sc_matrix.max() / i_list_data["ContactProbability_scale"]

            plot_matrix(
                sc_matrix,
                i_unique_barcodes,
                p["pixelSize"],
                c_m=i_list_data["ContactProbability_cm"],
                output_filename=output_filename,
                c_min=i_list_data["ContactProbability_cmin"],
                clim=c_scale,
                figtitle="HiM:" + dataset_name + i_tag,
                cmtitle="probability",
                n_cells=n_cells,
                cells_to_plot=cells_to_plot,
            )  # twilight_shifted_r terrain coolwarm
            write_string_to_file(
                markdown_filename,
                "![]({})\n".format(output_filename + "_HiMmatrix.png"),
                "a",
            )


def fuses_sc_matrix_collated_from_datasets(
    sc_matrix_collated, unique_barcodes, p, run_name, i_list_data
):
    # combines matrices from different embryos and calculates integrated contact probability matrix

    sc_matrix_all_datasets = []
    cells_to_plot = []
    n_cells_total = 0

    for i_sc_matrix_collated, i_unique_barcodes, mask, i_tag in zip(
        sc_matrix_collated, unique_barcodes, p["SClabeledCollated"], run_name
    ):
        n_cells_total += mask.shape[0]
        # selects cels according to label
        cells_to_plot = list_sc_to_keep(p, mask)

        if len(cells_to_plot) > 0:
            if max(cells_to_plot) > i_sc_matrix_collated.shape[2]:
                print(
                    "Error: max in cells2plot {} in dataset {} is larger than the number of available cells {}".format(
                        max(cells_to_plot), i_tag, i_sc_matrix_collated.shape[2]
                    )
                )
            else:
                if len(sc_matrix_all_datasets) > 0:
                    sc_matrix_all_datasets = np.concatenate(
                        (sc_matrix_all_datasets, i_sc_matrix_collated[:, :, cells_to_plot]),
                        axis=2,
                    )
                else:
                    sc_matrix_all_datasets = i_sc_matrix_collated[:, :, cells_to_plot]

                common_set_unique_barcodes = i_unique_barcodes

    if p["saveMatrix"]:
        # write out the ensemble PWD map
        n_barcodes = sc_matrix_all_datasets.shape[0]
        pixel_size = p["pixelSize"]
        mean_sc_matrix = np.zeros((n_barcodes, n_barcodes))
        for bin1 in range(n_barcodes):
            for bin2 in range(n_barcodes):
                if bin1 != bin2:
                    (
                        maximum_kernel_distribution,
                        _,
                        _,
                        _,
                    ) = distribution_maximum_kernel_density_estimation(
                        sc_matrix_all_datasets,
                        bin1,
                        bin2,
                        pixel_size,
                        optimize_kernel_width=False,
                    )
                    mean_sc_matrix[bin1, bin2] = maximum_kernel_distribution

        output_filename = (
            p["outputFolder"]
            + os.sep
            + "CombinedMatrix_PWD_KDE"
            + ":"
            + list(i_list_data.keys())[0]
            + "_Cells:"
            + p["action"]
            + ".dat"
        )

        print(">>> Saving fused sc_matrix to {}".format(output_filename))

        np.savetxt(
            output_filename,
            mean_sc_matrix,
            fmt="%.4f",
            delimiter=" ",
            newline="\n",
            header="Combined pairwise distance map (kernel density estimator)",
            footer="",
            comments="# ",
            encoding=None,
        )

    if p["getStructure"]:
        ## multi-dimensional scaling to get coordinates from PWDs
        # make sure mean_sc_matrix is symmetric
        mean_sc_matrix = 0.5 * (mean_sc_matrix + np.transpose(mean_sc_matrix))
        # run metric mds
        verbosity = 0  # default: 0, quite verbose: 2
        mds = manifold.MDS(
            n_components=3,
            metric=True,
            n_init=20,
            max_iter=3000,
            verbose=verbosity,
            eps=1e-9,
            n_jobs=1,
            random_state=1,
            dissimilarity="precomputed",  # euclidean | precomputed
        )
        xyz = mds.fit(mean_sc_matrix).embedding_
        print(xyz)
        output_filename_pdb = (
            p["outputFolder"]
            + os.sep
            + "CombinedMatrix_PWD_KDE"
            + ":"
            + list(i_list_data.keys())[0]
            + "_Cells:"
            + p["action"]
            + "_python.pdb"
        )
        write_xyz_2_pdb(output_filename_pdb, xyz)

    return sc_matrix_all_datasets, common_set_unique_barcodes, cells_to_plot, n_cells_total


def plot_ensemble_contact_probability_matrix(
    sc_matrix_collated,
    unique_barcodes,
    run_name,
    i_list_data,
    p,
    markdown_filename="tmp.md",
    dataset_name="",
):

    if "minNumberContacts" in i_list_data.keys():
        min_number_contacts = i_list_data["minNumberContacts"]
    else:
        min_number_contacts = 0

    # combines matrices from different samples/datasers and calculates integrated contact probability matrix
    (
        sc_matrix_all_datasets,
        common_set_unique_barcodes,
        cells_to_plot,
        n_cells_total,
    ) = fuses_sc_matrix_collated_from_datasets(
        sc_matrix_collated, unique_barcodes, p, run_name, i_list_data
    )

    print(
        "n_cells selected / processed: {}/{}".format(
            sc_matrix_all_datasets.shape[2], n_cells_total
        )
    )

    # calculates contact probability matrix from merged samples/datasets
    sc_matrix, n_cells = calculate_contact_probability_matrix(
        sc_matrix_all_datasets,
        common_set_unique_barcodes,
        p["pixelSize"],
        threshold=i_list_data["ContactProbability_distanceThreshold"],
        min_number_contacts=min_number_contacts,
        norm=p["HiMnormalization"],
    )  # norm: n_cells (default), nonNANs

    # outputs line for MD file and sets output filename
    c_scale = sc_matrix.max() / i_list_data["ContactProbability_scale"]
    output_filename = (
        p["outputFolder"]
        + os.sep
        + dataset_name
        + "_Cells:"
        + p["action"]
        + "_ensembleContactProbability"
    )
    write_string_to_file(
        markdown_filename, "![]({})\n".format(output_filename + "_HiMmatrix.png"), "a"
    )

    # plots results
    plot_matrix(
        sc_matrix,
        common_set_unique_barcodes,
        p["pixelSize"],
        c_m=i_list_data["ContactProbability_cm"],
        output_filename=output_filename,
        clim=c_scale,
        c_min=i_list_data["ContactProbability_cmin"],
        figtitle="HiM counts",
        cmtitle="probability",
        n_cells=n_cells,
    )  # twilight_shifted_r

    # saves SC matrix in text format
    np.savetxt(
        p["outputFolder"]
        + os.sep
        + "CombinedMatrix"
        + ":"
        + list(i_list_data.keys())[0]
        + "_Cells:"
        + p["action"]
        + ".dat",
        sc_matrix,
        fmt="%.4f",
        delimiter=" ",
        newline="\n",
        header="Combined contact probability matrix",
        footer="",
        comments="# ",
        encoding=None,
    )

    # saves barcodes in text format
    np.savetxt(
        p["outputFolder"]
        + os.sep
        + "UniqueBarcodes"
        + ":"
        + list(i_list_data.keys())[0]
        + "_Cells:"
        + p["action"]
        + ".dat",
        common_set_unique_barcodes,
        fmt="%.4f",
        delimiter=" ",
        newline="\n",
        header="unique barcodes",
        footer="",
        comments="# ",
        encoding=None,
    )

    return sc_matrix, common_set_unique_barcodes


def shuffle_matrix(matrix, index):

    new_size = len(index)
    new_matrix = np.zeros((new_size, new_size))

    if new_size <= matrix.shape[0]:
        for i in range(new_size):
            for j in range(new_size):
                if index[i] < matrix.shape[0] and index[j] < matrix.shape[0]:
                    new_matrix[i, j] = matrix[index[i], index[j]]
    else:
        print(
            "Error: shuffle size {} is larger than matrix dimensions {}".format(
                new_size, matrix.shape[0]
            )
        )
        print("Shuffle: {} ".format(index))

    return new_matrix


def comp_func(mat_a, n_w):
    n_1 = mat_a.shape[0]

    somme_short = np.zeros((n_1, 1))
    signal1 = np.zeros((n_1, 1))
    n_int = np.zeros((n_1, 1))

    for i in range(0, n_1):
        if i <= n_w:
            p_1 = n_1 + i - n_w
            p_2 = i + n_w
            for k in range(0, n_1):
                if k <= p_2 or k >= p_1:
                    somme_short[i] = somme_short[i] + mat_a[i, k]
                    n_int[i] = n_int[i] + 1

        elif (n_1 - n_w) <= i:
            p_1 = i - n_w
            p_2 = i + n_w - n_1
            for k in range(0, n_1):
                if k <= p_2 or k >= p_1:
                    somme_short[i] = somme_short[i] + mat_a[i, k]
                    n_int[i] = n_int[i] + 1

        else:
            p_1 = i - n_w
            p_2 = i + n_w
            for k in range(0, n_1):
                if p_1 <= k and k <= p_2:
                    somme_short[i] = somme_short[i] + mat_a[i, k]
                    n_int[i] = n_int[i] + 1

    signal1 = somme_short
    return signal1


def plot_scalogram(matrix2plot, output_filename=""):
    matrix_size = matrix2plot.shape[0]

    scales = range(0, 6)
    comp_scale1 = np.zeros((len(scales), matrix_size))
    for i_i, n_w in enumerate(scales):
        c = comp_func(matrix2plot, n_w)
        comp_scale1[i_i,] = c.T

    # definitions for the axes
    left, width = 0.065, 0.65
    bottom, height = 0.1, 0.8
    bottom_h = left_h = left + width + 0.02

    ax1_ = [0, 0.4, width, height]
    ax2_ = [left, bottom, width, 0.2]

    # fig = plt.figure()
    fig1 = plt.figure(constrained_layout=True)
    spec1 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig1)
    f_1 = fig1.add_subplot(spec1[1, 0])  # 16
    # ax2 = plt.axes(ax2_)#,sharex=ax1)
    f_1.contourf(
        comp_scale1,
        vmin=0.0,
        vmax=1.0,
        levels=[0, 0.15, 0.30, 0.45, 0.6, 0.75, 1.0],
        extent=[0, matrix2plot.shape[0], 0, 400],
        cmap="rainbow",
    )
    im_color_bar = contourf(
        comp_scale1,
        vmin=0.0,
        vmax=1.0,
        levels=[0, 0.15, 0.30, 0.45, 0.6, 0.75, 1.0],
        extent=[0, matrix2plot.shape[0], 0, 400],
        cmap="rainbow",
    )
    f_1.set_xlim([0, matrix_size])  # list(data['unique_barcodes']))
    f_1.set_ylabel("Scales")
    bounds = [0.0, 0.5, 1.0]
    cb1 = colorbar(
        im_color_bar,
        shrink=1,
        orientation="vertical",
        ticks=bounds,
        spacing="proportional",
        pad=0.04,
    )

    # ax1 = plt.axes(ax1_)
    # pos=ax1.imshow(matrix2plot,cmap='coolwarm')
    # ax1.set_xlabel("barcode #")
    # ax1.set_ylabel("barcode #")
    # ax1.set_xlim([-0.5 , matrix_size-.5])#list(data['unique_barcodes']))
    # ax1.set_ylim([matrix_size-.5,-0.5])#list(data['unique_barcodes']))
    # cbar = plt.colorbar(pos, fraction=0.046, pad=0.04)

    # ax2 = plt.axes(ax2_)#,sharex=ax1)
    # ax2.contourf(  comp_scale1, vmin=0.0, vmax= 1.0,levels=[0, .15, .30, .45, 0.6,0.75, 1.0],extent=[0,matrix2plot.shape[0],0,400],cmap="rainbow");
    # im_color_bar=contourf(  comp_scale1, vmin=0.0, vmax= 1.0,levels=[0, .15, .30, .45, 0.6,0.75, 1.0],extent=[0,matrix2plot.shape[0],0,400],cmap="rainbow");
    # ax2.set_xlim([-.5,matrix_size-.5])#list(data['unique_barcodes']))
    # ax2.set_ylabel("Scales");
    # bounds = [0.,0.5, 1.0];
    # cb1 = colorbar(im_color_bar,shrink = 1, orientation="vertical",ticks=bounds, spacing='proportional',pad=0.04);

    if output_filename:
        plt.savefig(output_filename)
        print("Output scalogram: {}".format(output_filename))


def write_xyz_2_pdb(file_name, xyz):
    # writes xyz coordinates to a PDB file wth pseudoatoms
    # file_name : string of output file path, e.g. '/foo/bar/test2.pdb'
    # xyz      : n-by-3 numpy array with atom coordinates

    n_atoms = xyz.shape[0]
    with open(file_name, mode="w+", encoding="utf-8") as fid:
        ## atom coordinates
        txt = "HETATM  {: 3d}  C{:02d} PSD P   1      {: 5.3f}  {: 5.3f}  {: 5.3f}  0.00  0.00      PSDO C  \n"
        for i in range(n_atoms):
            fid.write(txt.format(i + 1, i + 1, xyz[i, 0], xyz[i, 1], xyz[i, 2]))

        ## connectivity
        txt1 = "CONECT  {: 3d}  {: 3d}\n"
        txt2 = "CONECT  {: 3d}  {: 3d}  {: 3d}\n"
        # first line of connectivity
        fid.write(txt1.format(1, 2))
        # consecutive lines
        for i in range(2, n_atoms):
            fid.write(txt2.format(i, i - 1, i + 1))
        # last line
        fid.write(txt1.format(i + 1, i))

        print("Done writing {:s} with {:d} atoms.".format(file_name, n_atoms))


def plot_distance_histograms(
    sc_matrix_collated,
    pixel_size,
    output_filename="test",
    log_name_md="log.md",
    mode="hist",
    limit_n_plots=10,
    kernel_width=0.25,
    optimize_kernel_width=False,
    max_distance=4.0,
):

    if not is_notebook():
        n_plots_x = n_plots_y = sc_matrix_collated.shape[0]
    else:
        if limit_n_plots == 0:
            n_plots_x = n_plots_y = sc_matrix_collated.shape[0]
        else:
            n_plots_x = n_plots_y = min(
                [limit_n_plots, sc_matrix_collated.shape[0]]
            )  # sets a max of subplots if you are outputing to screen!

    bins = np.arange(0, max_distance, 0.25)

    size_x, size_y = n_plots_x * 4, n_plots_y * 4

    fig, axs = plt.subplots(
        figsize=(size_x, size_y), ncols=n_plots_x, nrows=n_plots_y, sharex=True
    )

    for i in trange(n_plots_x):
        for j in range(n_plots_y):
            if i != j:
                # print('Printing [{}:{}]'.format(i,j))
                if mode == "hist":
                    axs[i, j].hist(pixel_size * sc_matrix_collated[i, j, :], bins=bins)
                else:
                    (
                        max_kde,
                        distance_distribution,
                        kde,
                        x_d,
                    ) = distribution_maximum_kernel_density_estimation(
                        sc_matrix_collated,
                        i,
                        j,
                        pixel_size,
                        optimize_kernel_width=optimize_kernel_width,
                        kernel_width=kernel_width,
                        max_distance=max_distance,
                    )
                    axs[i, j].fill_between(x_d, kde, alpha=0.5)
                    axs[i, j].plot(
                        distance_distribution,
                        np.full_like(distance_distribution, -0.01),
                        "|k",
                        markeredgewidth=1,
                    )
                    axs[i, j].vlines(max_kde, 0, kde.max(), colors="r")

            axs[i, j].set_xlim(0, max_distance)
            axs[i, j].set_yticklabels([])

    plt.xlabel("distance, um")
    plt.ylabel("counts")

    file_extension = output_filename.split(".")[-1]

    if len(file_extension) == 3:
        file_name = output_filename + "_PWDhistograms." + file_extension
    else:
        file_name = output_filename + "_PWDhistograms.png"

    print("Output figure: {}\n".format(file_name))
    plt.savefig(file_name)

    if not is_notebook():
        plt.close()

    write_string_to_file(
        log_name_md, "![]({})\n".format(output_filename + "_PWDhistograms.png"), "a"
    )


def plot_matrix(
    sc_matrix_collated,
    unique_barcodes,
    pixel_size,
    number_rois=1,
    output_filename="test",
    log_name_md="log.md",
    clim=1.4,
    c_m="seismic",
    figtitle="PWD matrix",
    cmtitle="distance, um",
    n_cells=0,
    mode="median",
    inverse_matrix=False,
    c_min=0,
    cells_to_plot=[],
    filename_ending="_HiMmatrix.png",
):

    n_barcodes = sc_matrix_collated.shape[0]

    ######################################################
    # Calculates ensemble matrix from single cell matrices
    ######################################################

    if len(sc_matrix_collated.shape) == 3:

        # matrix is 3D and needs combining SC matrices into an ensemble matrix
        if len(cells_to_plot) == 0:
            cells_to_plot = range(sc_matrix_collated.shape[2])

        mean_sc_matrix, keep_plotting = calculate_ensemble_pwd_matrix(
            sc_matrix_collated, pixel_size, cells_to_plot, mode=mode
        )

    else:

        # matrix is 2D and does not need further treatment
        if mode == "counts":
            mean_sc_matrix = sc_matrix_collated
            keep_plotting = True
        else:
            mean_sc_matrix = pixel_size * sc_matrix_collated
            keep_plotting = True

    if keep_plotting:
        # no errors occurred

        # Calculates the inverse distance matrix if requested in the argument.
        if inverse_matrix:
            mean_sc_matrix = np.reciprocal(mean_sc_matrix)

        # plots figure
        fig = plt.figure(figsize=(10, 10))
        pos = plt.imshow(mean_sc_matrix, cmap=c_m)  # colormaps RdBu seismic
        plt.xlabel("barcode #")
        plt.ylabel("barcode #")
        plt.title(
            figtitle
            + " | "
            + str(mean_sc_matrix.shape[0])
            + " barcodes | n="
            + str(n_cells)
            + " | rois="
            + str(number_rois)
        )
        # print("matrix size: {} | barcodes:{}".format(sc_matrix_collated.shape[0],list(unique_barcodes)))
        plt.xticks(np.arange(sc_matrix_collated.shape[0]), unique_barcodes)
        plt.yticks(np.arange(sc_matrix_collated.shape[0]), unique_barcodes)
        cbar = plt.colorbar(pos, fraction=0.046, pad=0.04)
        cbar.minorticks_on()
        cbar.set_label(cmtitle)
        plt.clim(c_min, clim)

        if len(output_filename.split(".")) > 1:
            if output_filename.split(".")[1] != "png":
                if len(output_filename.split(".")[1]) == 3:
                    # keeps original extension
                    out_fn = output_filename
                    plt.savefig(output_filename)
                else:
                    # most likely the full filname contains other '.' in addition to that in the extension
                    out_fn = output_filename + filename_ending
                    plt.savefig(out_fn)
            else:
                out_fn = output_filename.split(".")[0] + filename_ending
                plt.savefig(out_fn)
        else:
            out_fn = output_filename + filename_ending
            plt.savefig(out_fn)

        if not is_notebook():
            plt.close()
        if "png" not in out_fn:
            out_fn += ".png"
        write_string_to_file(log_name_md, "![]({})\n".format(out_fn), "a")
    else:
        # errors during pre-processing
        print("Error plotting figure. Not executing script to avoid crash.")

    return meanSCmatrix

def calculate_contact_probability_matrix(
    i_sc_matrix_collated,
    i_unique_barcodes,
    pixel_size,
    threshold=0.25,
    norm="n_cells",
    min_number_contacts=0,
):

    n_x = n_y = i_sc_matrix_collated.shape[0]
    n_cells = i_sc_matrix_collated.shape[2]
    sc_matrix = np.zeros((n_x, n_y))

    for i in range(n_x):
        for j in range(n_y):
            if i != j:
                distance_distribution = pixel_size * i_sc_matrix_collated[i, j, :]

                number_contacts = distance_distribution.squeeze().shape[0] - len(
                    np.nonzero(np.isnan(distance_distribution))[0]
                )

                if number_contacts < min_number_contacts:
                    print(
                        "$ Rejected {}-{} because number contacts: {} < {}".format(
                            i, j, number_contacts, min_number_contacts
                        )
                    )
                    probability = 0.0
                else:
                    # normalizes # of contacts by the # of cells
                    if norm == "n_cells":
                        probability = (
                            len(np.nonzero(distance_distribution < threshold)[0])
                            / n_cells
                        )

                    # normalizes # of contacts by the # of PWD detected in each bin
                    elif norm == "nonNANs":
                        number_nans = len(np.nonzero(np.isnan(distance_distribution))[0])
                        if n_cells == number_nans:
                            probability = np.nan
                        else:
                            probability = len(
                                np.nonzero(distance_distribution < threshold)[0]
                            ) / (n_cells - number_nans)

                sc_matrix[i, j] = probability

    return sc_matrix, n_cells


# @jit(nopython=True)
def find_optimal_kernel_width(distance_distribution):
    bandwidths = 10 ** np.linspace(-1, 1, 100)
    grid = GridSearchCV(
        KernelDensity(kernel="gaussian"), {"bandwidth": bandwidths}, cv=LeaveOneOut()
    )
    grid.fit(distance_distribution[:, None])
    return grid.best_params_


# @jit(nopython=True)
def retrieve_kernel_density_estimator(
    distance_distribution_0, x_d, optimize_kernel_width=False, kernel_width=0.25
):
    """
    Gets the kernel density function and maximum from a distribution of PWD distances

    Parameters
    ----------
    distance_distribution_0 : nd array
        List of PWD distances.
    x_d : nd array
        x grid.
    optimize_kernel_width : Boolean, optional
        whether to optimize bandwidth. The default is False.

    Returns
    -------
    np array
        kde distribution.
    np array
        Original distribution without NaNs

    """

    nan_array = np.isnan(distance_distribution_0)

    not_nan_array = ~nan_array

    distance_distribution = distance_distribution_0[not_nan_array]

    # instantiate and fit the KDE model
    if optimize_kernel_width:
        kernel_width = find_optimal_kernel_width(distance_distribution)["bandwidth"]
    else:
        kernel_width = kernel_width

    kde = KernelDensity(bandwidth=kernel_width, kernel="gaussian")

    # makes sure the list is not full of NaNs.
    if distance_distribution.shape[0] > 0:
        kde.fit(distance_distribution[:, None])
    else:
        return np.array([0]), np.array([0])

    # score_samples returns the log of the probability density
    logprob = kde.score_samples(x_d[:, None])

    return logprob, distance_distribution


# @jit(nopython=True)
def distribution_maximum_kernel_density_estimation(
    sc_matrix_collated,
    bin1,
    bin2,
    pixel_size,
    optimize_kernel_width=False,
    kernel_width=0.25,
    max_distance=4.0,
):
    """
    calculates the kernel distribution and its maximum from a set of PWD distances

    Parameters
    ----------
    sc_matrix_collated : np array 3 dims
        SC PWD matrix.
    bin1 : int
        first bin.
    bin2 : int
        first bin.
    pixel_size : float
        pixel size in um
    optimize_kernel_width : Boolean, optional
        does kernel need optimization?. The default is False.

    Returns
    -------
    float
        maximum of kernel.
    np array
        list of PWD distances used.
    np array
        kernel distribution.
    x_d : np array
        x grid.

    """
    distance_distribution_unlimited = (
        pixel_size * sc_matrix_collated[bin1, bin2, :]
    )  # full distribution
    distance_distribution_unlimited = distance_distribution_unlimited[
        ~np.isnan(distance_distribution_unlimited)
    ]  # removes nans

    if bin1 == bin2:
        # protection agains bins in the diagonal
        distance_distribution_0 = distance_distribution_unlimited
    else:
        # removes values larger than max_distance
        distance_distribution_0 = distance_distribution_unlimited[
            np.nonzero(distance_distribution_unlimited < max_distance)
        ]
    x_d = np.linspace(0, max_distance, 2000)

    # checks that distribution is not empty
    if distance_distribution_0.shape[0] > 0:
        logprob, distance_distribution = retrieve_kernel_density_estimator(
            distance_distribution_0, x_d, optimize_kernel_width, kernel_width
        )
        if logprob.shape[0] > 1:
            kernel_distribution = 10 * np.exp(logprob)
            maximum_kernel_distribution = x_d[np.argmax(kernel_distribution)]
            return (
                maximum_kernel_distribution,
                distance_distribution,
                kernel_distribution,
                x_d,
            )
        else:
            return np.nan, np.zeros(x_d.shape[0]), np.zeros(x_d.shape[0]), x_d
    else:
        return np.nan, np.zeros(x_d.shape[0]), np.zeros(x_d.shape[0]), x_d


def get_rg_from_pwd(pwd_matrix_0, min_number_pwd=4, threshold=6):
    """
    Calculates the Rg from a 2D pairwise distance matrix
    while taking into account that some of the PWD might be NaN

    PWDmatrix:       numpy array, NxN
    minFracNotNaN:   require a minimal fraction of PWDs to be not NaN, return NaN otherwise

    for the math, see https://en.wikipedia.org/wiki/Radius_of_gyration#Molecular_applications
    """

    pwd_matrix = pwd_matrix_0.copy()

    # check that pwd_matrix is of right shape
    if pwd_matrix.ndim != 2:
        raise SystemExit(
            "get_rg_from_pwd: Expected 2D input but got {}D.".format(pwd_matrix.ndim)
        )
    if pwd_matrix.shape[0] != pwd_matrix.shape[1]:
        raise SystemExit("get_rg_from_pwd: Expected square matrix as input.")

    # make sure the diagonal is NaN
    np.fill_diagonal(pwd_matrix, np.NaN)

    # filters out PWD
    pwd_matrix[pwd_matrix > threshold] = np.nan

    # get the number of PWDs that are not NaN
    num_pwds = pwd_matrix.shape[0] * (pwd_matrix.shape[0] - 1) / 2
    num_not_nan = (
        np.sum(~np.isnan(pwd_matrix)) / 2
    )  # default is to compute the sum of the flattened array

    if num_not_nan < min_number_pwd:
        return np.NaN

    # calculate Rg
    sqr = np.square(pwd_matrix)
    sqr = np.nansum(sqr)  # default is to compute the sum of the flattened array

    rg_sq = sqr / (2 * (2 * num_not_nan + pwd_matrix.shape[0]))  # replaces 1/(2*N^2)

    radius_gyration = np.sqrt(rg_sq)

    return radius_gyration


def get_detection_eff_barcodes(sc_matrix_collated):
    """
    Return the detection efficiency of all barcodes.
    Assumes a barcode is detected as soon as one PWD with this barcode is detected.
    """

    # check that pwd_matrix is of right shape
    if sc_matrix_collated.ndim != 3:
        raise SystemExit(
            "getBarcodeEff: Expected 3D input but got {}D.".format(
                sc_matrix_collated.ndim
            )
        )
    if sc_matrix_collated.shape[0] != sc_matrix_collated.shape[1]:
        raise SystemExit(
            "getBarcodeEff: Expected axis 0 and 1 to have the same length."
        )

    # make sure the diagonal is NaN
    for i in range(sc_matrix_collated.shape[0]):
        sc_matrix_collated[i, i, :] = np.NaN

    # calculate barcode efficiency
    n_cells = sc_matrix_collated.shape[2]

    eff = np.sum(~np.isnan(sc_matrix_collated), axis=0)
    eff[eff > 1] = 1

    eff0 = eff.copy()
    n_cells_2 = np.nonzero(np.sum(eff0, axis=0) > 2)[0].shape[0]

    eff = np.sum(eff, axis=-1)  # sum over all cells

    eff = eff / n_cells_2

    print("\n\n *** n_cells={} | n_cells_2={}".format(n_cells, n_cells_2))
    return eff


def get_barcodes_per_cell(sc_matrix_collated):
    """
    Returns the number of barcodes that were detected in each cell of sc_matrix_collated.
    """

    # make sure the diagonal is NaN
    for i in range(sc_matrix_collated.shape[0]):
        sc_matrix_collated[i, i, :] = np.NaN

    num_barcodes = np.sum(~np.isnan(sc_matrix_collated), axis=0)
    num_barcodes[num_barcodes > 1] = 1
    num_barcodes = np.sum(num_barcodes, axis=0)

    return num_barcodes


def get_coordinates_from_pwd_matrix(matrix):
    ## multi-dimensional scaling to get coordinates from PWDs
    # make sure mean_sc_matrix is symmetric
    matrix = 0.5 * (matrix + np.transpose(matrix))
    # run metric mds
    verbosity = 0  # default: 0, quite verbose: 2
    mds = manifold.MDS(
        n_components=3,
        metric=True,
        n_init=20,
        max_iter=3000,
        verbose=verbosity,
        eps=1e-9,
        n_jobs=1,
        random_state=1,
        dissimilarity="precomputed",  # euclidean | precomputed
    )

    xyz = mds.fit(matrix).embedding_

    return xyz


def sort_cells_by_number_pwd(him_data):

    # sc_matrix = him_data.data["SCmatrixCollated"]
    sc_matrix = him_data.sc_matrix_selected

    n_cells = sc_matrix.shape[2]
    # print("Number of cells loaded: {}".format(n_cells))

    # finds the number of barcodes detected per cell.
    n_barcode_per_cell = []
    values = []
    dtype = [("cellID", int), ("n_pwd", int)]

    for i_cell in range(n_cells):
        sc_matrix_cell = sc_matrix[:, :, i_cell]
        n_pwd = int(np.count_nonzero(~np.isnan(sc_matrix_cell)) / 2)
        n_barcode_per_cell.append(n_pwd)
        values.append((i_cell, n_pwd))

    values_array = np.array(values, dtype=dtype)  # create a structured array
    sorted_values = np.sort(values_array, order="n_pwd")

    return sc_matrix, sorted_values, n_cells


def kde_fit(x, x_d, bandwidth=0.2, kernel="gaussian"):

    kde = KernelDensity(bandwidth=bandwidth, kernel="gaussian")
    kde.fit(x[:, None])

    logprob = kde.score_samples(x_d[:, None])

    return logprob, kde


def calculate_ensemble_pwd_matrix(sc_matrix, pixel_size, cells_to_plot, mode="median"):
    """
    performs a KDE or median to calculate the max of the PWD distribution

    Parameters
    ----------
    sc_matrix : TYPE
        DESCRIPTION.
    pixel_size : TYPE
        DESCRIPTION.

    Returns
    -------
    matrix = 2D npy array.

    """

    n_barcodes = sc_matrix.shape[0]
    # cells_to_plot = range(sc_matrix.shape[2])

    mean_sc_matrix = np.zeros((n_barcodes, n_barcodes))

    if mode == "median":
        # calculates the median of all values #
        #######################################
        if max(cells_to_plot) > sc_matrix.shape[2]:
            print(
                "Error with range in cells2plot {} as it is larger than the number of available cells {}".format(
                    max(cells_to_plot), sc_matrix.shape[2]
                )
            )
            keep_plotting = False
        else:
            mean_sc_matrix = pixel_size * np.nanmedian(sc_matrix[:, :, cells_to_plot], axis=2)
            keep_plotting = True

    elif mode == "KDE":
        keep_plotting = True

        if max(cells_to_plot) > sc_matrix.shape[2]:
            print(
                "Error with range in cells2plot {} as it is larger than the number of available cells {}".format(
                    max(cells_to_plot), sc_matrix.shape[2]
                )
            )
            keep_plotting = False
        else:
            for bin1 in trange(n_barcodes):
                for bin2 in range(n_barcodes):
                    if bin1 != bin2:
                        # print(f"cells_to_plot:{cells_to_plot}, ncells:{sc_matrix.shape}")
                        (
                            maximum_kernel_distribution,
                            _,
                            _,
                            _,
                        ) = distribution_maximum_kernel_density_estimation(
                            sc_matrix[:, :, cells_to_plot],
                            bin1,
                            bin2,
                            pixel_size,
                            optimize_kernel_width=False,
                        )
                        mean_sc_matrix[bin1, bin2] = maximum_kernel_distribution

    return mean_sc_matrix, keep_plotting
