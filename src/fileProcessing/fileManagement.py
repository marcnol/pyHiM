#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and functions for file management
"""

import glob
import json
import logging
import multiprocessing
import os
import re
from datetime import datetime
from os import path
from warnings import warn

import numpy as np
from dask.distributed import Client, LocalCluster, get_client

from core.pyhim_logging import print_log

# =============================================================================
# CLASSES
# =============================================================================


class Session:
    def __init__(self, root_folder, name="dummy"):
        now = datetime.now()
        session_root_name = now.strftime("%d%m%Y_%H%M%S")
        self.file_name = root_folder + os.sep + "Session_" + session_root_name + ".json"
        self.name = name
        self.data = {}

    # loads existing session
    def load(self):
        self.data = load_json(self.file_name)
        print(f"Session information read: {self.file_name}")

    # saves session to file
    def save(self):
        save_json(self.file_name, self.data)
        print(f"> Saved json session file to {self.file_name}")

    # add new task to session
    def add(self, key, value):
        if key not in self.data:
            self.data[key] = value
        else:
            self.data[key] = [self.data[key], value]

    def clear_data(self):
        self.data = {}


class Parameters:
    def __init__(
        self,
        root_folder="./",
        label="",
        file_name="infoList.json",
        stardist_basename=None,
    ):
        self.file_name = file_name
        self.label = label
        self.param_file = root_folder + os.sep + file_name
        self.files_to_process = []
        self.param_dict = {}

        self.initialize_standard_parameters()
        self.convert_parameter_file(self.param_file, self.label)
        if stardist_basename is not None:
            self.param_dict["segmentedObjects"]["stardist_basename"] = stardist_basename
        self.param_dict["rootFolder"] = root_folder
        self.file_parts = {}

    def get_param_section(self, param_section=""):
        if not param_section:
            return self.param_dict
        else:
            return self.param_dict[param_section]

    def initialize_standard_parameters(self):
        self.param_dict = {
            "common": {
                "acquisition": {
                    "positionROIinformation": 3,
                    "fileNameRegExp": "scan_(?P<runNumber>[0-9]+)_(?P<cycle>[\\w|-]+)_(?P<roi>[0-9]+)_ROI_converted_decon_(?P<channel>[\\w|-]+).tif",
                    "DAPI_channel": "ch00",
                    "fiducialDAPI_channel": "ch01",
                    "RNA_channel": "ch02",
                    "fiducialBarcode_channel": "ch00",
                    "fiducialMask_channel": "ch00",
                    "barcode_channel": "ch01",
                    "mask_channel": "ch01",
                    "label_channel": "ch00",  # in future this field will contain the ch for the label. This parameter will supersed the individual channel fields above.
                    "label_channel_fiducial": "ch01",  # in future this field will contain the ch for the label fiducial. This parameter will supersed the individual channel fields above.
                    "pixelSizeXY": 0.1,
                    "zBinning": 2,
                    "parallelizePlanes": False,  # if True it will parallelize inner loops (plane by plane). Otherwise outer loops (e.g. file by file)
                    "pixelSizeZ": 0.25,
                },  # barcode, fiducial
                "zProject": {
                    "folder": "zProject",  # output folder
                    "operation": "skip",  # overwrite, skip
                    "mode": "full",  # full, manual, automatic, laplacian
                    "blockSize": 256,
                    "display": True,
                    "saveImage": True,
                    "zmin": 1,
                    "zmax": 59,
                    "zwindows": 15,
                    "windowSecurity": 2,
                    "zProjectOption": "MIP",  # sum or MIP
                },
                "alignImages": {
                    "folder": "alignImages",  # output folder
                    "operation": "overwrite",  # overwrite, skip
                    "outputFile": "alignImages",
                    "referenceFiducial": "RT27",
                    "localAlignment": "block3D",  # options: None, mask2D, block3D
                    "alignByBlock": True,  # alignByBlock True will perform block alignment
                    "tolerance": 0.1,  # Used in blockAlignment to determine the % of error tolerated
                    "lower_threshold": 0.999,  # lower threshold to adjust image intensity levels before xcorrelation for alignment in 2D
                    "higher_threshold": 0.9999999,  # higher threshold to adjust image intensity levels before xcorrelation for alignment in 2D
                    "3D_lower_threshold": 0.9,  # lower threshold to adjust image intensity levels before xcorrelation for Alignment3D
                    "3D_higher_threshold": 0.9999,  # higher threshold to adjust image intensity levels before xcorrelation for Alignment3D
                    "background_sigma": 3.0,  # used to remove inhom background
                    "localShiftTolerance": 1,
                    "blockSize": 256,
                },
                "buildsPWDmatrix": {
                    "folder": "buildsPWDmatrix",  # output folder
                    "tracing_method": [
                        "masking",
                        "clustering",
                    ],  # available methods: masking, clustering
                    "mask_expansion": 8,  # Expands masks until they collide by a max of 'mask_expansion' pixels
                    "flux_min": 10,  # min flux to keeep object
                    "flux_min_3D": 0.1,  # min flux to keeep object
                    "kd_tree_distance_threshold_mum": 1,  # distance threshold used to build KDtree
                    "colormaps": {
                        "PWD_KDE": "terrain",
                        "PWD_median": "terrain",
                        "contact": "coolwarm",
                        "Nmatrix": "Blues",
                    },  # colormaps used for plotting matrices
                    "toleranceDrift": [
                        3,
                        1,
                        1,
                    ],  # zxy tolerance used for block drift correction, in px
                    "remove_uncorrected_localizations": True,  # if True it will removed uncorrected localizations, otherwise they will remain uncorrectd.
                },
                "segmentedObjects": {
                    "folder": "segmentedObjects",  # output folder
                    "operation": "2D,3D",  # options: 2D or 3D
                    "outputFile": "segmentedObjects",
                    "Segment3D": "overwrite",
                    "background_method": "inhomogeneous",  # flat or inhomogeneous or stardist
                    "stardist_basename": "/mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks",
                    "stardist_network": "stardist_nc14_nrays:64_epochs:40_grid:2",  # network for 2D barcode segmentation
                    "stardist_network3D": "stardist_nc14_nrays:64_epochs:40_grid:2",  # network for 3D barcode segmentation
                    "tesselation": True,  # tesselates masks
                    "background_sigma": 3.0,  # used to remove inhom background
                    "threshold_over_std": 1.0,  # threshold used to detect sources
                    "fwhm": 3.0,  # source size in px
                    "brightest": 1100,  # max number of sources segmented per FOV
                    "intensity_min": 0,  # min int to keep object
                    "intensity_max": 59,  # max int to keeep object
                    "area_min": 50,  # min area to keeep object
                    "area_max": 500,  # max area to keeep object
                    "3Dmethod": "thresholding",  # options: 'thresholding' or 'stardist', 'zASTROPY', 'zProfile'
                    "residual_max": 2.5,  # z-profile Fit: max residuals to keeep object
                    "sigma_max": 5,  # z-profile Fit: max sigma 3D fitting to keeep object
                    "centroidDifference_max": 5,  # z-profile Fit: max diff between Moment and z-gaussian fits to keeep object
                    "3DGaussianfitWindow": 3,  # z-profile Fit: window size to extract subVolume, px. 3 means subvolume will be 7x7.
                    "3dAP_window": 5,  # constructs a YZ image by summing from xPlane-window:xPlane+window
                    "3dAP_flux_min": 2,  # # threshold to keep a source detected in YZ
                    "3dAP_brightest": 100,  # number of sources sought in each YZ plane
                    "3dAP_distTolerance": 1,  # px dist to attribute a source localized in YZ to one localized in XY
                    "3D_threshold_over_std": 5,
                    "3D_sigma": 3,
                    "3D_boxSize": 32,
                    "3D_filter_size": 3,
                    "3D_area_min": 10,
                    "3D_area_max": 250,
                    "3D_nlevels": 64,
                    "3D_contrast": 0.001,
                    "3D_psf_z": 500,
                    "3D_psf_yx": 200,
                    "3D_lower_threshold": 0.99,
                    "3D_higher_threshold": 0.9999,
                    "reducePlanes": True,  # if true it will calculate focal plane and only use a region around it for segmentSources3D, otherwise will use the full stack
                },
            },
            "labels": {
                "fiducial": {"order": 1},
                "barcode": {"order": 2},
                "DAPI": {"order": 3},
                "mask": {"order": 5},
                "RNA": {"order": 4},
            },
        }

    def convert_parameter_file(self, param_file, label_selected):
        param_from_file = load_parameters_file(param_file)

        if param_from_file is None:
            raise ValueError("No infoList.json file found")

        converted_param = param_from_file["common"]

        labels = param_from_file["labels"]

        ordered_list = [" "] * len(labels.keys())
        for _, label in enumerate(labels.keys()):
            order = labels[label]["order"]
            ordered_list[order - 1] = label

        converted_param["labels"] = ordered_list

        # need to add keys not present in common dict
        if len(label_selected) > 0:
            number_keys_ammended = 0
            for key in param_from_file["labels"][label_selected].keys():
                if key == "order":
                    pass
                else:
                    for key2 in param_from_file["labels"][label_selected][key]:
                        # checks that key2 is in common
                        if key in param_from_file["common"].keys():
                            param_from_file["common"][key][key2] = param_from_file[
                                "labels"
                            ][label_selected][key][key2]
                            number_keys_ammended += 1
                        else:
                            print_log(
                                f"Did not find key <{key}> in common dictionary",
                                status="WARN",
                            )
            print_log(f"Amended {number_keys_ammended} keys for {label_selected}")

        # need to replace default keys by those in 'label' key
        converted_param["acquisition"]["label"] = label_selected
        self.param_dict = converted_param

    def set_channel(self, key, default):
        if key in self.param_dict["acquisition"].keys():
            channel = self.param_dict["acquisition"][key]
        else:
            channel = default

        return channel

    # method returns label-specific filenames from filename list
    def find_files_to_process(self, files_folder):
        # defines channel for DAPI, fiducials and barcodes
        channel_dapi = self.set_channel("DAPI_channel", "ch00")
        channel_barcode = self.set_channel("barcode_channel", "ch01")
        channel_mask = self.set_channel("mask_channel", "ch01")
        channel_barcode_fiducial = self.set_channel("fiducialBarcode_channel", "ch00")
        channel_mask_fiducial = self.set_channel("fiducialMask_channel", "ch00")

        # finds if there are 2 or 3 channels for DAPI acquisition
        dapi_files = [
            file
            for file in files_folder
            if self.decode_file_parts(path.basename(file))["channel"] == "ch02"
            and "DAPI" in path.basename(file).split("_")
        ]

        # defines channels for RNA and DAPI-fiducial
        if len(dapi_files) > 0:
            channel_dapi_fiducial = self.set_channel("fiducialDAPI_channel", "ch02")
            channel_dapi_rna = self.set_channel("RNA_channel", "ch01")
        else:
            channel_dapi_fiducial = self.set_channel("fiducialDAPI_channel", "ch01")
            channel_dapi_rna = self.set_channel("RNA_channel", "ch04")

        if channel_dapi_fiducial and len(dapi_files) == 0:
            warn(
                "\n\n****You are using ch02 for channel_dapi_fiducial but there are only 2 channels for DAPI!\n\n"
            )

        # selects DAPI files
        if self.param_dict["acquisition"]["label"] == "DAPI":
            self.files_to_process = [
                file
                for file in files_folder
                if self.decode_file_parts(path.basename(file))["channel"]
                == channel_dapi
                and "DAPI" in path.basename(file).split("_")
            ]

        # selects DAPIch2 files
        elif self.param_dict["acquisition"]["label"] == "RNA":
            self.files_to_process = [
                file
                for file in files_folder
                if self.decode_file_parts(path.basename(file))["channel"]
                == channel_dapi_rna
                and "DAPI" in path.basename(file).split("_")
            ]

        # selects barcode files
        elif self.param_dict["acquisition"]["label"] == "barcode":
            self.files_to_process = [
                file
                for file in files_folder
                if len([i for i in file.split("_") if "RT" in i]) > 0
                and self.decode_file_parts(path.basename(file))["channel"]
                == channel_barcode
            ]

        # selects mask files
        elif self.param_dict["acquisition"]["label"] == "mask":
            self.files_to_process = [
                file
                for file in files_folder
                if len([i for i in file.split("_") if "mask" in i]) > 0
                and self.decode_file_parts(path.basename(file))["channel"]
                == channel_mask
            ]

        # selects fiducial files
        elif self.param_dict["acquisition"]["label"] == "fiducial":
            self.files_to_process = [
                file
                for file in files_folder
                if (
                    len([i for i in file.split("_") if "RT" in i]) > 0
                    and self.decode_file_parts(path.basename(file))["channel"]
                    == channel_barcode_fiducial
                )
                or (
                    len([i for i in file.split("_") if "mask" in i]) > 0
                    and self.decode_file_parts(path.basename(file))["channel"]
                    == channel_mask_fiducial
                )
                or (
                    "DAPI" in file.split("_")
                    and self.decode_file_parts(path.basename(file))["channel"]
                    == channel_dapi_fiducial
                )
            ]

        else:
            self.files_to_process = []

        print_log(f"$ Files to process: {len(self.files_to_process)}")
        for i, file in enumerate(self.files_to_process):
            print_log(f"{i}\t{os.path.basename(file)}")

    def decode_file_parts(self, file_name):
        """
        decodes variables from an input file. typically, RE takes the form:

        "scan_(?P<runNumber>[0-9]+)_(?P<cycle>[\\w|-]+)_(?P<roi>[0-9]+)_ROI_converted_decon_(?P<channel>[\\w|-]+).tif" # pylint: disable=anomalous-backslash-in-string

        thus, by running decode_file_parts(current_param,file_name) you will get back either an empty dict if the RE were not present
        in your infoList...json file or a dict as follows if it all worked out fine:

        file_parts['runNumber']: runNumber number
        file_parts['cycle']: cycle string
        file_parts['roi']: roi number
        file_parts['channel']: channel string

        Parameters
        ----------
        file_name : string
            filename to decode

        Returns
        -------
        Dict with file_parts.

        """

        # decodes regular expressions
        if "fileNameRegExp" in self.param_dict["acquisition"].keys():
            file_parts = re.search(
                self.param_dict["acquisition"]["fileNameRegExp"], file_name
            )
            return file_parts
        else:
            return {}


class DaskCluster:
    def __init__(
        self, requested_number_nodes, maximum_load=0.6, memory_per_worker=12000
    ):
        self.requested_number_nodes = requested_number_nodes
        # self.n_threads will be created after exetution of initialize_cluster()
        self.maximum_load = maximum_load  # max number of workers that I can take
        self.memory_per_worker = memory_per_worker  # in Mb
        self.initialize_cluster()
        self.cluster = None
        self.client = None

    def initialize_cluster(self):
        number_cores_available = multiprocessing.cpu_count()

        # we want at least 12 GB per worker
        _, _, free_m = map(int, os.popen("free -t -m").readlines()[-1].split()[1:])

        max_number_threads = int(
            np.min(
                [
                    number_cores_available * self.maximum_load,
                    free_m / self.memory_per_worker,
                ]
            )
        )

        self.n_threads = int(np.min([max_number_threads, self.requested_number_nodes]))

        print(
            f"$ Cluster with {self.n_threads} workers started ({self.requested_number_nodes} requested)"
        )

    def create_distributed_client(self):
        client = try_get_client()
        if client is not None:
            print("# Shutting down existing cluster! ")
            client.shutdown()
        else:
            print("$ No running cluster detected. Will start one.")

        self.cluster = LocalCluster(
            n_workers=self.n_threads,
            threads_per_worker=1,
            memory_limit="64GB",
        )
        self.client = Client(self.cluster)

        print("$ Go to http://localhost:8787/status for information on progress...")


# =============================================================================
# FUNCTIONS
# =============================================================================


def load_parameters_file(file_name):
    if path.exists(file_name):
        with open(file_name, encoding="utf-8") as json_file:
            parameters = json.load(json_file)

        print(f"$ Parameters file read: {file_name}")
        return parameters
    else:
        return None


def save_json(file_name, data):
    with open(file_name, mode="w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, sort_keys=True, indent=4)


def load_json(file_name):
    if path.exists(file_name):
        with open(file_name, encoding="utf-8") as json_file:
            data = json.load(json_file)
    else:
        data = {}
    return data


def is_notebook():
    """
    This function detects if you are running on an ipython console or in the shell.
    It is used to either kill plots or leave them open.

    Returns
    -------
    TYPE Boolean
        true if running in Jupyter or Ipython consoles within spyder.
        false otherwise (terminal)

    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True  # Jupyter notebook or qtconsole
        elif shell == "TerminalInteractiveShell":
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False  # Probably standard Python interpreter


def rt_to_filename(current_param, reference_barcode):
    """
    Finds the files in a list that contain the ReferenceBarcode in their name
    Also returs the ROI of each file in this list


    Parameters
    ----------
    current_param : class
        parameters class.
    reference_barcode : string
        reference barcode name

    Returns
    -------
    filenames_with_ref_barcode : list
        list of files with reference barcode in their name
    roi_list : list
        list of rois.

    """
    filenames_with_ref_barcode = []
    roi_list = {}

    for file in current_param.files_to_process:
        if reference_barcode in file.split("_"):
            filenames_with_ref_barcode.append(file)
            # roi_list[file] = os.path.basename(file).split("_")[position_roi_information]
            file_parts = current_param.decode_file_parts(os.path.basename(file))
            roi_list[file] = file_parts["roi"]
    return filenames_with_ref_barcode, roi_list


def roi_to_fiducial_filename(current_param, file, barcode_name):
    """
    Produces list of fiducial files that need to be loaded from a specific mask/barcode image


    Parameters
    ----------
    current_param : class
        parameters class.
    file : string
        filename.
    barcode_name : string
        barcode name to assign a fiducial to.

    Returns
    -------
    candidates : TYPE
        DESCRIPTION.

    """
    # gets root_folder
    root_folder = os.path.dirname(file)
    roi = current_param.decode_file_parts(os.path.basename(file))["roi"]

    channel_fiducial = current_param.param_dict["acquisition"][
        "fiducialBarcode_channel"
    ]

    # looks for referenceFiducial file in folder
    list_files = glob.glob(root_folder + os.sep + "*.tif")

    candidates = [
        x
        for x in list_files
        if (barcode_name + "_" in x)
        and (roi == current_param.decode_file_parts(os.path.basename(x))["roi"])
        and (channel_fiducial in os.path.basename(x))
    ]

    return candidates


def load_alignment_dictionary(data_folder):
    dict_filename = (
        os.path.splitext(data_folder.output_files["dictShifts"])[0] + ".json"
    )

    dict_shifts = load_json(dict_filename)
    if len(dict_shifts) == 0:
        print_log(f"File with dictionary not found!: {dict_filename}")
        dict_shifts_available = False
    else:
        print_log(f"Dictionary File loaded: {dict_filename}")
        dict_shifts_available = True

    return dict_shifts, dict_shifts_available


def try_get_client():
    try:
        client = get_client()
        client.restart()
    except ValueError:
        client = None

    return client


def restart_client():
    # restarts client
    client = try_get_client()
    if client is not None:
        client.restart()
        print("$ Distributed network restarted")


def print_dict(dictionary):
    print("\n$ Parameters loaded:")
    for key in dictionary.keys():
        spacer = "\t" * (3 - int(len(key) / 8))
        print(f"\t{key}{spacer}{dictionary[key]}")
    print("\n")


def get_dictionary_value(dictionary, key, default=""):
    if key in dictionary.keys():
        value = dictionary[key]
    else:
        value = default

    return value


def loads_barcode_dict(file_name):
    import json

    # Check if the file exists
    if not os.path.exists(file_name):
        print("File does not exist")
        return dict()
    else:
        # Opening JSON file
        f = open(file_name)

        # returns JSON object as a dictionary
        barcode_type = json.load(f)

        # Closing file
        f.close()

        print("$ {} barcode dictionary loaded")
        return barcode_type
