#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and functions for file management
"""

import json
import os
import re
from os import path
from warnings import warn

from core.data_manager import load_json
from core.pyhim_logging import print_log


class Parameters:
    """Manage all pyHiM parameters."""

    def __init__(
        self,
        raw_dict,
        root_folder="./",
        label="",
        stardist_basename=None,
    ):
        # self.label = label
        self.files_to_process = []
        # self.param_dict = self.get_standard_parameters()
        self.param_dict = self.complete_with_default(raw_dict)
        if label:
            self.param_dict = self.get_labelled_params(label)
        # self.convert_parameter_file(raw_dict)
        self.set_stardist_basename(stardist_basename)
        self.param_dict["rootFolder"] = root_folder
        self.file_parts = {}

    def complete_with_default(self, raw_dict):
        default = self.get_standard_parameters()
        return deep_dict_update(default, raw_dict)

    def set_stardist_basename(self, stardist_basename: str):
        if stardist_basename is not None:
            if self.param_dict.get("common", False):
                self.param_dict["common"]["segmentedObjects"][
                    "stardist_basename"
                ] = stardist_basename
            else:
                self.param_dict["segmentedObjects"][
                    "stardist_basename"
                ] = stardist_basename

    # def convert_parameter_file(self, raw_dict):
    #     """Load and organize parameters from a JSON file

    #     Raises
    #     ------
    #     ValueError
    #         Error if file not found
    #     """
    #     param_from_file = raw_dict

    #     converted_param = param_from_file["common"]
    #     labels = param_from_file["labels"]

    #     ordered_list = [" "] * len(labels.keys())
    #     for _, label in enumerate(labels.keys()):
    #         order = labels[label]["order"]
    #         ordered_list[order - 1] = label

    #     converted_param["labels"] = ordered_list

    #     label_selected = self.label  # the label specific parameters that we want load

    #     # need to add keys not present in common dict
    #     if len(label_selected) > 0:
    #         number_keys_ammended = 0
    #         for key in param_from_file["labels"][label_selected].keys():
    #             if key == "order":
    #                 pass
    #             else:
    #                 for key2 in param_from_file["labels"][label_selected][key]:
    #                     # checks that key2 is in common
    #                     if key in param_from_file["common"].keys():
    #                         param_from_file["common"][key][key2] = param_from_file[
    #                             "labels"
    #                         ][label_selected][key][key2]
    #                         number_keys_ammended += 1
    #                     else:
    #                         print_log(
    #                             f"Did not find key <{key}> in common dictionary",
    #                             status="WARN",
    #                         )
    #         print_log(f"Amended {number_keys_ammended} keys for {label_selected}")

    #     # need to replace default keys by those in 'label' key
    #     converted_param["acquisition"]["label"] = label_selected
    #     self.param_dict = converted_param

    def get_sectioned_params(self, section_name: str):
        tempo = self.param_dict["common"].get(section_name)
        section_dict = {"common": tempo, "labels": {}}
        for key, value in self.param_dict["labels"].items():
            section_dict["labels"][key] = value.get(section_name)
        return section_dict

    def get_labelled_params(self, label_name: str):
        main_dict = self.param_dict["common"]
        main_dict["acquisition"]["label"] = label_name
        update = self.param_dict["labels"].get(label_name, {})
        return deep_dict_update(main_dict, update)

    def get_labeled_dict_value(self, section, param_name):
        labeled_dict = {}
        default = self.param_dict["common"][section][param_name]
        for key, value in self.param_dict["labels"].items():
            key_lower = key.lower()
            labeled_dict[key_lower] = value.get(section, {}).get(param_name, default)
        return labeled_dict

    def set_channel(self, key, default):
        """Set channel parameter with a default value

        Parameters
        ----------
        key : str
            Name like DAPI_channel, barcode_channel, fiducialMask_channel, ...
        default : str
            Like ch00, ch01, ...

        Returns
        -------
        str
            Channel value like 'ch02'
        """

        if key in self.param_dict["acquisition"].keys():
            return self.param_dict["acquisition"][key]
        return default

    # method returns label-specific filenames from filename list
    def find_files_to_process(self, files_folder):
        """Find label-specific filenames from filename list.
        Save these filenames in self.files_to_process.

        Parameters
        ----------
        files_folder : list
            List of files
        """
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
        if dapi_files:
            channel_dapi_fiducial = self.set_channel("fiducialDAPI_channel", "ch02")
            channel_dapi_rna = self.set_channel("RNA_channel", "ch01")
        else:
            channel_dapi_fiducial = self.set_channel("fiducialDAPI_channel", "ch01")
            channel_dapi_rna = self.set_channel("RNA_channel", "ch04")

        if channel_dapi_fiducial and not dapi_files:
            warn(
                "\n\n****You are using ch02 for channel_dapi_fiducial \
                    but there are only 2 channels for DAPI!\n\n"
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
        # sourcery skip: use-named-expression
        """
        decodes variables from an input file. typically, RE takes the form:

        "scan_(?P<runNumber>[0-9]+)_(?P<cycle>[\\w|-]+)_(?P<roi>[0-9]+)_ROI_converted_decon_(?P<channel>[\\w|-]+).tif" # pylint: disable=anomalous-backslash-in-string,line-too-long

        thus, by running decode_file_parts(current_param,file_name) you will get back
        either an empty dict if the RE were not present
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
        file_parts = {}
        # decodes regular expressions
        regex = self.param_dict.get("acquisition").get("fileNameRegExp")
        return re.search(regex, file_name) if regex else None

    @staticmethod
    def get_standard_parameters():
        """Reference of the standard parameters"""
        return {
            "common": {
                "acquisition": {
                    "positionROIinformation": 3,
                    "fileNameRegExp": "scan_(?P<runNumber>[0-9]+)_(?P<cycle>[\\w|-]+)_\
                        (?P<roi>[0-9]+)_ROI_converted_decon_(?P<channel>[\\w|-]+).tif",
                    "DAPI_channel": "ch00",
                    "fiducialDAPI_channel": "ch01",
                    "RNA_channel": "ch02",
                    "fiducialBarcode_channel": "ch00",
                    "fiducialMask_channel": "ch00",
                    "barcode_channel": "ch01",
                    "mask_channel": "ch01",
                    # in future this field will contain the ch for the label.
                    # This parameter will supersed the individual channel fields above.
                    "label_channel": "ch00",
                    # in future this field will contain the ch for the label fiducial.
                    # This parameter will supersed the individual channel fields above.
                    "label_channel_fiducial": "ch01",
                    "pixelSizeXY": 0.1,
                    "zBinning": 2,
                    # if True it will parallelize inner loops (plane by plane).
                    # Otherwise outer loops (e.g. file by file)
                    "parallelizePlanes": False,
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
                    # Used in blockAlignment to determine the % of error tolerated
                    "tolerance": 0.1,
                    # lower threshold to adjust image intensity levels
                    # before xcorrelation for alignment in 2D
                    "lower_threshold": 0.999,
                    # higher threshold to adjust image intensity levels
                    # before xcorrelation for alignment in 2D
                    "higher_threshold": 0.9999999,
                    # lower threshold to adjust image intensity levels
                    # before xcorrelation for Alignment3D
                    "3D_lower_threshold": 0.9,
                    # higher threshold to adjust image intensity levels
                    # before xcorrelation for Alignment3D
                    "3D_higher_threshold": 0.9999,
                    "background_sigma": 3.0,  # used to remove inhom background
                    "localShiftTolerance": 1,
                    "blockSize": 256,
                },
                "buildsPWDmatrix": {
                    "folder": "buildsPWDmatrix",  # output folder
                    # available methods: masking, clustering
                    "tracing_method": ["masking", "clustering"],
                    # Expands masks until they collide by a max of 'mask_expansion' pixels
                    "mask_expansion": 8,
                    "flux_min": 10,  # min flux to keeep object
                    "flux_min_3D": 0.1,  # min flux to keeep object
                    "kd_tree_distance_threshold_mum": 1,  # distance threshold used to build KDtree
                    # colormaps used for plotting matrices
                    "colormaps": {
                        "PWD_KDE": "terrain",
                        "PWD_median": "terrain",
                        "contact": "coolwarm",
                        "Nmatrix": "Blues",
                    },
                    # zxy tolerance used for block drift correction, in px
                    "toleranceDrift": [3, 1, 1],
                    # if True it will removed uncorrected localizations,
                    # otherwise they will remain uncorrectd.
                    "remove_uncorrected_localizations": True,
                },
                "segmentedObjects": {
                    "folder": "segmentedObjects",  # output folder
                    "operation": "2D,3D",  # options: 2D or 3D
                    "outputFile": "segmentedObjects",
                    "Segment3D": "overwrite",
                    "background_method": "inhomogeneous",  # flat or inhomogeneous or stardist
                    "stardist_basename": "/mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks",
                    # network for 2D barcode segmentation
                    "stardist_network": "stardist_nc14_nrays:64_epochs:40_grid:2",
                    # network for 3D barcode segmentation
                    "stardist_network3D": "stardist_nc14_nrays:64_epochs:40_grid:2",
                    "tesselation": True,  # tesselates masks
                    "background_sigma": 3.0,  # used to remove inhom background
                    "threshold_over_std": 1.0,  # threshold used to detect sources
                    "fwhm": 3.0,  # source size in px
                    "brightest": 1100,  # max number of sources segmented per FOV
                    "intensity_min": 0,  # min int to keep object
                    "intensity_max": 59,  # max int to keeep object
                    "area_min": 50,  # min area to keeep object
                    "area_max": 500,  # max area to keeep object
                    # options: 'thresholding' or 'stardist', 'zASTROPY', 'zProfile'
                    "3Dmethod": "thresholding",
                    "residual_max": 2.5,  # z-profile Fit: max residuals to keeep object
                    "sigma_max": 5,  # z-profile Fit: max sigma 3D fitting to keeep object
                    # z-profile Fit: max diff between Moment and z-gaussian fits to keeep object
                    "centroidDifference_max": 5,
                    # z-profile Fit: window size to extract subVolume, px.
                    # 3 means subvolume will be 7x7.
                    "3DGaussianfitWindow": 3,
                    # constructs a YZ image by summing from xPlane-window:xPlane+window
                    "3dAP_window": 5,
                    "3dAP_flux_min": 2,  # # threshold to keep a source detected in YZ
                    "3dAP_brightest": 100,  # number of sources sought in each YZ plane
                    # px dist to attribute a source localized in YZ to one localized in XY
                    "3dAP_distTolerance": 1,
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
                    # if reducePlanes==True it will calculate focal plane and only use a region
                    # around it for segmentSources3D, otherwise will use the full stack
                    "reducePlanes": True,
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


def load_alignment_dict(data_folder):
    """Load a JSON file with 'dictShifts' in the file name.

    Parameters
    ----------
    data_folder : Folders
        Folders object

    Returns
    -------
    (dict,dict)
        Shift dictionaries
    """
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


def print_dict(dictionary):
    """Print parameters in your shell terminal

    Parameters
    ----------
    dictionary : dict
        Parameters dictionary
    """
    print("\n$ Parameters loaded:")
    for key in dictionary.keys():
        spacer = "\t" * (3 - len(key) // 8)
        print(f"\t{key}{spacer}{dictionary[key]}")
    print("\n")


def get_dictionary_value(dictionary, key, default=""):
    """Get dict value with a default option if key doesn't exist.

    Parameters
    ----------
    dictionary : dict
        dictionary object
    key : str
        key for the dict
    default : str, optional
        default value is key doesn't exist, by default ""

    Returns
    -------
    str
        value or default
    """
    return dictionary[key] if key in dictionary.keys() else default


def loads_barcode_dict(file_name):
    """Loads a barcode type dict JSON

    Parameters
    ----------
    file_name : str
        JSON file name

    Returns
    -------
    dict
        Barcode dictionary
    """
    bc_dict = {}
    # Check if the file exists
    if not os.path.exists(file_name):
        print("File does not exist")
    else:
        # Opening JSON file
        with open(file_name, encoding="utf-8") as json_f:
            # returns JSON object as a dictionary
            barcode_type = json.load(json_f)
            print("$ {} barcode dictionary loaded")
            bc_dict = barcode_type

    return bc_dict


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
            file_parts = current_param.decode_file_parts(os.path.basename(file))
            roi_list[file] = file_parts["roi"]
    return filenames_with_ref_barcode, roi_list


def deep_dict_update(main_dict: dict, update: dict):
    for key, value in update.items():
        if isinstance(value, dict):
            main_dict[key] = deep_dict_update(main_dict.get(key, {}), value)
        else:
            main_dict[key] = value
    return main_dict
