#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data manager module

Manage writing, reading and checking data.
"""

import json
import os

from skimage import io


def extract_files(root: str):
    """Extract recursively file informations of all files into a given directory.
    Note: filename is the name without extension

    Parameters
    ----------
    root : str
        The name of root directory

    Returns
    -------
    List[Tuple(str,str,str)]
        List of file informations: (filepath, filename, extension)
    """
    files = []
    # Iterate into dirname and each subdirectories dirpath, dirnames, filenames
    for dirpath, dirnames, filenames in os.walk(root):
        for filename in filenames:
            split_filename = filename.split(".")
            extension = split_filename.pop()
            short_filename = ".".join(split_filename)
            filepath = os.path.join(dirpath, filename)
            files.append((filepath, short_filename, extension))

        if len(dirnames) > 0:
            print("[INFO] Inside:")
            print(dirpath)
            print("[INFO] Subdirectories detected:")
            print(dirnames)

    return files


class DataManager:
    """Single party responsible for communicating data with the system"""

    def __init__(
        self,
        data_path: str,
        stardist_basename: str = "",
        params_filename: str = "infoList",
    ):
        self.m_data_path = self.__set_data_path(data_path)
        self.m_stardist_basename = str(stardist_basename)
        self.m_filename_params = params_filename
        self.all_files = extract_files(self.m_data_path)
        self.param_file_path = self.find_param_file(params_filename)
        self.data_images = []
        self.data_tables = []
        self.filename_regex = ""
        self.channels = {
            "dapi_acq": {
                "ch00": "dapi",
                "ch01": "fiducial",
                "ch02": "rna",
            },
            "mask_acq": {
                "ch00": "mask",
                "ch01": "fiducial",
            },
            "barcode_acq": {
                "ch00": "barcode",
                "ch01": "fiducial",
            },
        }
        self.img_info = {
            "parallelize_planes": False,
            "pixel_size_XY": 0.1,
            "pixel_size_Z": 0.25,
            "z_binning": 2,
        }
        # self.dispatch_files()

    @staticmethod
    def __set_data_path(data_path):
        if data_path:
            return str(data_path)
        return os.getcwd()

    def find_param_file(self, params_filename):
        for path, name, ext in self.all_files:
            if ext == "json" and name == self.m_filename_params:
                return str(path)
        # If we loop over all files, parameter file aren't detected.
        raise ValueError(
            f"Parameters file NOT FOUND, expected filename: {params_filename}.json"
        )

    def dispatch_files(self):
        """Get all input files and sort by extension type"""
        img_ext = ["tif", "tiff", "npy", "png", "jpg"]
        table_ext = ["csv", "ecsv", "dat"]
        for path, name, ext in self.all_files:
            if ext in img_ext:
                label = self.find_label(name)
                self.data_images.append(
                    ImageFile(path, parts["cycle"], parts["roi"], parts["channel"])
                )
            elif ext in table_ext:
                self.data_tables.append((path, name, ext))
            # elif ext == "json" and name == self.m_filename_params:
            #     self.user_parameter = str(path)
            else:
                print(f"Unrecognized data file: {path}")

    def find_label(self, filename):
        parts = self.decode_file_parts(filename)
        channel = parts["channel"]

        if "DAPI" in path.basename(file).split("_"):
            return self.channels["dapi_acq"][channel]
        elif "RT" in path.basename(file).split("_"):
            return self.channels["barcode_acq"][channel]
        elif "mask" in path.basename(file).split("_"):
            return self.channels["mask_acq"][channel]
        else:
            raise ValueError(f"Label NOT FOUND for this filename: {filename}")

        return label

    def load_user_param(self):
        """Load user parameter JSON file like a Python dict

        Returns
        -------
        dict
            Python dict

        Raises
        ------
        ValueError
            file not found
        """
        params = load_json(self.user_parameter)
        if params is None:
            raise ValueError(f"Parameters file NOT FOUND: {self.user_parameter}")
        print(f"$ Parameters file read: {self.user_parameter}")
        return params

    def set_up(acquisition_dict: dict):
        acq = acquisition_dict["common"]
        # Regular expression
        self.filename_regex = remove_extension(acq["fileNameRegExp"])
        # Channels
        self.channels["dapi_acq"][acq.get("DAPI_channel", "ch00")] = "dapi"
        self.channels["dapi_acq"][acq.get("fiducialDAPI_channel", "ch01")] = "fiducial"
        self.channels["dapi_acq"][acq.get("RNA_channel", "ch02")] = "rna"
        self.channels["mask_acq"][acq.get("mask_channel", "ch00")] = "mask"
        self.channels["mask_acq"][acq.get("fiducialMask_channel", "ch01")] = "fiducial"
        self.channels["barcode_acq"][acq.get("barcode_channel", "ch00")] = "barcode"
        self.channels["barcode_acq"][
            acq.get("fiducialBarcode_channel", "ch01")
        ] = "fiducial"
        # Image informations
        self.img_info["parallelize_planes"] = acq["parallelizePlanes"]
        self.img_info["pixel_size_XY"] = acq["pixelSizeXY"]
        self.img_info["pixel_size_Z"] = acq["pixelSizeZ"]
        self.img_info["z_binning"] = acq["zBinning"]

    def decode_file_parts(self, file_name):
        """
        decodes variables from an input file. typically, RE takes the form:

        "scan_(?P<runNumber>[0-9]+)_(?P<cycle>[\\w|-]+)_(?P<roi>[0-9]+)_ROI_converted_decon_(?P<channel>[\\w|-]+)" # pylint: disable=anomalous-backslash-in-string,line-too-long

        thus, by running decode_file_parts(file_name) you will get back
        either an empty dict if the RE were not present
        in your infoList.json file or a dict as follows if it all worked out fine:

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
        if self.filename_regex:
            return re.search(self.filename_regex, file_name)
        return None

    def get_inputs(self, labels: List[str]):
        inputs = []
        for img_file in self.data_images:
            if img_file.label in labels:
                inputs.append(img_file)
        return inputs


def load_json(file_name):
    """Load a JSON file like a python dict

    Parameters
    ----------
    file_name : str
        JSON file name

    Returns
    -------
    dict
        Python dict
    """
    if os.path.exists(file_name):
        with open(file_name, encoding="utf-8") as json_file:
            return json.load(json_file)
    return None


def write_string_to_file(file_name, text_to_output, attribute="a"):
    """write a line of text into a file

    Parameters
    ----------
    file_name : str
        log file
    text_to_output : str
        text to write in file
    attribute : str, optional
        Open file mode option, by default "a"
    """
    with open(file_name, mode=attribute, encoding="utf-8") as file_handle:
        file_handle.write(str(text_to_output) + "\n")


def save_json(data, file_name):
    """Save a python dict as a JSON file

    Parameters
    ----------
    data : dict
        Data to save
    file_name : str
        Output JSON file name
    """
    with open(file_name, mode="w", encoding="utf-8") as json_f:
        json.dump(data, json_f, ensure_ascii=False, sort_keys=True, indent=4)


# class File:
#     def __init__(self, path, name, ext):
#         self.path = ""
#         self.type = []  # Fiducial, barcode, mask, dapi, rna, localization, trace, ...
#         self.status = []  # projected, registered, filtered
#         self.reference_to = ""  # path ?

#     def load(self):
#         pass

#     def save(self):
#         pass


class ImageFile(File):
    def __init__(self, path, label):
        self.m_path = path
        # self.acquisition = ""  # DAPI, Barcode, Mask
        # self.channel = channel
        self.m_label = label
        # self.roi = roi
        # self.cycle = cycle

    def load(self):
        io.imread(self.m_path).squeeze()


def remove_extension(filename: str):
    return ".".join(filename.split(".").pop())
