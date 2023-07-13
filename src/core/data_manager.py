#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data manager module

Manage writing, reading and checking data.
"""

import json
import os
import re

import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from skimage import io

from core.pyhim_logging import Logger, print_log, write_string_to_file
from core.saving import image_show_with_values


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


class ImageFile:
    def __init__(self, path, img_name, ext, label):
        self.all_path = path
        self.name = img_name
        self.extension = ext
        self.root = self.get_root()
        # self.acquisition = ""  # DAPI, Barcode, Mask
        # self.channel = channel
        self.m_label = label
        # self.roi = roi
        # self.cycle = cycle

    def get_root(self):
        length = len(self.all_path) - len(self.name) - 1 - len(self.extension)
        return self.all_path[:length]

    def load(self):
        return io.imread(self.all_path).squeeze()

    # def save(self, result, folder_name: str, tag: str):
    #     path_name = self.root + os.sep + folder_name + os.sep + self.name + tag
    #     if result.shape > (1, 1):
    #         np.save(path_name, result)
    #         print_log(f"$ Image saved to disk: {path_name}.npy", "info")
    #     else:
    #         print_log("# Warning, image is empty", "Warning")


def remove_extension(filename: str):
    fl_split = filename.split(".")
    fl_split.pop()
    return ".".join(fl_split)


class DataManager:
    """Single party responsible for communicating data with the system"""

    def __init__(
        self,
        data_path: str,
        logger: Logger,
        stardist_basename: str = "",
        params_filename: str = "infoList",
    ):
        self.m_data_path = self.__set_data_path(data_path)
        self.out_path = self.m_data_path
        self.m_logger = logger
        self.m_stardist_basename = stardist_basename
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

    def create_folder(self, folder_name: str):
        """Create folder with `makedirs` from os module.
        It's a recursive directory creation function.

        Parameters
        ----------
        folder_name : str
            Relative path name of folder
        """
        folder_path = self.out_path + os.sep + folder_name
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Folder '{folder_path}' created successfully.")
        else:
            print(f"Folder '{folder_path}' already exists.")

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
        img_ext = ["tif", "tiff"]
        # img_ext = ["tif", "tiff", "npy", "png", "jpg"]
        table_ext = ["csv", "ecsv", "dat"]
        for path, name, ext in self.all_files:
            if ext in img_ext:
                label = self.find_label(name)
                self.data_images.append(ImageFile(path, name, ext, label))
            elif ext in table_ext:
                self.data_tables.append((path, name, ext))
            elif ext in ["log", "md"] or (
                ext == "json" and name == self.m_filename_params
            ):
                pass
            else:
                print(f"Unrecognized data file: {path}")

    def find_label(self, filename):
        parts = self.decode_file_parts(filename)
        channel = parts["channel"]

        if "DAPI" in filename.split("_"):
            label = self.channels["dapi_acq"][channel]
        elif "RT" in filename:
            label = self.channels["barcode_acq"][channel]
        elif "mask" in filename.split("_"):
            label = self.channels["mask_acq"][channel]
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
        params = load_json(self.param_file_path)
        if params is None:
            raise ValueError(f"Parameters file NOT FOUND: {self.param_file_path}")
        print(f"$ Parameters file read: {self.param_file_path}")
        return params

    def set_up(self, acquisition_dict: dict):
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

        self.dispatch_files()

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

    def get_inputs(self, labels: list[str]):
        inputs = []
        for img_file in self.data_images:
            if img_file.m_label in labels:
                inputs.append(img_file)
        return inputs

    def save_data(
        self,
        results: list,
        tags: list[str],
        out_folder: str,
        associated_file: ImageFile,
    ):
        if len(results) != len(tags):
            # TODO: Improuve this possibility, you may be can have just one image but different tags
            # Because we want plot this image with different ways
            raise SystemExit(f"len(results) {len(results)} != len(tags) {len(tags)}")

        partial_path = (
            self.out_path + os.sep + out_folder + os.sep + associated_file.name
        )
        for res, tag in zip(results, tags):
            if tag == "_2d":
                self._save_2d_npy(res, partial_path)
                self._save_2d_png(res, partial_path)
            elif tag == "_focalPlaneMatrix":
                self._save_focal_plane_matrix((res, partial_path))
            else:
                raise SystemExit(f"tag UNRECOGNIZED: {tag}")

    def _save_2d_npy(self, data, partial_path):
        path_name = partial_path + "_2d"
        if data.shape > (1, 1):
            np.save(path_name, data)
            print_log(f"$ Image saved to disk: {path_name}.npy", "info")
        else:
            raise ValueError(f"Image is empty! Original file: {partial_path}.tif")

    def _save_2d_png(self, data, partial_path):
        out_path = partial_path + "_2d.png"

        fig = plt.figure()
        size = (10, 10)
        fig.set_size_inches(size)
        ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
        ax.set_axis_off()
        norm = ImageNormalize(stretch=SqrtStretch())

        ax.set_title("2D Data")

        fig.add_axes(ax)
        ax.imshow(data, origin="lower", cmap="Greys_r", norm=norm)
        fig.savefig(out_path)
        plt.close(fig)
        original_filename = os.path.basename(partial_path + ".tif")
        write_string_to_file(
            self.m_logger.md_filename,
            f"{original_filename}\n ![]({out_path})\n",
            "a",
        )

    def _save_focal_plane_matrix(self, data: tuple, partial_path: str):
        focal_plane_matrix, focus_plane = data
        out_path = partial_path + "_focalPlaneMatrix.png"
        image_show_with_values(
            [focal_plane_matrix],
            title="focal plane = " + f"{focus_plane:.2f}",
            output_name=out_path,
        )
        original_filename = os.path.basename(partial_path + ".tif")
        write_string_to_file(
            self.m_logger.md_filename,
            f"{original_filename}\n ![]({out_path})\n",
            "a",
        )


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
