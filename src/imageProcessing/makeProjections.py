#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 23:17:58 2020

@author: marcnol

This file contains functions to project 3D images to 2D

Operation will be defined in the parameters file. Options are:
    - user-defined range
    - all z range
    - optimal range based on detection of focal plane and use of user defined window around it


"""
# =============================================================================
# IMPORTS
# =============================================================================

import glob
import os
import warnings

import numpy as np
from apifish.stack import projection
from dask.distributed import get_client, wait

from core.data_manager import write_string_to_file
from core.folder import Folders
from core.pyhim_logging import print_log, print_session_name
from core.saving import image_show_with_values
from imageProcessing.imageProcessing import Image

warnings.filterwarnings("ignore")

np.seterr(divide="ignore", invalid="ignore")


class Feature:
    def __init__(self, param):
        self.m_param = param
        self.required_data = []
        self.required_ref = []
        self.required_table = []
        self.out_folder = ""

    def required_inputs(self):
        return self.required_data, self.required_ref, self.required_table


class Project(Feature):
    def __init__(self, param):
        super().__init__(param)
        self.required_data = ["barcode", "mask", "dapi", "mask", "rna"]
        self.out_folder = "zProject"

        self.block_size = 256
        self.display = True
        self.folder_name = "zProject"
        self.mode = "full"
        self.operation = "skip"
        self.save_image = True
        self.window_security = 2
        self.z_project_option = "MIP"
        self.zmax = 59
        self.zmin = 1
        self.zwindows = 15

    def run(self, img, option):
        # find the correct range for the projection
        img_reduce = self.precise_z_planes(img, option)
        self.print_image_properties()
        img_projected = self.projection_2d(img_reduce, option)
        return img_projected

    def check_zmax(self, img_size):
        if self.zmax > img_size[0]:
            print_log("$ Setting z max to the last plane")
            self.zmax = img_size[0]

    def precise_z_planes(self, img, option):
        img_size = img.shape
        self.check_zmax(img_size)
        return img[z_range[1][0] : z_range[1][-1]]

    def projection_2d(self, img, mode):
        # sums images
        i_collapsed = None

        if "MIP" in mode:
            # Max projection of selected planes
            i_collapsed = projection.maximum_projection(img)
        elif "sum" in mode:
            # Sums selected planes
            i_collapsed = projection.sum_projection(img)
        else:
            print_log(f"ERROR: mode not recognized. Expected: MIP or sum. Read: {mode}")

        return i_collapsed

    def print_image_properties(self):
        pass

    # processes sum image in axial direction given range
    # @jit(nopython=True)
    def z_projection_range(self):
        # find the correct range for the projection
        # if self.current_param.param_dict["zProject"]["zmax"] > self.image_size[0]:
        #     print_log("$ Setting z max to the last plane")
        #     self.current_param.param_dict["zProject"]["zmax"] = self.image_size[0]

        if self.current_param.param_dict["zProject"]["mode"] == "automatic":
            print_log("> Calculating planes...")
            z_range = calculate_zrange(self.data, self.current_param)

        elif self.current_param.param_dict["zProject"]["mode"] == "full":
            (zmin, zmax) = (0, self.image_size[0])
            z_range = (round((zmin + zmax) / 2), range(zmin, zmax))

        if self.current_param.param_dict["zProject"]["mode"] == "laplacian":
            print_log("Stacking using Laplacian variance...")
            (
                self.data_2d,
                self.focal_plane_matrix,
                self.z_range,
            ) = reinterpolate_focal_plane(self.data, self.current_param.param_dict)
            self.focus_plane = self.z_range[0]

        else:
            # Manual: reads from parameters file
            (zmin, zmax) = (
                self.current_param.param_dict["zProject"]["zmin"],
                self.current_param.param_dict["zProject"]["zmax"],
            )
            if zmin >= zmax:
                raise SystemExit(
                    "zmin is equal or larger than zmax in configuration file. Cannot proceed."
                )
            z_range = (round((zmin + zmax) / 2), range(zmin, zmax))

        if "laplacian" not in self.current_param.param_dict["zProject"]["mode"]:
            self.data_2d = project_image_2d(
                self.data,
                z_range,
                self.current_param.param_dict["zProject"]["zProjectOption"],
            )
            self.focus_plane = z_range[0]
            self.z_range = z_range[1]

        print_log(f"> Processing z_range:{self.z_range}")


# =============================================================================
# FUNCTIONS
# =============================================================================


def make_2d_projections_file(file_name, current_param, current_session, data_folder):
    if file_name in current_session.data:
        print_log(f"# File already projected: {os.path.basename(file_name)}")
    else:
        print_log(f"\n> Analysing file: {os.path.basename(file_name)}")

        # creates image object
        im_obj = Image(current_param)

        # loads image
        im_obj.load_image(file_name)

        # makes actual 2d projection
        im_obj.z_projection_range()

        # outputs information from file
        if im_obj.file_name:
            im_obj.print_image_properties()

        # saves output 2d zProjection as png
        if current_param.param_dict["zProject"]["display"]:
            png_file_name = (
                data_folder.output_folders["zProject"]
                + os.sep
                + os.path.basename(file_name)
                + "_2d.png"
            )
            im_obj.show_image(
                save=current_param.param_dict["zProject"]["saveImage"],
                output_name=png_file_name,
            )
            write_string_to_file(
                current_param.param_dict["fileNameMD"],
                f"{os.path.basename(file_name)}\n ![]({png_file_name})\n",
                "a",
            )  # initialises MD file

            if current_param.param_dict["zProject"]["mode"] == "laplacian":
                output_name = im_obj.get_image_filename(
                    data_folder.output_folders["zProject"], "_focalPlaneMatrix.png"
                )
                image_show_with_values(
                    [im_obj.focal_plane_matrix],
                    title="focal plane = " + f"{im_obj.focus_plane:.2f}",
                    output_name=output_name,
                )

                write_string_to_file(
                    current_param.param_dict["fileNameMD"],
                    f"{os.path.basename(file_name)}\n ![]({output_name})\n",
                    "a",
                )  # initialises MD file
        # saves output 2d zProjection as matrix
        im_obj.save_image_2d(data_folder.output_folders["zProject"])

        del im_obj


def make_projections(current_param, current_session, file_name=None):
    session_name = "makesProjections"

    # processes folders and files
    print_session_name(session_name)
    data_folder = Folders(current_param.param_dict["rootFolder"])
    print_log(f"> Folders read: {len(data_folder.list_folders)}")
    write_string_to_file(
        current_param.param_dict["fileNameMD"],
        f"""## {session_name}: {current_param.param_dict["acquisition"]["label"]}\n""",
        "a",
    )  # initialises MD file

    for current_folder in data_folder.list_folders:
        files_folder = glob.glob(current_folder + os.sep + "*.tif")
        data_folder.create_folders(current_folder, current_param)

        # generates lists of files to process
        current_param.find_files_to_process(files_folder)
        print_log(f"> Processing Folder: {current_folder}")
        print_log(f"> About to process {len(current_param.files_to_process)} files\n")

        if current_param.param_dict["parallel"]:
            threads = []
            files_to_process_filtered = [
                x
                for x in current_param.files_to_process
                if (file_name is None)
                or (
                    file_name is not None
                    and (
                        os.path.basename(x)
                        in [os.path.basename(x1) for x1 in file_name]
                    )
                )
            ]

            if len(files_to_process_filtered) > 0:
                # dask
                client = get_client()
                threads = [
                    client.submit(
                        make_2d_projections_file,
                        x,
                        current_param,
                        current_session,
                        data_folder,
                    )
                    for x in files_to_process_filtered
                ]

                print_log(f"$ Waiting for {len(threads)} threads to complete ")
                for _, _ in enumerate(threads):
                    wait(threads)

        else:
            for _, filename_to_process in enumerate(current_param.files_to_process):
                if (file_name is None) or (
                    file_name is not None
                    and (
                        os.path.basename(filename_to_process)
                        in [os.path.basename(x) for x in file_name]
                    )
                ):
                    make_2d_projections_file(
                        filename_to_process, current_param, current_session, data_folder
                    )
                    current_session.add(filename_to_process, session_name)
                else:
                    pass


# @jit(nopython=True)
def gaussian(x, a=1, mean=0, std=0.5):
    """Gaussian function

    Parameters
    ----------
    x : _type_
        _description_
    a : int, optional
        _description_, by default 1
    mean : int, optional
        _description_, by default 0
    std : float, optional
        _description_, by default 0.5

    Returns
    -------
    _type_
        _description_
    """
    return (
        a
        * (1 / (std * (np.sqrt(2 * np.pi))))
        * (np.exp(-((x - mean) ** 2) / ((2 * std) ** 2)))
    )


# =============================================================================
# FOCAL PLANE INTERPOLATION
# =============================================================================


def _remove_z_planes(image_3d, z_range):
    """
    Removes planes in input image.
    For instance, if you provide a z_range = range(0,image_3d.shape[0],2)

    then the routine will remove any other plane. Number of planes skipped
    can be programmed by tuning z_range.

    Parameters
    ----------
    image_3d : numpy array
        input 3D image.
    z_range : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """
    output = np.zeros((len(z_range), image_3d.shape[1], image_3d.shape[2]))
    for i, index in enumerate(z_range):
        output[i, :, :] = image_3d[index, :, :]

    return output


def _average_z_planes(image_3d, z_range):
    """
    Removes z-planes by calculating the average between successive planes

    Parameters
    ----------
    image_3d : numpy array
        input 3D image.
    z_range : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """

    output = np.zeros((len(z_range), image_3d.shape[1], image_3d.shape[2]))
    for i, index in enumerate(z_range):
        average = (
            image_3d[index, :, :].astype(np.float)
            + image_3d[index + 1, :, :].astype(np.float)
        ) / 2
        output[i, :, :] = average.astype(np.uint16)

    return output


def _interpolate_z_planes(image_3d, z_range):
    """
    Removes z planes by reinterpolation

    TODO

    Parameters
    ----------
    image_3d : numpy array
        input 3D image.
    z_range : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """

    output = np.zeros((len(z_range), image_3d.shape[1], image_3d.shape[2]))

    # need to code using interpn
    output = image_3d

    return output


def reinterpolate_z(image_3d, z_range, mode="average"):
    """
    wrapper function for any kind of z-interpolation
    to reduce the number of planes in an image

    Parameters
    ----------
    image_3d : numpy array
        input 3D image.
    z_range : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """
    if "interpolate" in mode:
        output = _interpolate_z_planes(image_3d, z_range)
    elif "remove" in mode:
        output = _remove_z_planes(image_3d, z_range)
    elif "average" in mode:
        output = _average_z_planes(image_3d, z_range)

    print_log(f"$ Reduced Z-planes from {image_3d.shape[0]} to {output.shape[0]}")

    return output
