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
import scipy.optimize as spo
from apifish.stack import projection
from dask.distributed import get_client, wait

from core.folder import Folders
from core.pyhim_logging import print_log, write_string_to_file
from core.saving import image_show_with_values
from imageProcessing.imageProcessing import Image

warnings.filterwarnings("ignore")

np.seterr(divide="ignore", invalid="ignore")

# =============================================================================
# FUNCTIONS
# =============================================================================


def make_2d_projections_file(file_name, current_param, current_session, data_folder):
    if file_name in current_session.data:
        print_log("# File already projected: {}".format(os.path.basename(file_name)))
    else:
        print_log("\n> Analysing file: {}".format(os.path.basename(file_name)))

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
    print_log("\n===================={}====================\n".format(session_name))
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
        print_log("> Processing Folder: {}".format(current_folder))
        print_log(
            "> About to process {} files\n".format(len(current_param.files_to_process))
        )

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

                print_log("$ Waiting for {} threads to complete ".format(len(threads)))
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


def project_image_2d(img, z_range, mode):
    # sums images
    i_collapsed = None

    if "MIP" in mode:
        # Max projection of selected planes
        i_collapsed = projection.maximum_projection(img[z_range[1][0] : z_range[1][-1]])
    elif "sum" in mode:
        # Sums selected planes
        i_collapsed = projection.sum_projection(
            img[z_range[1][0] : (z_range[1][-1] + 1)]
        )
    else:
        print_log(f"ERROR: mode not recognized. Expected: MIP or sum. Read: {mode}")

    return i_collapsed


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


# @jit(nopython=True)
def calculate_zrange(idata, parameters):
    """
    Calculates the focal planes based max standard deviation
    Finds best focal plane by determining the max of the std deviation vs z curve
    """
    num_planes = (
        parameters.param_dict["zProject"]["zmax"]
        - parameters.param_dict["zProject"]["zmin"]
    )
    std_matrix = np.zeros(num_planes)
    mean_matrix = np.zeros(num_planes)

    # calculate STD in each plane
    for i in range(0, num_planes):
        std_matrix[i] = np.std(idata[i])
        mean_matrix[i] = np.mean(idata[i])

    max_std = np.max(std_matrix)
    i_focus_plane = np.where(std_matrix == max_std)[0][0]
    # Select a window to avoid being on the edges of the stack

    if i_focus_plane < parameters.param_dict["zProject"]["windowSecurity"] or (
        i_focus_plane > num_planes - parameters.param_dict["zProject"]["windowSecurity"]
    ):
        focus_plane = i_focus_plane
    else:
        # interpolate zfocus
        axis_z = range(
            max(
                parameters.param_dict["zProject"]["zmin"],
                i_focus_plane - parameters.param_dict["zProject"]["windowSecurity"],
                min(
                    parameters.param_dict["zProject"]["zmax"],
                    i_focus_plane + parameters.param_dict["zProject"]["windowSecurity"],
                ),
            )
        )

        std_matrix -= np.min(std_matrix)
        std_matrix /= np.max(std_matrix)

        try:
            fitgauss = spo.curve_fit(
                gaussian, axis_z, std_matrix[axis_z[0] : axis_z[-1] + 1]
            )
            # print_log("Estimation of focal plane (px): ", int(fitgauss[0][1]))
            focus_plane = int(fitgauss[0][1])
        except RuntimeError:
            print_log("Warning, too many iterations")
            focus_plane = i_focus_plane

    zmin = max(
        parameters.param_dict["zProject"]["windowSecurity"],
        focus_plane - parameters.param_dict["zProject"]["zwindows"],
    )
    zmax = min(
        num_planes,
        parameters.param_dict["zProject"]["windowSecurity"] + num_planes,
        focus_plane + parameters.param_dict["zProject"]["zwindows"],
    )
    zrange = range(zmin, zmax + 1)

    return focus_plane, zrange


# =============================================================================
# FOCAL PLANE INTERPOLATION
# =============================================================================


def reinterpolate_focal_plane(data, param_dict):
    if "blockSize" in param_dict["zProject"]:
        block_size_xy = param_dict["zProject"]["blockSize"]
    else:
        block_size_xy = 128

    if "zwindows" in param_dict["zProject"]:
        window = param_dict["zProject"]["zwindows"]
    else:
        window = 0

    focal_plane_matrix, z_range, block = projection.reinterpolate_focal_plane(
        data, block_size_xy=block_size_xy, window=window
    )

    # reassembles image
    output = projection.reassemble_images(focal_plane_matrix, block, window=window)

    return output, focal_plane_matrix, z_range


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
