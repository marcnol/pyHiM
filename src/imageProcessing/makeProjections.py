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
from core.parameters import Parameters
from core.pyhim_logging import print_log, print_session_name, write_string_to_file
from core.saving import image_show_with_values
from imageProcessing.imageProcessing import Image

warnings.filterwarnings("ignore")

np.seterr(divide="ignore", invalid="ignore")


class Feature:
    def __init__(self, params: Parameters):
        self.m_params = params
        self.required_data = []
        self.required_ref = []
        self.required_table = []
        self.out_folder = ""

    def get_required_inputs(self):
        return self.required_data, self.required_ref, self.required_table


class Project(Feature):
    def __init__(self, params: Parameters):
        super().__init__(params)
        self.required_data = ["barcode", "mask", "dapi", "fiducial", "rna"]
        self.out_folder = "zProject"
        self.out_tag = "_2d"

        self.block_size = params.get_labeled_dict_value("zProject", "blockSize")
        self.display = params.get_labeled_dict_value("zProject", "display")
        self.folder_name = params.get_labeled_dict_value("zProject", "folder")
        self.mode = params.get_labeled_dict_value("zProject", "mode")
        self.operation = params.get_labeled_dict_value("zProject", "operation")
        self.save_image = params.get_labeled_dict_value("zProject", "saveImage")
        self.window_security = params.get_labeled_dict_value(
            "zProject", "windowSecurity"
        )
        self.z_project_option = params.get_labeled_dict_value(
            "zProject", "zProjectOption"
        )
        self.zmax = params.get_labeled_dict_value("zProject", "zmax")
        self.zmin = params.get_labeled_dict_value("zProject", "zmin")
        self.zwindows = params.get_labeled_dict_value("zProject", "zwindows")

    def run(self, img, label: str):
        mode = self.mode[label]
        if mode == "laplacian":
            return self._projection_laplacian(img, label)
        # find the correct range for the projection
        img_reduce = self.precise_z_planes(img, mode, label)
        img_projected = self.projection_2d(img_reduce, label)
        return img_projected

    def check_zmax(self, img_size, label):
        if self.zmax[label] > img_size[0]:
            print_log("$ Setting z max to the last plane")
            self.zmax[label] = img_size[0]

    def precise_z_planes(self, img, mode, label):
        img_size = img.shape
        self.check_zmax(img_size, label)
        if mode == "automatic":
            focus_plane, z_range = self._precise_z_planes_auto(img, label)
        elif mode == "full":
            focus_plane, z_range = self._precise_z_planes_full(img_size)
        elif mode == "manual":
            focus_plane, z_range = self._precise_z_planes_manual(label)
        else:
            raise ValueError(
                f"Projection mode UNRECOGNIZED: {mode}\n> Available mode: automatic,full,manual,laplacian"
            )
        print_log(f"$ Image Size={img_size}")
        print_log(f"$ Focal plane={focus_plane}")
        print_log(f"> Processing z_range:{z_range}")
        return img[z_range[0] : (z_range[-1] + 1)]

    def _precise_z_planes_auto(self, img, label):
        """
        Calculates the focal planes based max standard deviation
        Finds best focal plane by determining the max of the std deviation vs z curve
        """
        win_sec = self.window_security[label]

        print_log("> Calculating planes...")

        nb_of_planes = self.zmax[label] - self.zmin[label]
        std_matrix = np.zeros(nb_of_planes)
        mean_matrix = np.zeros(nb_of_planes)

        # calculate STD in each plane
        for i in range(0, nb_of_planes):
            std_matrix[i] = np.std(img[i])
            mean_matrix[i] = np.mean(img[i])

        max_std = np.max(std_matrix)
        i_focus_plane = np.where(std_matrix == max_std)[0][0]
        # Select a window to avoid being on the edges of the stack

        if i_focus_plane < win_sec or (i_focus_plane > nb_of_planes - win_sec):
            focus_plane = i_focus_plane
        else:
            # interpolate zfocus
            axis_z = range(
                max(
                    self.zmin[label],
                    i_focus_plane - win_sec,
                    min(self.zmax[label], i_focus_plane + win_sec),
                )
            )

            std_matrix -= np.min(std_matrix)
            std_matrix /= np.max(std_matrix)

            try:
                fitgauss = spo.curve_fit(
                    gaussian, axis_z, std_matrix[axis_z[0] : axis_z[-1] + 1]
                )
                focus_plane = int(fitgauss[0][1])
            except RuntimeError:
                print_log("Warning, too many iterations")
                focus_plane = i_focus_plane

        zmin = max(win_sec, focus_plane - win_sec)
        zmax = min(
            nb_of_planes, win_sec + nb_of_planes, focus_plane + self.zwindows[label]
        )
        zrange = range(zmin, zmax + 1)

        return focus_plane, zrange

    def _precise_z_planes_full(self, img_size):
        (zmin, zmax) = (0, img_size[0])
        focus_plane = round((zmin + zmax) / 2)
        z_range = range(zmin, zmax)
        return focus_plane, z_range

    def _precise_z_planes_manual(self, label):
        # Manual: reads from parameters file
        if self.zmin[label] >= self.zmax[label]:
            raise SystemExit(
                "zmin is equal or larger than zmax in configuration file. Cannot proceed."
            )
        focus_plane = round((self.zmin[label] + self.zmax[label]) / 2)
        z_range = range(self.zmin[label], self.zmax[label])
        return focus_plane, z_range

    def _projection_laplacian(self, img, label):
        print_log("Stacking using Laplacian variance...")
        focal_plane_matrix, z_range, block = projection.reinterpolate_focal_plane(
            img, block_size_xy=self.block_size[label], window=self.zwindows[label]
        )
        # reassembles image
        output = projection.reassemble_images(
            focal_plane_matrix, block, window=self.zwindows[label]
        )

        return output, focal_plane_matrix, z_range

    def projection_2d(self, img, label):
        # sums images
        i_collapsed = None
        option = self.z_project_option[label]
        if "MIP" == option:
            # Max projection of selected planes
            i_collapsed = projection.maximum_projection(img)
            # i_collapsed = projection.maximum_projection(img[:-1])
        elif "sum" == option:
            # Sums selected planes
            i_collapsed = projection.sum_projection(img)
        else:
            print_log(
                f"ERROR: mode not recognized. Expected: MIP or sum. Read: {option}"
            )

        return i_collapsed

    # def z_projection_range(self):
    # find the correct range for the projection
    # if self.current_param.param_dict["zProject"]["zmax"] > self.image_size[0]:
    #     print_log("$ Setting z max to the last plane")
    #     self.current_param.param_dict["zProject"]["zmax"] = self.image_size[0]

    # if self.current_param.param_dict["zProject"]["mode"] == "automatic":
    #     print_log("> Calculating planes...")
    #     z_range = calculate_zrange(self.data, self.current_param)

    # elif self.current_param.param_dict["zProject"]["mode"] == "full":
    #     (zmin, zmax) = (0, self.image_size[0])
    #     z_range = (round((zmin + zmax) / 2), range(zmin, zmax))

    # if self.current_param.param_dict["zProject"]["mode"] == "laplacian":
    # print_log("Stacking using Laplacian variance...")
    # (
    #     self.data_2d,
    #     self.focal_plane_matrix,
    #     self.z_range,
    # ) = reinterpolate_focal_plane(self.data, self.current_param.param_dict)
    # self.focus_plane = self.z_range[0]

    # else:
    # # Manual: reads from parameters file
    # (zmin, zmax) = (
    #     self.current_param.param_dict["zProject"]["zmin"],
    #     self.current_param.param_dict["zProject"]["zmax"],
    # )
    # if zmin >= zmax:
    #     raise SystemExit(
    #         "zmin is equal or larger than zmax in configuration file. Cannot proceed."
    #     )
    # z_range = (round((zmin + zmax) / 2), range(zmin, zmax))

    # if "laplacian" not in self.current_param.param_dict["zProject"]["mode"]:
    #     self.data_2d = project_image_2d(
    #         self.data,
    #         z_range,
    #         self.current_param.param_dict["zProject"]["zProjectOption"],
    #     )
    #     self.focus_plane = z_range[0]
    #     self.z_range = z_range[1]

    # print_log(f"> Processing z_range:{self.z_range}")


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
