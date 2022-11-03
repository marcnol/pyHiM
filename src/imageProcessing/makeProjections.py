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

from dask.distributed import get_client, wait

from fileProcessing.fileManagement import Folders, print_log, write_string_to_file
from imageProcessing.imageProcessing import Image

# =============================================================================
# FUNCTIONS
# =============================================================================


def make_2d_projections_file(file_name, current_param, current_session, data_folder):

    if file_name in current_session.data:
        # creates image object
        im_obj = Image(current_param)
        im_obj.load_image_2d(file_name, data_folder.output_folders["zProject"])
        if current_param.param_dict["zProject"]["display"]:
            im_obj.show_image()
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
                "{}\n ![]({})\n".format(os.path.basename(file_name), png_file_name),
                "a",
            )  # initialises MD file

            if current_param.param_dict["zProject"]["mode"] == "laplacian":
                output_name = im_obj.get_image_filename(
                    data_folder.output_folders["zProject"], "_focalPlaneMatrix.png"
                )
                im_obj.image_show_with_values(output_name)

                write_string_to_file(
                    current_param.param_dict["fileNameMD"],
                    "{}\n ![]({})\n".format(os.path.basename(file_name), output_name),
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
    print_log("> Folders read: {}".format(len(data_folder.list_folders)))
    write_string_to_file(
        current_param.param_dict["fileNameMD"],
        "## {}: {}\n".format(
            session_name, current_param.param_dict["acquisition"]["label"]
        ),
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
