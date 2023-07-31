#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Manage folders for outputs"""

import glob
import os
import re

# =============================================================================
# CLASSES
# =============================================================================


class Folders:
    """Used to create and access to outputs folders."""

    def __init__(self, master_folder=r"/home/marcnol/Documents/Images"):
        self.master_folder = master_folder
        self.list_folders = []

        # list of sub-folders in rootFilder with images
        self.z_project_folder = ""
        self.output_folders = {}
        self.output_files = {}

        self.set_folders()

    def set_folders(self):
        """returns list of directories"""
        self.list_folders = [self.master_folder]

    def create_folder_with_key(self, folder_key_name: str):
        """Create one folder for one type of pyHiM outputs.

        Parameters
        ----------
        folder_key_name : str
            The key word to access at the folder name inside the parameter file
        """
        folder_path = self.output_folders[folder_key_name]
        create_single_folder(folder_path)

    def create_folders(self, files_folder, current_param):
        """
        Creates folders for outputs.
        this function will create all the folders required for processingPipeline

        Parameters
        ----------
        files_folder : string
            root_folder
        current_param : Parameters Class
            with filenames of folders to be created

        Returns
        -------
        None.

        """
        self.output_folders["zProject"] = (
            files_folder + os.sep + current_param.param_dict["zProject"]["folder"]
        )
        self.create_folder_with_key("zProject")

        self._create_folder(files_folder, current_param, "alignImages")

        self.output_files["dictShifts"] = (
            self.output_folders["alignImages"]
            + os.sep
            + current_param.param_dict["alignImages"]["outputFile"]
        )

        if "segmentedObjects" in current_param.param_dict.keys():
            self._create_folder(files_folder, current_param, "segmentedObjects")

        # backwards compatibility
        if "buildsPWDmatrix" in current_param.param_dict.keys():
            self.output_folders["buildsPWDmatrix"] = (
                files_folder
                + os.sep
                + current_param.param_dict["buildsPWDmatrix"]["folder"]
            )
        else:
            self.output_folders["buildsPWDmatrix"] = (
                files_folder + os.sep + "buildsPWDmatrix"
            )
        self.create_folder_with_key("buildsPWDmatrix")
        self.output_files["buildsPWDmatrix"] = (
            self.output_folders["buildsPWDmatrix"] + os.sep + "buildsPWDmatrix"
        )

    # TODO: apply this method for each folder creation inside `create_folders`
    def _create_folder(self, files_folder, current_param, arg2):
        self.output_folders[arg2] = (
            files_folder + os.sep + current_param.param_dict[arg2]["folder"]
        )
        self.create_folder_with_key(arg2)
        self.output_files[arg2] = (
            self.output_folders[arg2] + os.sep
        ) + current_param.param_dict[arg2]["outputFile"]


def create_single_folder(folder_path):
    """Create folder with `makedirs` from os module.
    It's a recursive directory creation function.

    Parameters
    ----------
    folder_path : str
        Relative or absolute path + name of folder
    """
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"$ Folder '{folder_path}' created successfully.")
    else:
        print(f"! Folder '{folder_path}' already exists.")


def retrieve_number_rois_folder(root_folder, reg_exp, ext="tif"):
    """
    given a directory and a regular expression, it returns the number of unique rois detected

    Parameters
    ----------
    root_folder : string
    ext : string, optional
        File extension. The default is 'tif'.

    Returns
    -------
    list
        list of unique rois.

    """
    files = glob.glob(root_folder + os.sep + "*" + ext)
    rois = [re.search(reg_exp, x)["roi"] for x in files]
    return unique(rois)


def unique(list1):
    """function to get unique values"""
    # intilize a null list
    unique_list = []

    # traverse for all elements
    for val in list1:
        # check if exists in unique_list or not
        if val not in unique_list:
            unique_list.append(val)

    return unique_list
