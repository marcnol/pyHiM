#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Manage folders for outputs"""

import glob
import os
import re

from core.pyhim_logging import print_log
from core.data_manager import create_folder

# =============================================================================
# CLASSES
# =============================================================================


class Folders:
    """Used to create and access to outputs folders."""

    def __init__(self, master_folder=r"/home/marcnol/Documents/Images"):
        self.master_folder = master_folder

        # list of sub-folders in rootFilder with images
        self.z_project_folder = ""
        self.output_folders = {}
        self.output_files = {}

    def create_folder_with_key(self, folder_key_name: str):
        """Create one folder for one type of pyHiM outputs.

        Parameters
        ----------
        folder_key_name : str
            The key word to access at the folder name inside the parameter file
        """
        folder_path = self.output_folders[folder_key_name]
        create_folder(folder_path)
        create_folder(folder_path + os.sep + "data")

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
            + "data"
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


def unique(list_to_sort: list):
    """function to get unique values"""
    return list(set(list_to_sort))
