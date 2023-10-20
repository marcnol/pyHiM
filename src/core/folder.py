#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Manage folders for outputs"""

import os

from core.data_manager import create_folder

# =============================================================================
# CLASSES
# =============================================================================


class Folders:
    """Used to create and access to outputs folders."""

    def __init__(self, master_folder=r"/home/marcnol/Documents/Images"):
        self.master_folder = master_folder

        # list of sub-folders in rootFilder with images
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

    def create_folders(self, data_path, current_param):
        """
        Creates folders for outputs.
        this function will create all the folders required for processingPipeline

        Parameters
        ----------
        data_path : string
            root_folder
        current_param : Parameters Class
            with filenames of folders to be created

        Returns
        -------
        None.

        """

        self._create_folder(data_path, current_param, "alignImages")

        if "segmentedObjects" in current_param.param_dict.keys():
            self._create_folder(data_path, current_param, "segmentedObjects")

        # backwards compatibility
        if "buildsPWDmatrix" in current_param.param_dict.keys():
            self.output_folders["buildsPWDmatrix"] = (
                data_path
                + os.sep
                + current_param.param_dict["buildsPWDmatrix"]["folder"]
            )
        else:
            self.output_folders["buildsPWDmatrix"] = (
                data_path + os.sep + "buildsPWDmatrix"
            )
        self.create_folder_with_key("buildsPWDmatrix")
        self.output_files["buildsPWDmatrix"] = (
            self.output_folders["buildsPWDmatrix"] + os.sep + "buildsPWDmatrix"
        )

    # TODO: apply this method for each folder creation inside `create_folders`
    def _create_folder(self, data_path, current_param, arg2):
        self.output_folders[arg2] = (
            data_path + os.sep + current_param.param_dict[arg2]["folder"]
        )
        self.create_folder_with_key(arg2)
        self.output_files[arg2] = (
            self.output_folders[arg2] + os.sep
        ) + current_param.param_dict[arg2]["outputFile"]
