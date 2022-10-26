#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 17:16:07 2020

@author: marcnol

cleans directory structure

In the command line, run as
$ cleanHiM_run.py

to erase all the directories produced recursively

and

$ cleanHiM_run.py --all

to erase also the output MD, Log, and Session files

"""
import argparse
import glob
import os
import shutil

from fileProcessing.fileManagement import Folders, Parameters

# =============================================================================
# MAIN
# =============================================================================

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images, default: .")
    parser.add_argument(
        "-P", "--fileParameters", help="parameters file, default: infoList.json"
    )
    parser.add_argument(
        "-A", "--all", help="Deletes folders, MD files, LOG files", action="store_true"
    )

    args = parser.parse_args()

    if args.rootFolder:
        ROOT_FOLDER = args.rootFolder
    else:
        ROOT_FOLDER = "."

    if args.fileParameters:
        FILE_PARAMETERS = args.fileParameters
    else:
        FILE_PARAMETERS = "infoList.json"

    # removes files in rootFolder
    if args.all:
        markdown_files = glob.glob(
            ROOT_FOLDER + os.sep + "HiM_analysis*.md", recursive=True
        )
        md_log_files = glob.glob(
            ROOT_FOLDER + os.sep + "HiM_analysis*.log", recursive=True
        )
        log_files = glob.glob(ROOT_FOLDER + os.sep + "log*.txt", recursive=True)
        session_files = glob.glob(
            ROOT_FOLDER + os.sep + "Session*.json", recursive=True
        )

        for f in markdown_files + log_files + session_files + md_log_files:
            try:
                os.remove(f)
                print(f"File deleted: {f} ")
            except OSError as e:
                print(f"Error: {f} : {e.strerror}")

    # Removes directories produced during previous runs
    current_param = Parameters(
        root_folder=ROOT_FOLDER, label="", file_name=FILE_PARAMETERS
    )

    data_folder = Folders(current_param.param_dict["rootFolder"])

    for current_folder in data_folder.list_folders:

        folders_to_remove = []
        folders_to_remove.append(
            current_folder + os.sep + current_param.param_dict["zProject"]["folder"]
        )
        folders_to_remove.append(
            current_folder + os.sep + current_param.param_dict["alignImages"]["folder"]
        )
        folders_to_remove.append(
            current_folder
            + os.sep
            + current_param.param_dict["segmentedObjects"]["folder"]
        )
        folders_to_remove.append(current_folder + os.sep + "buildsPWDmatrix")

        for new_folder in folders_to_remove:
            if os.path.isdir(new_folder):
                shutil.rmtree(new_folder)
                print(f"{new_folder} removed")
            else:
                print(f"{new_folder} does not exist")

if __name__ == "__main__":
    main()