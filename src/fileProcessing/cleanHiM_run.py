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
import os
import glob
import argparse
from fileManagement import Parameters, folders
import shutil

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images, default: .")
    parser.add_argument("-P", "--fileParameters", help="parameters file, default: infoList.json")
    parser.add_argument("-A", "--all", help="Deletes folders, MD files, LOG files", action="store_true")

    args = parser.parse_args()

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = "."
        # rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset'

    if args.fileParameters:
        fileParameters = args.fileParameters
    else:
        fileParameters = "infoList.json"

    # removes files in rootFolder
    if args.all:
        filesMD = glob.glob(rootFolder + os.sep + "HiM_analysis*.md", recursive=True)
        filesLOGMD = glob.glob(rootFolder + os.sep + "HiM_analysis*.log", recursive=True)
        filesLOG = glob.glob(rootFolder + os.sep + "log*.txt", recursive=True)
        filesSession = glob.glob(rootFolder + os.sep + "Session*.json", recursive=True)

        for f in filesMD + filesLOG + filesSession + filesLOGMD:
            try:
                os.remove(f)
                print("File deleted: {} ".format(f))
            except OSError as e:
                print("Error: {} : {}".format(f, e.strerror))

    # Removes directories produced during previous runs
    param = Parameters(rootFolder=rootFolder, label='', fileName = fileParameters)

    dataFolder = folders(param.param["rootFolder"])

    for currentFolder in dataFolder.listFolders:

        folders2Remove = []
        folders2Remove.append(currentFolder + os.sep + param.param["zProject"]["folder"])
        folders2Remove.append(currentFolder + os.sep + param.param["alignImages"]["folder"])
        folders2Remove.append(currentFolder + os.sep + param.param["segmentedObjects"]["folder"])
        folders2Remove.append(currentFolder + os.sep + "buildsPWDmatrix")
        folders2Remove.append(currentFolder + os.sep + param.param["projectsBarcodes"]["folder"])

        for newFolder in folders2Remove:
            if os.path.isdir(newFolder):
                shutil.rmtree(newFolder)
                print("{} removed".format(newFolder))
            else:
                print("{} does not exist".format(newFolder))
