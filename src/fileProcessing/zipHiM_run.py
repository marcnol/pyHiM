#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 11:14:47 2020

@author: marcnol

zips all png, MD and output files from a HiM run. It excludes .npy and Tiff images.

In the command line, run as 
$ zipHiM_run.py 

to zip all the directories recursively

"""
import os
import glob
import argparse
from fileManagement import Parameters, folders
import shutil
import tarfile

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images, default: .")
    parser.add_argument("-P", "--fileParameters", help="parameters file, default: infoList_barcode.json")

    args = parser.parse_args()

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = os.getcwd()

    if args.fileParameters:
        fileParameters = args.fileParameters
    else:
        fileParameters = "infoList.json"

    # opens tarfile
    os.chdir(rootFolder)
    tarFileName = "HiMrun.tar"
    print("creating archive: {}".format(tarFileName))

    # tar files in rootFolder
    filesMD = [os.path.basename(f) for f in glob.glob(rootFolder + os.sep + "HiM_analysis*.md", recursive=True)]
    filesLOGMD = [os.path.basename(f) for f in glob.glob(rootFolder + os.sep + "log*.txt", recursive=True)]
    filesLOG = [os.path.basename(f) for f in glob.glob(rootFolder + os.sep + "HiM_analysis*.log", recursive=True)]
    filesSession = [os.path.basename(f) for f in glob.glob(rootFolder + os.sep + "Session*.json", recursive=True)]

    tarcmd = "tar -cvf " + tarFileName + " " + " ".join(filesMD + filesLOG + filesSession + filesLOGMD)
    print("Archiving:\n{}".format("\n".join(filesMD + filesLOG + filesSession + filesLOGMD)))

    os.system(tarcmd)

    # tars directories produced during previous runs
    param = Parameters(rootFolder, fileParameters)

    dataFolder = folders(param.param["rootFolder"])

    for currentFolder in dataFolder.listFolders:

        folders2Remove = []
        folders2Remove.append(currentFolder + os.sep + param.param["zProject"]["folder"])
        folders2Remove.append(currentFolder + os.sep + param.param["alignImages"]["folder"])
        folders2Remove.append(currentFolder + os.sep + param.param["segmentedObjects"]["folder"])
        folders2Remove.append(currentFolder + os.sep + "buildsPWDmatrix")
        folders2Remove.append(currentFolder + os.sep + param.param["projectsBarcodes"]["folder"])

        for newFolder in folders2Remove:
            if rootFolder == ".":
                newFolderRelative = "." + newFolder.split(os.getcwd())[1]
            else:
                newFolderRelative = "." + newFolder.split(rootFolder)[1]

            fileExtensions = ["/*.png", "/*.dat", "/*.ecsv", "/buildsPWDmatrix*.npy"]
            for newFileExtension in fileExtensions:

                newFiles = newFolderRelative + newFileExtension

                if len(glob.glob(newFiles)) > 0:
                    tarcmd = "tar -rf " + tarFileName + " " + newFiles
                    os.system(tarcmd)
                    print("Archiving: {}".format(newFiles))

    if os.path.exists(tarFileName):
        print("Zipping {}".format(tarFileName))
        os.system("gzip " + tarFileName)
