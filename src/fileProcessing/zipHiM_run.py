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
import shutil
import tarfile

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images, default: .")
    parser.add_argument("-P", "--fileParameters", help="parameters file, default: infoList_barcode.json")
    parser.add_argument("-R", "--recursive", help="One more depth of folders will be explored and zipped", action="store_true")

    args = parser.parse_args()

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = os.getcwd()

    if args.fileParameters:
        fileParameters = args.fileParameters
    else:
        fileParameters = "infoList.json"

    if args.recursive:
        recursive= args.recursive
    else:
        recursive= False

    print("RootFolders: {}".format(rootFolder))

    # opens tarfile
    os.chdir(rootFolder)
    tarFileName = "HiMrun.tar"
    print("creating archive: {} in {}".format(tarFileName,rootFolder))

    # tar files in rootFolder
    filesMD = [os.path.basename(f) for f in glob.glob(rootFolder + os.sep + "*.md")]
    filesLOG = [os.path.basename(f) for f in glob.glob(rootFolder + os.sep + "*.log")]
    filesSession = [os.path.basename(f) for f in glob.glob(rootFolder + os.sep + "*.json")]

    tarcmd = "tar -cvf " + tarFileName + " " + " ".join(filesMD + filesLOG + filesSession)
    print("Archiving:\n{}".format("\n".join(filesMD + filesLOG + filesSession )))

    os.system(tarcmd)

    # tars directories produced during previous runs
    if recursive:
        folders = glob.glob(rootFolder + os.sep + "*")
        folders = [x for x in folders if os.path.isdir(x)] # keeps only folders
        folders = [x for x in folders if os.path.exists(x+os.sep+"infoList.json")]

    else:
        folders  = [rootFolder]

    print("Folders to zip:\n{}".format(folders))
    print("="*30)
    for currentFolder in folders:

        folders2zip = []
        folders2zip.append(currentFolder + os.sep + "zProject")
        folders2zip.append(currentFolder + os.sep + "alignImages")
        folders2zip.append(currentFolder + os.sep + "segmentedObjects")
        folders2zip.append(currentFolder + os.sep + "buildsPWDmatrix")
        folders2zip.append(currentFolder + os.sep + "projectsBarcodes")

        print("sub-folders to zip:\n{}".format(folders2zip))

        for newFolder in folders2zip:
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

        print("-"*30)

    if os.path.exists(tarFileName):
        print("Zipping {}".format(tarFileName))
        os.system("gzip " + tarFileName)
