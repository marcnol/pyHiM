#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 13:40:38 2022

@author: marcnol

This script will 2D project the 3D labeled numpy arrays produced by segmentMasks3D and replace those produced by segmentMasks
s
"""

import os
import numpy as np
import glob
import argparse
from datetime import datetime
from tifffile import imsave


def saves_projections(files,data2D):
    output_files = [x.split('_3Dmasks.npy')[0].rstrip('.')+"_Masks.npy" for x in files]
    output_files_TIFF = [x.split('_3Dmasks.npy')[0].rstrip('.')+"_Masks.tif" for x in files]

    print(f"output files: {output_files}\n\n")

    for output_file,_data2D, output_file_TIFF in zip(output_files,data2D,output_files_TIFF):
        if os.path.exists(output_file):
            print(f"----Warning!----\nRenaming {output_file} as it exists already!\n")
            os.rename(output_file, output_file+"._2Dmasks.npy")

        print(f"> Saving:\n--> {output_file}\n--> {output_file_TIFF} \n")
        np.save(output_file,_data2D)
        imsave(output_file_TIFF,_data2D)

def projects_3D_volumes(data,files):
    numberObjects = [np.amax(x) for x in data]

    data2D = [np.max(x,axis=0) for x in data]

    for numberObject,file in zip(numberObjects,files):
        print(f"\nFile: {os.path.basename(file)}")
        print(f"Number of objects detected: {numberObject}")

    return data2D

def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    args = parser.parse_args()

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = "."

    folder_segmentedObjects = "segmentedObjects"

    folder = rootFolder + os.sep + folder_segmentedObjects

    print(f"Folder: {folder}")
    files = glob.glob(folder+os.sep+"*_3Dmasks.npy")

    return files

def _main():

    begin_time = datetime.now()

    # parses input files
    files = parse_arguments()
    print(f"\n{len(files)} files found matching criteria:\n\n{files}")

    # loads data
    data = [np.load(x) for x in files]

    # projects labeled image
    data2D = projects_3D_volumes(data,files)

    # saves projections
    saves_projections(files,data2D)

    print("Elapsed time: {}".format(datetime.now() - begin_time))



_main()