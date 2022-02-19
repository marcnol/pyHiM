#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 10:47:29 2022

@author: marcnol

This script will load a trace file and a number of numpy masks and assign them to each individual trace

$ trace_selector.py

outputs

chromatin_trace_table() object and output .ecsv formatted table file with assembled trace tables.


"""

# =============================================================================
# IMPORTS
# =============================================================================q

import numpy as np
import os, sys
import json
from datetime import datetime
import argparse
import csv
import glob

from matrixOperations.chromatin_trace_table import chromatin_trace_table
from imageProcessing.imageProcessing import Image

# =============================================================================
# FUNCTIONS
# =============================================================================q

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument(
        "-P", "--parameters", help="Provide name of parameter files. folders2Load.json assumed as default",
    )
    parser.add_argument("-A", "--label", help="Add name of label (e.g. doc)")
    parser.add_argument("-W", "--action", help="Select: [all], [labeled] or [unlabeled] cells plotted ")
    parser.add_argument("--saveMatrix", help="Use to load matlab formatted data", action="store_true")
    parser.add_argument("--ndims", help="Dimensions of trace")
    parser.add_argument("--method", help="Method or mask ID used for tracing: KDtree, mask, mask0")

    p={}

    args = parser.parse_args()
    if args.rootFolder:
        p["rootFolder"] = args.rootFolder
    else:
        p["rootFolder"] = "."

    if args.parameters:
        p["parametersFileName"] = args.parameters
    else:
        p["parametersFileName"] = "folders2Load.json"

    if args.label:
        p["label"] = args.label
    else:
        p["label"] = "M"

    if args.action:
        p["action"] = args.action
    else:
        p["action"] = "all"

    if args.saveMatrix:
        p["saveMatrix"] = True
    else:
        p["saveMatrix"] = False

    if args.ndims:
        p["ndims"] = args.ndims
    else:
        p["ndims"] = 3

    if args.method:
        p["method"] = args.method
    else:
        p["method"] = "mask"

    return p


def assign_masks(trace, folder_masks, pixel_size = 0.1):

    numberFilesProcessed = 0

    # [checks if DAPI mask exists for the file to process] # NOT SURE IF THIS WILL ALWAYS FIND THE MASK FILES
    mask_files = glob.glob(folder_masks + os.sep + "_Masks.npy")

    if len(mask_files) < 1:
        print(f"No mask file found in folder: {folder_masks}")
        return

    for mask_file in mask_files:
        label = mask_file.split("_")[-1].split(".")[0] # NOT SURE IF THIS WILL GET THE RIGHT LABEL OFF THE FILENAME

        # load mask
        print(f"\nWill attemp to match mask {label} from: {mask_file}")
        mask = Image()
        mask.data_2D = np.load(mask_file, allow_pickle=False).squeeze()

        # matches traces and masks
        for trace_row in trace.data:
            x_int = int(trace_row['x']/pixel_size)
            y_int = int(trace_row['y']/pixel_size)

            if mask.data_2D[x_int,y_int] == 1:
                trace_row['label'] = trace_row['label'] + "," + label

def process_traces(folder):

    trace_folder = folder.rstrip('/') + os.sep + 'buildsPWDmatrix' + os.sep
    masks_folder = folder.rstrip('/') + os.sep + 'segmentedObjects' + os.sep

    trace_files = [x for x in glob.glob(trace_folder + 'Trace*ecsv') if 'uniqueBarcodes' not in x]


    print(f"\nTrace files to process= {trace_files}")

    if len(trace_files) > 0:
        # iterates over traces in folder
        for trace_file in trace_files:

            trace = chromatin_trace_table()
            trace.initialize()

            #reads new trace
            trace.load(trace_file)

            assign_masks(trace, masks_folder, pixel_size = 0.1)

            outputfile = trace_file.rstrip('.ecsv') + "_labeled" + '.ecsv'

            traces.save(outputfile,traces.data, comments = 'labeled')

    return traces

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    # [parsing arguments]
    p=parseArguments()
    positionROIinformation = p["positionROIinformation"]

    # [loops over lists of datafolders]
    folder = p['rootFolder']
    traces = process_traces(folder)

    print("Finished execution")
