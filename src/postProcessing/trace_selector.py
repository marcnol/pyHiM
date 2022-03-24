#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 10:47:29 2022

@author: marcnol

This script will load a trace file and a number of numpy masks and assign them labels

$ trace_selector.py

outputs

ChromatinTraceTable() object and output .ecsv trace table file .


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
import select

from matrixOperations.chromatin_trace_table import ChromatinTraceTable
from imageProcessing.imageProcessing import Image

# =============================================================================
# FUNCTIONS
# =============================================================================q


def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("--pixel_size", help="Lateral pixel size un microns. Default = 0.1")
    parser.add_argument("--pipe", help="inputs Trace file list from stdin (pipe)", action = 'store_true')

    p = {}

    args = parser.parse_args()
    if args.rootFolder:
        p["rootFolder"] = args.rootFolder
    else:
        p["rootFolder"] = "."

    if args.pixel_size:
        p["pixel_size"] = args.pixel_size
    else:
        p["pixel_size"] = 0.1

    p["trace_files"] = list()
    if args.pipe:
        p["pipe"] = True
        if select.select([sys.stdin, ], [], [], 0.0)[0]:
            p["trace_files"] = [line.rstrip("\n") for line in sys.stdin]
        else:
            print("Nothing in stdin")
    else:
        p["pipe"] = False

    return p


def assign_masks(trace, folder_masks, pixel_size=0.1):

    # [checks if DAPI mask exists for the file to process]
    mask_files = glob.glob(folder_masks.rstrip("/") + os.sep + "*.npy")
    mask_files = [x for x in mask_files if "SNDmask" in x.split("_")]

    if len(mask_files) < 1:
        print(f"No mask file found in folder: {folder_masks}")
        return

    for mask_file in mask_files:
        label = mask_file.split("_")[-1].split(".")[0]

        # load mask
        print(f"\nWill attemp to match mask {label} from: {mask_file}")
        mask = Image()
        mask.data_2D = np.load(mask_file, allow_pickle=False).squeeze()

        # matches traces and masks
        index = 0
        labeled_trace = list()
        for trace_row in trace.data:
            x_int = int(trace_row["x"] / pixel_size)
            y_int = int(trace_row["y"] / pixel_size)
            if "x" in trace_row["label"]:
                trace_row["label"] = "_"

            # labels are appended as comma separated lists. Thus a localization can have multiple labels
            if mask.data_2D[x_int, y_int] == 1:
                trace_row["label"] = trace_row["label"] + "," + label
                index += 1
                labeled_trace.append(trace_row["Trace_ID"])
                # print("label assigned: {}, {}".format(trace_row['label'],label))

        unique_traces_labeled = set(labeled_trace)
        print(
            f"\n> {index} trace rows out of {len(trace.data)} were associated to mask {label}. Unique traces: {len(unique_traces_labeled)}"
        )

    return trace


def process_traces(folder, pixel_size=0.1, trace_files = list()):

    trace_folder = folder.rstrip("/") + os.sep + "buildsPWDmatrix" + os.sep
    masks_folder = folder.rstrip("/") + os.sep + "segmentedObjects" + os.sep

    if len(trace_files)<1:
        trace_files = [x for x in glob.glob(trace_folder + "Trace*ecsv") if "uniqueBarcodes" not in x]

    # removes already labeled trace files
    trace_files = [x for x in trace_files if "labeled" not in x]

    print("\n{} trace files to process= {}".format(len(trace_files), "\n".join(map(str, trace_files))))

    if len(trace_files) > 0:
        # iterates over traces in folder
        for trace_file in trace_files:

            trace = ChromatinTraceTable()
            trace.initialize()

            # reads new trace
            trace.load(trace_file)

            trace = assign_masks(trace, masks_folder, pixel_size=pixel_size)

            outputfile = trace_file.rstrip(".ecsv") + "_labeled" + ".ecsv"

            trace.save(outputfile, trace.data, comments="labeled")

    return trace


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    # [parsing arguments]
    p = parseArguments()
    # [loops over lists of datafolders]
    folder = p["rootFolder"]
    traces = process_traces(folder, pixel_size=p["pixel_size"], trace_files = p["trace_files"])

    print("Finished execution")
