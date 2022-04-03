#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 09:57:36 2022

@author: marcnol

This script will load a trace file and a filter traces based on a series of user arguments

$ trace_filter.py

outputs

.ecsv trace table file.

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
    parser.add_argument("--N_barcodes", help="minimum_number_barcodes. Default = 2")
    parser.add_argument("--pipe", help="inputs Trace file list from stdin (pipe)", action = 'store_true')

    p = {}

    args = parser.parse_args()
    if args.rootFolder:
        p["rootFolder"] = args.rootFolder
    else:
        p["rootFolder"] = "."

    if args.N_barcodes:
        p["N_barcodes"] = args.N_barcodes
    else:
        p["N_barcodes"] = 2

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

def filter_traces_by_n(trace_table, minimum_number_barcodes = 2):
    """
    Removes rows in trace table with less than `minimum_number_barcodes` barcodes

    Parameters
    ----------
    trace_table : ASTROPY Table
        input trace table.
    minimum_number_barcodes : TYPE, optional
        minimum number of barcodes in trace. The default is 1.

    Returns
    -------
    trace_table : ASTROPY Table
        output trace table.

    """

    # indexes trace file
    trace_table_indexed = trace_table.group_by("Trace_ID")

    # iterates over traces
    print(f"\n$ WIll keep traces with {minimum_number_barcodes } spots")
    print(f"$ Number of original spots / traces: {len(trace_table)} / {len(trace_table_indexed.groups)}")

    barcodes_to_remove = list()

    for idx, trace in enumerate(trace_table_indexed.groups):

        if len(trace["Trace_ID"].data) < minimum_number_barcodes:
            barcodes_to_remove.append(list(trace["Spot_ID"].data))

    print(f"$ Number of barcodes to remove: {len(barcodes_to_remove)}")

    list_barcode_to_remove=list()
    for barcodes in barcodes_to_remove:
        for x in barcodes:
            list_barcode_to_remove.append(x)

    rows_to_remove = list()

    for idx, row in enumerate(trace_table):
        spot_id = row["Spot_ID"]

        if spot_id in list_barcode_to_remove:
            rows_to_remove.append(idx)


    trace_table.remove_rows(rows_to_remove)
    trace_table_indexed = trace_table.group_by("Trace_ID")
    print(f"$ Number of spots / traces left: {len(trace_table)} / {len(trace_table_indexed.groups)}")

    return trace_table


def runtime(folder, N_barcodes = 2, trace_files = list()):

    # gets trace files from buildsPWDmatrix folder
    if len(trace_files) < 1:
        trace_folder = folder.rstrip("/") + os.sep + "buildsPWDmatrix" + os.sep
        trace_files = [x for x in glob.glob(trace_folder + "Trace*ecsv") if "uniqueBarcodes" not in x]

    print("\n{} trace files to process= {}".format(len(trace_files), "\n".join(map(str, trace_files))))

    if len(trace_files) > 0:

        # iterates over traces in folder
        for trace_file in trace_files:

            trace = ChromatinTraceTable()
            trace.initialize()

            # reads new trace
            trace.load(trace_file)

            # filters trace
            trace = filter_traces_by_n(trace_file, minimum_number_barcodes = N_barcodes)

            # saves output trace
            outputfile = trace_file.rstrip(".ecsv") + "_filtered" + ".ecsv"
            trace.save(outputfile, trace.data, comments="filtered")
    else:
        print("No trace file found to process!")

    return len(trace_files)


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    # [parsing arguments]
    p = parseArguments()
    # [loops over lists of datafolders]
    folder = p["rootFolder"]
    n_traces_processed = runtime(folder, N_barcodes=p["N_barcodes"], trace_files = p["trace_files"])

    print(f"Processed {n_traces_processed} traces")
    print("Finished execution")
