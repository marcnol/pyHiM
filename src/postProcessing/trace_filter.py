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

import argparse
import csv
import glob
import json
import os
import select
import sys
from datetime import datetime

import numpy as np

from imageProcessing.imageProcessing import Image
from matrixOperations.chromatin_trace_table import ChromatinTraceTable

# =============================================================================
# FUNCTIONS
# =============================================================================q


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("--N_barcodes", help="minimum_number_barcodes. Default = 2")
    parser.add_argument(
        "--pipe", help="inputs Trace file list from stdin (pipe)", action="store_true"
    )

    p = {}

    args = parser.parse_args()
    if args.rootFolder:
        p["rootFolder"] = args.rootFolder
    else:
        p["rootFolder"] = "."

    if args.N_barcodes:
        p["N_barcodes"] = int(args.N_barcodes)
    else:
        p["N_barcodes"] = 2

    p["trace_files"] = []
    if args.pipe:
        p["pipe"] = True
        if select.select([sys.stdin,], [], [], 0.0)[0]:
            p["trace_files"] = [line.rstrip("\n") for line in sys.stdin]
        else:
            print("Nothing in stdin")
    else:
        p["pipe"] = False

    return p


def runtime(folder, N_barcodes=2, trace_files=[]):

    # gets trace files from buildsPWDmatrix folder
    if len(trace_files) < 1:
        trace_folder = folder.rstrip("/") + os.sep + "buildsPWDmatrix" + os.sep
        trace_files = [
            x
            for x in glob.glob(trace_folder + "Trace*ecsv")
            if "uniqueBarcodes" not in x
        ]

    print(
        "\n{} trace files to process= {}".format(
            len(trace_files), "\n".join(map(str, trace_files))
        )
    )

    if len(trace_files) > 0:

        # iterates over traces in folder
        for trace_file in trace_files:

            trace = ChromatinTraceTable()
            trace.initialize()

            # reads new trace
            trace.load(trace_file)

            # filters trace
            trace.filter_traces_by_n(minimum_number_barcodes=N_barcodes)

            # saves output trace
            outputfile = trace_file.rstrip(".ecsv") + "_filtered" + ".ecsv"
            trace.save(outputfile, trace.data, comments="filtered")
    else:
        print("No trace file found to process!")

    return len(trace_files)


# =============================================================================
# MAIN
# =============================================================================

def main():
    begin_time = datetime.now()

    # [parsing arguments]
    p = parse_arguments()
    # [loops over lists of datafolders]
    folder = p["rootFolder"]
    n_traces_processed = runtime(
        folder, N_barcodes=p["N_barcodes"], trace_files=p["trace_files"]
    )

    print(f"Processed <{n_traces_processed}> trace(s)")
    print("Finished execution")

if __name__ == "__main__":
    main()