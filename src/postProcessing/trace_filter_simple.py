#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 09:57:36 2022

@author: marcnol

This script will load a trace file and a filter traces based on a series of user arguments

--> Usage

$ trace_filter.py --input Trace.ecsv --z_min 4 --z_max 5 --y_max 175 --output 'zy_filtered' --N_barcodes 3

will analyze 'Trace.ecsv' and remove spots with 4>z>5 amd z>175 and less than 3 barcodes

--> outputs

.ecsv trace table file.

"""

# =============================================================================
# IMPORTS
# =============================================================================q

import argparse
import os
import select
import sys
from datetime import datetime

import numpy as np

from matrixOperations.chromatin_trace_table import ChromatinTraceTable

# =============================================================================
# FUNCTIONS
# =============================================================================q


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-O", "--output", help="Tag to add to the output file. Default = filtered")
    parser.add_argument(
        "--pipe", help="inputs Trace file list from stdin (pipe)", action="store_true"
    )
    parser.add_argument("--input", help="Name of input trace file.")
    parser.add_argument("--N_barcodes", help="minimum_number_barcodes. Default = 2")
    parser.add_argument("--z_min", help="Z minimum for a localization. Default = 0")
    parser.add_argument("--z_max", help="Z maximum for a localization. Default = np.inf")
    parser.add_argument("--y_min", help="Y minimum for a localization. Default = 0")
    parser.add_argument("--y_max", help="Y maximum for a localization. Default = np.inf")
    parser.add_argument("--x_min", help="X minimum for a localization. Default = 0")
    parser.add_argument("--x_max", help="X maximum for a localization. Default = np.inf")

    p = {}

    args = parser.parse_args()
    if args.output:
        p["output"] = args.output
    else:
        p["output"] = "filtered"

    if args.input:
        p["input"] = args.input
    else:
        p["input"] = None

    if args.N_barcodes:
        p["N_barcodes"] = int(args.N_barcodes)
    else:
        p["N_barcodes"] = 2

    if args.z_min:
        p["z_min"] = float(args.z_min)
    else:
        p["z_min"] = 0

    if args.z_max:
        p["z_max"] = float(args.z_max)
    else:
        p["z_max"] = np.inf

    if args.y_min:
        p["y_min"] = float(args.y_min)
    else:
        p["y_min"] = 0

    if args.y_max:
        p["y_max"] = float(args.y_max)
    else:
        p["y_max"] = np.inf

    if args.x_min:
        p["x_min"] = float(args.x_min)
    else:
        p["x_min"] = 0
        
    if args.x_max:
        p["x_max"] = float(args.x_max)
    else:
        p["x_max"] = np.inf
        
        
    p["trace_files"] = []
    if args.pipe:
        p["pipe"] = True
        if select.select([sys.stdin,], [], [], 0.0)[0]:
            p["trace_files"] = [line.rstrip("\n") for line in sys.stdin]
        else:
            print("Nothing in stdin")
    else:
        p["pipe"] = False
        p["trace_files"] = [p["input"]]
        
    return p


def runtime(trace_files=[], N_barcodes=2, coor_limits = dict(), tag = "filtered"):

    # checks number of trace files
    if len(trace_files) < 1:
        print("! Error: no trace file provided. Please either use pipe or the --input option to provide a filename.")
        return 0
    elif len(trace_files) == 1:
        print("\n$ trace files to process= {}".format(trace_files))
    else:
        print(
            "\n{} trace files to process= {}".format(
                len(trace_files), "\n".join(map(str, trace_files))
            )
        )

    coors = ['x','y','z']
    if len(trace_files) > 0:

        # iterates over traces 
        for trace_file in trace_files:

            trace = ChromatinTraceTable()
            trace.initialize()

            # reads new trace
            trace.load(trace_file)

            # filters trace by minimum number of barcodes
            trace.filter_traces_by_n(minimum_number_barcodes=N_barcodes)

            # filters trace by coordinate
            for coor in coors:
                trace.filter_traces_by_coordinate(coor = coor, coor_min = coor_limits[coor+'_min'], coor_max = coor_limits[coor+'_max'])

            # saves output trace
            outputfile = trace_file.rstrip(".ecsv") + "_" + tag + ".ecsv"
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
    n_traces_processed = runtime(
        trace_files=p["trace_files"], N_barcodes=p["N_barcodes"], coor_limits = p, tag=p["output"]
    )

    print(f"Processed <{n_traces_processed}> trace(s)")
    print("Finished execution")


if __name__ == "__main__":
    main()
