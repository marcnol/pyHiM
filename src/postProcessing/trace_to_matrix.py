#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 08:42:12 2023

@author: marcnol

uses the core routines of pyHiM to convert a trace file to a matrix in a standalone script

"""


# =============================================================================
# IMPORTS
# =============================================================================q

import argparse
import select
import sys
from datetime import datetime
import numpy as np

from matrixOperations.build_matrix import BuildMatrix

# =============================================================================
# FUNCTIONS
# =============================================================================q


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--outputFolder", help="Output folder, Default: PWD")
    parser.add_argument("--input", help="Name of input trace file.")
    parser.add_argument("--distance_threshold", help="Threshold for the maximum distance allowed. Default: np.inf")
   
    parser.add_argument(
        "--pipe", help="inputs Trace file list from stdin (pipe)", action="store_true"
    )

    p = {}

    args = parser.parse_args()
    if args.outputFolder:
        p["rootFolder"] = args.outputFolder
    else:
        p["rootFolder"] = "./"

    if args.input:
        p["input"] = args.input
    else:
        p["input"] = None
    
    if args.distance_threshold:
        p["distance_threshold"] = float(args.distance_threshold)
    else:
        p["distance_threshold"] = np.inf

    p["trace_files"] = []
    if args.pipe:
        p["pipe"] = True
        if select.select(
            [
                sys.stdin,
            ],
            [],
            [],
            0.0,
        )[0]:
            p["trace_files"] = [line.rstrip("\n") for line in sys.stdin]
        else:
            print("Nothing in stdin")
    else:
        p["pipe"] = False
        p["trace_files"] = [p["input"]]

    p["colormaps"] = {
        "Nmatrix": "Blues",
        "PWD_KDE": "terrain",
        "PWD_median": "terrain",
        "contact": "coolwarm",
    }

    return p


def runtime(folder, N_barcodes=2, trace_files=[], colormaps=dict(), distance_threshold = np.inf):
    if len(trace_files) < 1:
        print(
            "! Error: no trace file provided. Please either use pipe or the --input option to provide a filename."
        )
        return 0
    elif len(trace_files) == 1:
        print("\n$ trace file to process= {}".format(trace_files))
    else:
        print(
            "\n{} trace files to process= {}".format(
                len(trace_files), "\n".join(map(str, trace_files))
            )
        )

    if len(trace_files) > 0:
        # iterates over traces in folder
        for trace_file in trace_files:
            # converts trace to matrix

            param = dict()
            new_matrix = BuildMatrix(param, colormaps=colormaps)
            new_matrix.launch_analysis(trace_file,distance_threshold=distance_threshold)

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
        folder, trace_files=p["trace_files"], colormaps=p["colormaps"],distance_threshold=p["distance_threshold"]
    )

    print(f"Processed <{n_traces_processed}> trace(s)")
    print("Finished execution")


if __name__ == "__main__":
    main()
