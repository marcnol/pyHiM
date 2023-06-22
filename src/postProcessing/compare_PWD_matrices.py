#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 21:04:10 2023

@author: marcnol

script to compare PWD matrices from two experiments
*- pearson correlation in ensemble *
- make sure to map barcodes to allow for experiments with different barcode combinations [TODO]
- same but single cell [TODO]


"""


# =============================================================================
# IMPORTS
# =============================================================================q

import argparse
import select
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import seaborn as sns

from matrixOperations.HIMmatrixOperations import (
    calculate_contact_probability_matrix,
    calculate_ensemble_pwd_matrix,
)

# =============================================================================
# FUNCTIONS
# =============================================================================q


def usage():
    print("-" * 80)
    print("Example usage:\n")
    print("$ ls buildsPWDmatrix/*HiM*npy")
    print(
        "buildsPWDmatrix/buildsPWDmatrix_DAPI_HiMscMatrix.npy  buildsPWDmatrix/buildsPWDmatrix_mask0_HiMscMatrix.npy"
    )
    print("$ ls buildsPWDmatrix/*HiM*npy | compare_PWD_matrices.py\n")
    print("You should now see a file called scatter_plot.png\n")
    print("-" * 80)


def parse_arguments():
    parser = argparse.ArgumentParser(
        "This script will compare two PWD matrices. Please provide input matrices by --input1 and --input2.\n Example: $ compare_PWD_matrices --input1 matrix1.npy --input2 matrix2.npy"
    )
    parser.add_argument("--input1", help="Name of first input trace file.")
    parser.add_argument("--input2", help="Name of second input trace file.")
    parser.add_argument(
        "--mode",
        help="Mode used to calculate the mean distance. Can be either 'median', 'KDE' or 'proximity'. Default: median",
    )
    parser.add_argument("--x_min", help="X minimum for a localization. Default = 0")
    parser.add_argument(
        "--x_max", help="X maximum for a localization. Default = np.inf"
    )
    parser.add_argument(
        "--max_distance",
        help="Upper distance threshold in micro-meters. Default = np.inf",
    )
    parser.add_argument("--scale", help="Scale: log or linear. Default: linear")

    p = {}

    args = parser.parse_args()

    if args.input1:
        p["input1"] = args.input1
    else:
        p["input1"] = None

    if args.input2:
        p["input2"] = args.input2
    else:
        p["input2"] = None

    if args.mode:
        p["mode"] = args.mode
    else:
        p["mode"] = "median"

    if args.x_min:
        p["x_min"] = float(args.x_min)
    else:
        p["x_min"] = None

    if args.x_max:
        p["x_max"] = float(args.x_max)
    else:
        p["x_max"] = None

    if args.max_distance:
        p["max_distance"] = float(args.max_distance)
    else:
        p["max_distance"] = np.inf

    if args.scale:
        p["scale"] = args.scale
    else:
        p["scale"] = "linear"

    p["input_files"] = list()
    p["input_files"].append(p["input1"])
    p["input_files"].append(p["input2"])

    print("Input parameters\n" + "=" * 16)
    for item in p.keys():
        print("{}-->{}".format(item, p[item]))

    return p


def plot_result(x, y, r, p):
    if p["x_min"] is None:
        x_min = np.min([x, y])
    else:
        x_min = p["x_min"]

    if p["x_max"] is None:
        x_max = np.max([x, y])
    else:
        x_max = p["x_max"]

    print(f"$ limits: {x_min}-->{x_max}")
    plt.figure(figsize=(15, 15))
    plt.rcParams.update({"font.size": 20})

    plt.scatter(
        x,
        y,
        label="Pearson = {:.2f}".format(r),
        color="m",
        marker="o",
        s=30,
        linewidth=5,
    )

    if p["scale"] == "log":
        plt.xscale("log")
        plt.yscale("log")

    plt.title("Correlation of {} values".format(p["mode"]))
    plt.xlabel("dataset 1")
    plt.ylabel("dataset 2")
    plt.legend()

    axisX = np.linspace(x_min, x_max, 100)
    plt.plot(axisX, axisX, color="red", linewidth=3)
    plt.xlim([x_min, x_max])
    plt.ylim([x_min, x_max])

    # plt.show()
    filename = "scatter_plot.png"
    plt.savefig(filename)
    print(f"> Output image saved as : {filename}")


def calculates_pearson_correlation(x, y):
    r, p = scipy.stats.pearsonr(x, y)
    return r


def parses_matrix_to_vector(matrix):
    matrix_shape = matrix.shape
    matrix = np.nan_to_num(matrix)
    return matrix.reshape(matrix_shape[0] * matrix_shape[1])


def load_matrix(file):
    return np.load(file)


def calculates_ensemble_matrices(matrices, mode="median", max_distance=2):
    mean_sc_matrices = list()
    for matrix in matrices:
        matrix[matrix > max_distance] = np.nan
        if "proximity" in mode:
            mean_sc_matrix, n_cells = calculate_contact_probability_matrix(
                matrix,
                list(),
                1.0,
                norm="n_cells",
            )
        else:
            cells_to_plot = range(matrix.shape[2])
            mean_sc_matrix, _ = calculate_ensemble_pwd_matrix(
                matrix, 1.0, cells_to_plot, mode=mode
            )

        mean_sc_matrices.append(mean_sc_matrix)

    return mean_sc_matrices


def main_script(p):
    print("Processing files\n" + "=" * 16)

    mode = p["mode"]
    files = p["input_files"]
    max_distance = p["max_distance"]

    matrices = [load_matrix(file) for file in files]

    mean_sc_matrices = calculates_ensemble_matrices(
        matrices, mode=mode, max_distance=max_distance
    )

    [x, y] = [parses_matrix_to_vector(matrix) for matrix in mean_sc_matrices]

    r = calculates_pearson_correlation(x, y)

    print("Pearson Correlation Coefficient: ", r)

    plot_result(x, y, r, p)


# =============================================================================
# MAIN
# =============================================================================


def main():
    usage()

    # [parsing arguments]
    p = parse_arguments()

    # [loops over lists of datafolders]

    n_files = len(p["input_files"])
    print(f"> Number of input files: {n_files}")
    if n_files < 1:
        print("Please provide 2 input matrices")
        sys.exit()
    elif n_files > 2:
        print("Only two matrices can be processed at once.\n")
        sys.exit()

    print("Input files: ")
    for file in p["input_files"]:
        print(f"{file}")
    print("-" * 80)

    main_script(p)

    print("Finished execution")
    print("-" * 80)


if __name__ == "__main__":
    main()
