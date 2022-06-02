#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 09:24:51 2020

@author: marcnol
"""


import argparse
import csv
import json

#%% imports and plotting settings
import os

import matplotlib.gridspec as gridspec

# import matplotlib as plt
import matplotlib.pyplot as plt
import numpy as np

from matrixOperations.alignBarcodesMasks import plot_distance_histograms, plot_matrix

# import scaleogram as scg
from matrixOperations.HIMmatrixOperations import (
    AnalysisHiMMatrix,
    calculate_3_way_contact_matrix,
    get_multi_contact,
    plot_ensemble_3_way_contact_matrix,
    shuffle_matrix,
)

#%% define and loads datasets


def parse_arguments():
    # [parsing arguments]
    parser = argparse.ArgumentParser()
    parser.add_argument("-F1", "--rootFolder1", help="Folder with dataset 1")
    parser.add_argument("-F2", "--rootFolder2", help="Folder with dataset 2")
    parser.add_argument("-O", "--outputFolder", help="Folder for outputs")

    parser.add_argument(
        "-P",
        "--parameters",
        help="Provide name of parameter files. folders_to_load.json assumed as default",
    )
    parser.add_argument(
        "-A1", "--label1", help="Add name of label for dataset 1 (e.g. doc)"
    )
    parser.add_argument(
        "-W1",
        "--action1",
        help="Select: [all], [labeled] or [unlabeled] cells plotted for dataset 1 ",
    )
    parser.add_argument(
        "-A2", "--label2", help="Add name of label for dataset 1  (e.g. doc)"
    )
    parser.add_argument(
        "-W2",
        "--action2",
        help="Select: [all], [labeled] or [unlabeled] cells plotted for dataset 1 ",
    )
    parser.add_argument("--fontsize", help="Size of fonts to be used in matrix")
    parser.add_argument(
        "--axisLabel", help="Use if you want a label in x and y", action="store_true"
    )
    parser.add_argument(
        "--axisTicks", help="Use if you want axes ticks", action="store_true"
    )
    parser.add_argument(
        "--ratio",
        help="Does ratio between matrices. Default: difference",
        action="store_true",
    )
    parser.add_argument("--cAxis", help="absolute cAxis value for colormap")
    parser.add_argument(
        "--plottingFileExtension", help="By default: svg. Other options: pdf, png"
    )
    parser.add_argument(
        "--shuffle1",
        help="Provide shuffle vector: 0,1,2,3... of the same size or smaller than the original matrix. No spaces! comma-separated!",
    )
    parser.add_argument(
        "--shuffle2",
        help="Provide shuffle vector: 0,1,2,3... of the same size or smaller than the original matrix. No spaces! comma-separated!",
    )
    parser.add_argument(
        "--cMinMax",
        help="Provide min and max value for the colormap. Comma-separated, no spaces: 0,0.5 Overwrites --cAxis.",
    )

    args = parser.parse_args()

    run_parameters = {}
    run_parameters["pixelSize"] = 0.1

    if args.rootFolder1:
        rootFolder1 = args.rootFolder1
    else:
        rootFolder1 = "."

    if args.rootFolder2:
        rootFolder2 = args.rootFolder2
        run_parameters["run2Datasets"] = True
    else:
        rootFolder2 = "."
        run_parameters["run2Datasets"] = False

    if args.outputFolder:
        output_folder = args.outputFolder
    else:
        output_folder = "none"

    if args.parameters:
        run_parameters["parametersFileName"] = args.parameters
    else:
        run_parameters["parametersFileName"] = "folders_to_load.json"

    if args.label1:
        run_parameters["label1"] = args.label1
    else:
        run_parameters["label1"] = "doc"

    if args.label2:
        run_parameters["label2"] = args.label2
    else:
        run_parameters["label2"] = "NE"

    if args.action1:
        run_parameters["action1"] = args.action1
    else:
        run_parameters["action1"] = "labeled"

    if args.action2:
        run_parameters["action2"] = args.action2
    else:
        run_parameters["action2"] = "labeled"

    if args.fontsize:
        run_parameters["fontsize"] = args.fontsize
    else:
        run_parameters["fontsize"] = 12

    if args.axisLabel:
        run_parameters["axisLabel"] = args.axisLabel
    else:
        run_parameters["axisLabel"] = False

    if args.axisTicks:
        run_parameters["axisTicks"] = args.axisTicks
    else:
        run_parameters["axisTicks"] = False

    if args.ratio:
        run_parameters["ratio"] = args.ratio
    else:
        run_parameters["ratio"] = False

    if args.cAxis:
        run_parameters["cAxis"] = float(args.cAxis)
    else:
        run_parameters["cAxis"] = 0.6

    if args.plottingFileExtension:
        run_parameters["plottingFileExtension"] = "." + args.plottingFileExtension
    else:
        run_parameters["plottingFileExtension"] = ".svg"

    if args.shuffle1:
        run_parameters["shuffle1"] = args.shuffle1
    else:
        run_parameters["shuffle1"] = 0

    if args.shuffle2:
        run_parameters["shuffle2"] = args.shuffle2
    else:
        run_parameters["shuffle2"] = 0

    if args.cMinMax:
        run_parameters["cMinMax"] = args.cMinMax
    else:
        run_parameters["cMinMax"] = 0

    print("Input Folders:{}, {}".format(rootFolder1, rootFolder2))
    print("Input parameters:{}".format(run_parameters))

    return rootFolder1, rootFolder2, output_folder, run_parameters


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    rootFolder1, rootFolder2, output_folder, run_parameters = parse_arguments()

    him_data_1 = AnalysisHiMMatrix(run_parameters, rootFolder1)
    him_data_1.run_parameters["action"] = him_data_1.run_parameters["action1"]
    him_data_1.run_parameters["label"] = him_data_1.run_parameters["label1"]
    him_data_1.load_data()
    n_cells = him_data_1.n_cells_loaded()

    him_data_2 = AnalysisHiMMatrix(run_parameters, rootFolder2)
    him_data_2.run_parameters["action"] = him_data_2.run_parameters["action2"]
    him_data_2.run_parameters["label"] = him_data_2.run_parameters["label2"]
    him_data_2.load_data()
    n_cells_2 = him_data_2.n_cells_loaded()

    # cScale1 = him_data_1.data['ensembleContactProbability'].max() / run_parameters['cAxis']
    # cScale2 = him_data_2.data['ensembleContactProbability'].max() / run_parameters['scalingParameter']
    # print('scalingParameters={}'.format(run_parameters["scalingParameter"] ))

    if output_folder == "none":
        output_folder = him_data_1.data_folder

    outputFileName1 = (
        output_folder
        + os.sep
        + "Fig_ratio2HiMmatrices"
        + "_dataset1:"
        + him_data_1.dataset_name
        + "_label1:"
        + run_parameters["label1"]
        + "_action1:"
        + run_parameters["action1"]
        + "_dataset2:"
        + him_data_2.dataset_name
        + "_label2:"
        + run_parameters["label2"]
        + "_action2:"
        + run_parameters["action2"]
        + run_parameters["plottingFileExtension"]
    )

    outputFileName2 = (
        output_folder
        + os.sep
        + "Fig_mixedHiMmatrices"
        + "_dataset1:"
        + him_data_1.dataset_name
        + "_label1:"
        + run_parameters["label1"]
        + "_action1:"
        + run_parameters["action1"]
        + "_dataset2:"
        + him_data_2.dataset_name
        + "_label2:"
        + run_parameters["label2"]
        + "_action2:"
        + run_parameters["action2"]
        + run_parameters["plottingFileExtension"]
    )

    if (
        him_data_1.data["ensembleContactProbability"].shape
        == him_data_2.data["ensembleContactProbability"].shape
    ):
        ### Fig1: difference or ratio of the two matrices
        fig1 = plt.figure(constrained_layout=True)
        spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
        f_1 = fig1.add_subplot(spec1[0, 0])  # 16
        m1 = him_data_1.data["ensembleContactProbability"]
        m2 = him_data_2.data["ensembleContactProbability"]

        if run_parameters["shuffle1"] != 0:
            index1 = [int(i) for i in run_parameters["shuffle1"].split(",")]
            m1 = shuffle_matrix(m1, index1)

        if run_parameters["shuffle2"] != 0:
            index2 = [int(i) for i in run_parameters["shuffle2"].split(",")]
            m2 = shuffle_matrix(m2, index2)

        m1 = m1 / m1.max()
        m2 = m2 / m2.max()

        if run_parameters["ratio"] == True:
            matrix = np.log(m1 / m2)
            cmtitle = "log(ratio)"
        else:
            matrix = m1 - m2
            cmtitle = "difference"

        f1_ax1_im = him_data_1.plot_2d_matrix_simple(
            f_1,
            matrix,
            list(him_data_1.data["uniqueBarcodes"]),
            run_parameters["axisLabel"],
            run_parameters["axisLabel"],
            cmtitle=cmtitle,
            c_min=-run_parameters["cAxis"],
            c_max=run_parameters["cAxis"],
            fontsize=run_parameters["fontsize"],
            colorbar=True,
            axis_ticks=run_parameters["axisTicks"],
            c_m="RdBu",
        )
        plt.savefig(outputFileName1)
        print("Output figure: {}".format(outputFileName1))
        # plt.close()

        ### Fig2: "mixed matrix"
        fig2 = plt.figure(constrained_layout=True)
        spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
        f2 = fig2.add_subplot(spec2[0, 0])  # 16

        # load data once more
        matrix1 = him_data_1.data["ensembleContactProbability"]
        matrix2 = him_data_2.data["ensembleContactProbability"]

        if run_parameters["shuffle1"] != 0:
            index1 = [int(i) for i in run_parameters["shuffle1"].split(",")]
            matrix1 = shuffle_matrix(matrix1, index1)

        if run_parameters["shuffle2"] != 0:
            index2 = [int(i) for i in run_parameters["shuffle2"].split(",")]
            matrix2 = shuffle_matrix(matrix2, index2)

        for i in range(matrix1.shape[0]):
            for j in range(0, i):
                matrix1[i, j] = matrix2[i, j]

        if run_parameters["cMinMax"] == 0:
            c_min = 0
            c_max = run_parameters["cAxis"]
        else:
            index = [float(i) for i in run_parameters["cMinMax"].split(",")]
            c_min = index[0]
            c_max = index[1]

        him_data_1.plot_2d_matrix_simple(
            f2,
            matrix1,
            list(him_data_1.data["uniqueBarcodes"]),
            run_parameters["axisLabel"],
            run_parameters["axisLabel"],
            cmtitle="probability",
            c_min=c_min,
            c_max=c_max,
            fontsize=run_parameters["fontsize"],
            colorbar=True,
            axis_ticks=run_parameters["axisTicks"],
            c_m="coolwarm",
        )
        plt.savefig(outputFileName2)
        print("Output figure: {}".format(outputFileName2))

        # save also the npy
        outputFileName3 = outputFileName2.replace(
            run_parameters["plottingFileExtension"], ".npy"
        )
        np.save(outputFileName3, matrix1)
    else:
        print("Error: matrices do not have the same dimensions!")
