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

# import scaleogram as scg
from matrixOperations.HIMmatrixOperations import (
    AnalysisHiMMatrix,
    calculate_3_way_contact_matrix,
    calculate_ensemble_pwd_matrix,
    get_multi_contact,
    list_sc_to_keep,
    plot_distance_histograms,
    plot_matrix,
    plot_ensemble_3_way_contact_matrix,
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
        "--normalize",
        help="Matrix normalization factor: maximum, none, single value (normalize 2nd matrix by), bin pair e.g. 1,2",
    )
    parser.add_argument(
        "--inputMatrix", help="Source of input matrix: contact (default), PWD, iPWD"
    )
    parser.add_argument("--pixelSize", help="pixelSize in microns")

    run_parameters = {}

    args = parser.parse_args()

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

    if args.pixelSize:
        run_parameters["pixelSize"] = args.pixelSize
    else:
        run_parameters["pixelSize"] = 0.1

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
        run_parameters["cAxis"] = [float(i) for i in args.cAxis.split(",")]
    else:
        run_parameters["cAxis"] = 0.6

    if args.plottingFileExtension:
        run_parameters["plottingFileExtension"] = "." + args.plottingFileExtension
    else:
        run_parameters["plottingFileExtension"] = ".svg"

    if args.normalize:
        run_parameters["normalize"] = args.normalize
    else:
        run_parameters["normalize"] = "none"

    if args.inputMatrix:
        run_parameters["inputMatrix"] = args.inputMatrix
    else:
        run_parameters["inputMatrix"] = "contact"

    print("Input Folders:{}, {}".format(rootFolder1, rootFolder2))
    print("Input parameters:{}".format(run_parameters))

    return rootFolder1, rootFolder2, output_folder, run_parameters


def normalize_matrix(m1, m2, mode):

    print("Normalization: {}".format(mode))

    if "maximum" in mode:  # normalizes by maximum
        m1_norm = m1.max()
        m2_norm = m2.max()
    elif len(mode.split(",")) > 1:  # normalizes by bin
        N = mode.split(",")
        m1_norm = 1
        m2_norm = m2[int(N[0]), int(N[1])] / m1[int(N[0]), int(N[1])]
    elif "none" in mode:  # no normalization
        m1_norm = 1
        m2_norm = 1
    else:  # normalizes by given factor
        norm_factor = float(mode)
        m1_norm = 1
        m2_norm = norm_factor

    print("Normalizations: m1= {} | m2={}".format(m1_norm, m2_norm))

    m1 = m1 / m1_norm
    m2 = m2 / m2_norm

    return m1, m2


# =============================================================================
# MAIN
# =============================================================================

def main():

    rootFolder1, rootFolder2, output_folder, run_parameters = parse_arguments()

    him_data_1 = AnalysisHiMMatrix(run_parameters, rootFolder1)
    him_data_1.run_parameters["action"] = him_data_1.run_parameters["action1"]
    him_data_1.run_parameters["label"] = him_data_1.run_parameters["label1"]
    him_data_1.load_data()
    n_cells = him_data_1.n_cells_loaded()
    him_data_1.retrieve_sc_matrix()

    him_data_2 = AnalysisHiMMatrix(run_parameters, rootFolder2)
    him_data_2.run_parameters["action"] = him_data_2.run_parameters["action2"]
    him_data_2.run_parameters["label"] = him_data_2.run_parameters["label2"]
    him_data_2.load_data()
    n_cells_2 = him_data_2.n_cells_loaded()
    him_data_2.retrieve_sc_matrix()

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
    )

    if "contact" in run_parameters["inputMatrix"]:
        m1 = him_data_1.data["ensembleContactProbability"]
        m2 = him_data_2.data["ensembleContactProbability"]
    elif "iPWD" in run_parameters["inputMatrix"]:
        m1 = him_data_1.sc_matrix_selected
        m2 = him_data_2.sc_matrix_selected
        cells2Plot1 = list_sc_to_keep(run_parameters, him_data_1.data["SClabeledCollated"])
        cells2Plot2 = list_sc_to_keep(run_parameters, him_data_2.data["SClabeledCollated"])
        dataset1 = list(him_data_1.list_data.keys())[0]
        dataset2 = list(him_data_2.list_data.keys())[0]
        m1, _ = calculate_ensemble_pwd_matrix(
            m1,
            run_parameters["pixelSize"],
            cells2Plot1,
            mode=him_data_1.list_data[dataset1]["PWD_mode"],
        )
        m2, _ = calculate_ensemble_pwd_matrix(
            m2,
            run_parameters["pixelSize"],
            cells2Plot2,
            mode=him_data_2.list_data[dataset2]["PWD_mode"],
        )
        m1 = np.reciprocal(m1)
        m2 = np.reciprocal(m2)
    elif "PWD" in run_parameters["inputMatrix"]:
        m1 = him_data_1.sc_matrix_selected
        m2 = him_data_2.sc_matrix_selected
        cells2Plot1 = list_sc_to_keep(run_parameters, him_data_1.data["SClabeledCollated"])
        cells2Plot2 = list_sc_to_keep(run_parameters, him_data_2.data["SClabeledCollated"])
        dataset1 = list(him_data_1.list_data.keys())[0]
        dataset2 = list(him_data_2.list_data.keys())[0]
        m1, _ = calculate_ensemble_pwd_matrix(
            m1,
            run_parameters["pixelSize"],
            cells2Plot1,
            mode=him_data_1.list_data[dataset1]["PWD_mode"],
        )
        m2, _ = calculate_ensemble_pwd_matrix(
            m2,
            run_parameters["pixelSize"],
            cells2Plot2,
            mode=him_data_2.list_data[dataset2]["PWD_mode"],
        )

    if m1.shape == m2.shape:

        fig1 = plt.figure(constrained_layout=True)
        spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
        f_1 = fig1.add_subplot(spec1[0, 0])  # 16

        if "none" in run_parameters["normalize"]:  # sets default operation
            mode = "maximum"
        else:
            mode = run_parameters["normalize"]

        _m1, _m2 = m1.copy(), m2.copy()

        _m1, _m2 = normalize_matrix(_m1, _m2, mode)

        if run_parameters["ratio"] == True:
            matrix = np.log(_m1 / _m2)
            cmtitle = "log(ratio)"
        else:
            matrix = _m1 - _m2
            cmtitle = "difference"

        if len(run_parameters["cAxis"]) == 2:
            c_scale = run_parameters["cAxis"][1]
        else:
            c_scale = run_parameters["cAxis"][0]
        print("Clim used: {}\n".format(c_scale))

        f1_ax1_im = him_data_1.plot_2d_matrix_simple(
            f_1,
            matrix,
            list(him_data_1.data["uniqueBarcodes"]),
            run_parameters["axisLabel"],
            run_parameters["axisLabel"],
            cmtitle=cmtitle,
            c_min=-c_scale,
            c_max=c_scale,
            fontsize=run_parameters["fontsize"],
            colorbar=True,
            axis_ticks=run_parameters["axisTicks"],
            c_m="RdBu",
        )
        plt.savefig(outputFileName1)
        print("Output figure: {}".format(outputFileName1))
        # plt.close()

        # plots mixed matrix
        fig2 = plt.figure(constrained_layout=True)
        spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
        f2 = fig2.add_subplot(spec2[0, 0])  # 16

        # plots mixed matrix
        _m1, _m2 = m1.copy(), m2.copy()

        if "none" in run_parameters["normalize"]:  # sets default operation
            mode = "none"
        else:
            mode = run_parameters["normalize"]
        _m1, _m2 = normalize_matrix(_m1, _m2, mode)
        matrix2 = _m1

        for i in range(matrix2.shape[0]):
            for j in range(0, i):
                matrix2[i, j] = _m2[i, j]

        him_data_1.plot_2d_matrix_simple(
            f2,
            matrix2,
            list(him_data_1.data["uniqueBarcodes"]),
            run_parameters["axisLabel"],
            run_parameters["axisLabel"],
            cmtitle="probability",
            c_min=0,
            c_max=run_parameters["cAxis"][0],
            fontsize=run_parameters["fontsize"],
            colorbar=True,
            axis_ticks=run_parameters["axisTicks"],
            c_m="coolwarm",
        )
        plt.savefig(outputFileName2 + run_parameters["plottingFileExtension"])
        np.save(outputFileName2 + ".npy", matrix2)
        print(
            "Output figure: {}".format(
                outputFileName2 + run_parameters["plottingFileExtension"]
            )
        )

    else:
        print("Error: matrices do not have the same dimensions!")

if __name__ == "__main__":
    main()