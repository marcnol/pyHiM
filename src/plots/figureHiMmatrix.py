#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 09:04:10 2020

@author: marcnol
 Produces 
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

from matrixOperations.HIMmatrixOperations import (
    AnalysisHiMMatrix,
    calculate_ensemble_pwd_matrix,
    list_sc_to_keep,
    normalize_matrix,
    plot_distance_histograms,
    plot_matrix,
    plot_scalogram,
    shuffle_matrix,
)

# import scaleogram as scg


#%% define and loads datasets


def parse_arguments():
    # [parsing arguments]
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with dataset")
    parser.add_argument("-O", "--outputFolder", help="Folder for outputs")

    parser.add_argument(
        "-P",
        "--parameters",
        help="Provide name of parameter files. folders_to_load.json assumed as default",
    )
    parser.add_argument("-A", "--label", help="Add name of label (e.g. doc)")
    parser.add_argument(
        "-W", "--action", help="Select: [all], [labeled] or [unlabeled] cells plotted "
    )
    parser.add_argument("--fontsize", help="Size of fonts to be used in matrix")
    parser.add_argument(
        "--axisLabel", help="Use if you want a label in x and y", action="store_true"
    )
    parser.add_argument(
        "--axisTicks", help="Use if you want axes ticks", action="store_true"
    )
    parser.add_argument(
        "--barcodes",
        help="Use if you want barcode images to be displayed",
        action="store_true",
    )
    parser.add_argument(
        "--scalingParameter",
        help="Normalizing scaling parameter of colormap. Max will matrix.max()/scalingParameter",
    )
    parser.add_argument("--cScale", help="Colormap absolute scale")
    parser.add_argument(
        "--plottingFileExtension", help="By default: svg. Other options: pdf, png"
    )
    parser.add_argument(
        "--shuffle",
        help="Provide shuffle vector: 0,1,2,3... of the same size or smaller than the original matrix. No spaces! comma-separated!",
    )
    parser.add_argument(
        "--scalogram",
        help="Use if you want scalogram image to be displayed",
        action="store_true",
    )
    parser.add_argument(
        "--inputMatrix", help="contact, PWD, or iPWD. Default is contact"
    )
    parser.add_argument("--pixelSize", help="pixel size in um")
    parser.add_argument("--cmap", help="Colormap. Default: coolwarm")
    parser.add_argument(
        "--PWDmode",
        help="Mode used to calculate the mean distance. Can be either 'median' or KDE. Default: median",
    )

    args = parser.parse_args()

    run_parameters = {}

    if args.rootFolder:
        root_folder = args.rootFolder
    else:
        root_folder = "."
        # root_folder='/home/marcnol/data'+os.sep+'Experiment_18'

    if args.outputFolder:
        output_folder = args.outputFolder
    else:
        output_folder = "none"

    if args.parameters:
        run_parameters["parametersFileName"] = args.parameters
    else:
        run_parameters["parametersFileName"] = "folders_to_load.json"

    if args.label:
        run_parameters["label"] = args.label
    else:
        run_parameters["label"] = "doc"

    if args.action:
        run_parameters["action"] = args.action
    else:
        run_parameters["action"] = "labeled"

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

    if args.barcodes:
        run_parameters["barcodes"] = args.barcodes
    else:
        run_parameters["barcodes"] = False

    if args.scalingParameter:
        run_parameters["scalingParameter"] = float(args.scalingParameter)
    else:
        run_parameters["scalingParameter"] = 1.0

    if args.cScale:
        run_parameters["cScale"] = float(args.cScale)
    else:
        run_parameters["cScale"] = 0.0

    if args.plottingFileExtension:
        run_parameters["plottingFileExtension"] = "." + args.plottingFileExtension
    else:
        run_parameters["plottingFileExtension"] = ".png"

    if args.shuffle:
        run_parameters["shuffle"] = args.shuffle
    else:
        run_parameters["shuffle"] = 0

    if args.scalogram:
        run_parameters["scalogram"] = args.scalogram
    else:
        run_parameters["scalogram"] = False

    if args.inputMatrix:
        run_parameters["inputMatrix"] = args.inputMatrix
    else:
        run_parameters["inputMatrix"] = "contact"

    if args.pixelSize:
        run_parameters["pixelSize"] = args.pixelSize
    else:
        run_parameters["pixelSize"] = 0.1

    if args.cmap:
        run_parameters["cmap"] = args.cmap
    else:
        run_parameters["cmap"] = "coolwarm"

    if args.PWDmode:
        run_parameters["PWDmode"] = args.PWDmode
    else:
        run_parameters["PWDmode"] = "median"

    return root_folder, output_folder, run_parameters


#%%

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    print(">>> Producing HiM matrix")
    root_folder, output_folder, run_parameters = parse_arguments()

    him_data = AnalysisHiMMatrix(run_parameters, root_folder)

    him_data.load_data()

    n_cells = him_data.n_cells_loaded()

    him_data.retrieve_sc_matrix()

    if run_parameters["inputMatrix"] == "contact":
        # contact probability matrix
        matrix = him_data.data["ensembleContactProbability"]
        c_scale = matrix.max() / run_parameters["scalingParameter"]
    elif run_parameters["inputMatrix"] == "PWD":
        # PWD matrix
        # sc_matrix = him_data.sc_matrix_selected
        # print("SC matrix size: {}".format(sc_matrix.shape))
        # matrix, keep_plotting = calculate_ensemble_pwd_matrix(sc_matrix, run_parameters["pixelSize"],mode=run_parameters["PWDmode"])
        sc_matrix = him_data.sc_matrix_selected
        cells_to_plot = list_sc_to_keep(run_parameters, him_data.data["SClabeledCollated"])
        print("N cells to plot: {}/{}".format(len(cells_to_plot), sc_matrix.shape[2]))
        matrix, keep_plotting = calculate_ensemble_pwd_matrix(
            sc_matrix,
            run_parameters["pixelSize"],
            cells_to_plot,
            mode=run_parameters["PWDmode"],
        )
        if run_parameters["cScale"] == 0:
            c_scale = (
                matrix[~np.isnan(matrix)].max() / run_parameters["scalingParameter"]
            )
        else:
            c_scale = run_parameters["cScale"]
    elif run_parameters["inputMatrix"] == "iPWD":
        sc_matrix = him_data.sc_matrix_selected
        cells_to_plot = list_sc_to_keep(run_parameters, him_data.data["SClabeledCollated"])
        print("N cells to plot: {}/{}".format(len(cells_to_plot), sc_matrix.shape[2]))
        matrix, keep_plotting = calculate_ensemble_pwd_matrix(
            sc_matrix,
            run_parameters["pixelSize"],
            cells_to_plot,
            mode=run_parameters["PWDmode"],
        )
        matrix = np.reciprocal(matrix)
        c_scale = run_parameters["cScale"]

    print(
        "scalingParameters, scale={}, {}".format(
            run_parameters["scalingParameter"], c_scale
        )
    )

    n_cells = him_data.n_cells_loaded()

    n_datasets = len(him_data.data["runName"])

    if output_folder == "none":
        output_folder = him_data.data_folder

    output_filename = (
        output_folder
        + os.sep
        + "Fig_HiMmatrix"
        + "_dataset1:"
        + him_data.dataset_name
        + "_label:"
        + run_parameters["label"]
        + "_action:"
        + run_parameters["action"]
        + run_parameters["plottingFileExtension"]
    )

    if run_parameters["barcodes"]:
        fig1 = plt.figure(figsize=(10, 10), constrained_layout=False)
        gs1 = fig1.add_gridspec(
            nrows=19, ncols=22, left=0.05, right=0.95, wspace=0.05, hspace=0.05
        )
        f_1 = fig1.add_subplot(gs1[0:-1, 5:-1])
        f2 = fig1.add_subplot(gs1[:-1, 3], sharey=f_1)
        f3 = fig1.add_subplot(gs1[-1, 5:-1], sharex=f_1)
        ATACseqMatrix = (
            np.array(him_data.list_data[him_data.dataset_name]["BarcodeColormap"]) / 10
        )
        ATACseqMatrixV = np.copy(ATACseqMatrix).reshape((-1, 1))
        pos1 = f2.imshow(
            np.atleast_2d(ATACseqMatrixV), cmap="tab10"
        )  # colormaps RdBu seismic
        f2.set_xticklabels(())
        f2.set_yticklabels(())
        pos1.set_clim(vmin=-1, vmax=1)

        pos2 = f3.imshow(
            np.atleast_2d(ATACseqMatrix), cmap="tab10"
        )  # colormaps RdBu seismic
        f3.set_xticklabels(())
        f3.set_yticklabels(())
        pos2.set_clim(vmin=-1, vmax=1)

        barcodeLabels = np.arange(1, ATACseqMatrix.shape[0] + 1)
        for j in range(len(ATACseqMatrix)):
            text = f3.text(
                j,
                0,
                barcodeLabels[j],
                ha="center",
                va="center",
                color="w",
                fontsize=int((14.0 / 22.0) * float(run_parameters["fontsize"])),
            )
            text = f2.text(
                0,
                j,
                barcodeLabels[j],
                ha="center",
                va="center",
                color="w",
                fontsize=int((14.0 / 22.0) * float(run_parameters["fontsize"])),
            )

        colorbar = False
    else:
        fig1 = plt.figure(constrained_layout=True)
        spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)
        f_1 = fig1.add_subplot(spec1[0, 0])  # 16
        colorbar = True

    if run_parameters["shuffle"] == 0:
        index = range(matrix.shape[0])
    else:
        index = [int(i) for i in run_parameters["shuffle"].split(",")]
        matrix = shuffle_matrix(matrix, index)

    f1_ax1_im = him_data.plot_2d_matrix_simple(
        f_1,
        matrix,
        list(him_data.data["uniqueBarcodes"]),
        run_parameters["axisLabel"],
        run_parameters["axisLabel"],
        cmtitle="probability",
        c_min=0,
        c_max=c_scale,
        c_m=run_parameters["cmap"],
        fontsize=run_parameters["fontsize"],
        colorbar=colorbar,
        axis_ticks=run_parameters["axisTicks"],
        n_cells=n_cells,
        n_datasets=n_datasets,
        show_title=True,
    )

    # him_data.update_clims(0, c_scale, f_1)
    print("Output written to {}".format(output_filename))
    plt.savefig(output_filename)
    np.save(output_filename + ".npy", matrix)
    title_text = "N = {} | n = {}".format(n_cells, n_datasets)
    print("Title: {}".format(title_text))
    print("Output figure: {}".format(output_filename))

    # if run_parameters["scalogram"]:
    #     outputFileNameScalogram = (
    #         output_folder
    #         + os.sep
    #         + "Fig_HiMmatrix_scalogram"
    #         + "_dataset1:"
    #         + him_data.dataset_name
    #         + "_label:"
    #         + run_parameters["label"]
    #         + "_action:"
    #         + run_parameters["action"]
    #         + run_parameters["plottingFileExtension"]
    #     )

    #     plot_scalogram(matrix, outputFileNameScalogram)

    print("\nDone\n\n")
