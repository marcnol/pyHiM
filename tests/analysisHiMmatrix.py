#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 15:36:31 2020

@author: marcnol

Figure 1

"""

#%% imports and plotting settings
import os
import numpy as np
import argparse

# import matplotlib as plt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json, csv
from alignBarcodesMasks import plot_distance_histograms, plot_matrix
import scaleogram as scg
from HIMmatrixOperations import plot_ensemble_3_way_contact_matrix, calculate_3_way_contact_matrix, get_multi_contact


# font = {'family' : 'normal',
#         'weight' : 'normal',
#         'size'   : 18}

# plt.rc('font', **font)

plottingFileExtension = ".pdf"


class AnalysisHiMMatrix:
    def __init__(self, run_parameters, root_folder="."):
        self.data_folder = root_folder + os.sep + "scHiMmatrices"
        self.run_parameters = run_parameters
        self.root_folder = root_folder

    def load_data(self):

        # loads datasets
        filename_list_data_json = self.root_folder + os.sep + self.run_parameters["parametersFileName"]

        with open(filename_list_data_json, encoding="utf-8") as json_file:
            list_data = json.load(json_file)

        dataset_name = list(list_data.keys())[0]
        print("Dataset: {}".format(dataset_name))

        output_filename = self.data_folder + os.sep + dataset_name + "_Cells:" + self.run_parameters["action"]

        filename_parameters_json = output_filename + "_parameters.json"
        with open(filename_parameters_json, encoding="utf-8") as json_file:
            p = json.load(json_file)

        data_files = {}
        data_files["ensembleContactProbability"] = "_ensembleContactProbability.npy"
        # data_files['unique_barcodes']="_uniqueBarcodes.npy"
        data_files["sc_matrix_collated"] = "_SCmatrixCollated.npy"
        data_files["SClabeledCollated"] = "_SClabeledCollated.npy"

        data_files["anchor:4"] = "_anchor:4_ensemble3wayContacts.npy"
        data_files["anchor:6"] = "_anchor:6_ensemble3wayContacts.npy"
        data_files["anchor:9"] = "_anchor:9_ensemble3wayContacts.npy"
        data_files["anchor:10"] = "_anchor:10_ensemble3wayContacts.npy"
        data_files["anchor:13"] = "_anchor:13_ensemble3wayContacts.npy"
        data_files["anchor:16"] = "_anchor:16_ensemble3wayContacts.npy"

        data = {}

        for i_data_file in data_files.keys():
            print("Loaded: {}".format(i_data_file))
            data[i_data_file] = np.load(output_filename + data_files[i_data_file]).squeeze()

        run_name = load_list(output_filename + "_runName.csv")
        data["unique_barcodes"] = load_list(output_filename + "_uniqueBarcodes.csv")

        print("Loaded: {}".format("run_name"))
        data["run_name"] = run_name

        self.data = data
        self.data_files = data_files

    # functions

    def plot_2d_matrix_simple(
        self,
        ifigure,
        matrix,
        unique_barcodes,
        yticks,
        xticks,
        cmtitle="probability",
        c_min=0,
        c_max=1,
        c_m="coolwarm",
        fontsize=12,
        colorbar=False,
        axis_ticks=False,
    ):

        pos = ifigure.imshow(matrix, cmap=c_m)  # colormaps RdBu seismic
        # plots figure
        if xticks:
            ifigure.set_xlabel("barcode #", fontsize=fontsize)
            if not axis_ticks:
                ifigure.set_xticklabels(())
            else:
                print("barcodes:{}".format(unique_barcodes))
                # ifigure.set_xticks(np.arange(matrix.shape[0]),unique_barcodes)
                ifigure.set_xticklabels(unique_barcodes)

        else:
            ifigure.set_xticklabels(())
        if yticks:
            ifigure.set_ylabel("barcode #", fontsize=fontsize)
            if not axis_ticks:
                ifigure.set_yticklabels(())
            else:
                # ifigure.set_yticks(np.arange(matrix.shape[0]), unique_barcodes)
                ifigure.set_yticklabels(unique_barcodes)
        else:
            ifigure.set_yticklabels(())

        for xtick, ytick in zip(ifigure.xaxis.get_majorticklabels(), ifigure.yaxis.get_majorticklabels()):
            xtick.set_fontsize(fontsize)
            ytick.set_fontsize(fontsize)

        if colorbar:
            cbar = plt.colorbar(pos, fraction=0.046, pad=0.04)
            cbar.minorticks_on()
            cbar.set_label(cmtitle)
            pos.set_clim(vmin=c_min, vmax=c_max)
        return pos

    def update_clims(self, c_min, c_max, axes):
        for ax in axes:
            ax.set_clim(vmin=c_min, vmax=c_max)


def load_list(file_name):
    with open(file_name, newline="", encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ", quotechar="|")
        run_name = []
        for row in spamreader:
            print(", ".join(row))
            if len(run_name) > 0:
                run_name.append(row)
            else:
                run_name = row

    return run_name
