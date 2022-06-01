#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:03:05 2020

@author: marcnol
"""


from scipy.io import loadmat
import os
import numpy as np
from astropy.table import Table, vstack

root_folder = "/home/marcnol/data/Experiment_Julian"


def loadsSCdataMATLAB(list_data, dataset_name, p):

    print("Dataset to load: {}\n\n".format(list(list_data.keys())[0]))

    sc_matrix_collated, unique_barcodes = [], []
    build_pwd_matrix_collated, run_name, SClabeledCollated = [], [], []

    for root_folder in list_data[dataset_name]["Folders"]:
        # [finds sample name]
        run_name.append(os.path.basename(os.path.dirname(root_folder)))
        # [loads and accumulates barcodes and scHiM matrix]
        filename_matrix = root_folder + os.sep + "HiMscMatrix.mat"
        filename_barcodes = root_folder + os.sep + "buildsPWDmatrix_uniqueBarcodes.ecsv"

        # loads SC matrix
        if os.path.exists(filename_matrix):
            data = loadmat(filename_matrix)
            sc_matrix_1 = data["distanceMatrixCumulative"]
            sc_matrix_collated.append(sc_matrix_1)
        else:
            print("*** Error: could not find {}".format(filename_matrix))

        # loads barcodes
        if os.path.exists(filename_barcodes):
            unique_barcodes.append(np.loadtxt(filename_barcodes).astype(int))
        else:
            print("*** Error: could not find {}".format(filename_barcodes))

        # loads cell attributes
        cell_attributes_matrix = data["cell_attributes_matrix"]
        results_table = cell_attributes_matrix[0, :]

        sc_labeled = np.zeros(len(results_table))
        index_cells_with_label = [i_row for i_row, row in enumerate(results_table) if row > 0]
        sc_labeled[index_cells_with_label] = 1
        SClabeledCollated.append(sc_labeled)

        print("\n>>>Merging root_folder: {}".format(root_folder))
        print("Cells added after merge: {}\n".format(sc_matrix_1.shape[2]))

    print("{} datasets loaded\n".format(len(sc_matrix_collated)))

    return (
        sc_matrix_collated,
        unique_barcodes,
        run_name,
        SClabeledCollated,
    )
