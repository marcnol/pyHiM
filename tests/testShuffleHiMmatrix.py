#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 20:48:16 2021

@author: marcnol
"""

import numpy as np
import matplotlib.pylab as plt

def shuffle_matrix(matrix, index):

    new_size = len(index)
    new_matrix = np.zeros((new_size, new_size))

    if new_size <= matrix.shape[0]:
        for i in range(new_size):
            for j in range(new_size):
                if index[i] < matrix.shape[0] and index[j] < matrix.shape[0]:
                    new_matrix[i, j] = matrix[index[i], index[j]]
    else:
        print("Error: shuffle size {} is larger than matrix dimensions {}".format(new_size, matrix.shape[0]))
        print("Shuffle: {} ".format(index))

    return new_matrix

def calculate_contact_probability_matrix(i_sc_matrix_collated, i_unique_barcodes, pixel_size, threshold=0.25, norm="n_cells"):

    n_x = n_y = i_sc_matrix_collated.shape[0]
    n_cells = i_sc_matrix_collated.shape[2]
    sc_matrix = np.zeros((n_x, n_y))

    for i in range(n_x):
        for j in range(n_y):
            if i != j:
                distance_distribution = pixel_size * i_sc_matrix_collated[i, j, :]

                # normalizes # of contacts by the # of cells
                if norm == "n_cells":
                    probability = len(np.nonzero(distance_distribution < threshold)[0]) / n_cells
                    # print('Using n_cells normalisation')

                # normalizes # of contacts by the # of PWD detected in each bin
                elif norm == "nonNANs":
                    number_nans = len(np.nonzero(np.isnan(distance_distribution))[0])
                    if n_cells == number_nans:
                        probability = 1
                    else:
                        probability = len(np.nonzero(distance_distribution < threshold)[0]) / (n_cells - number_nans)
                    # print('Using NonNANs normalisation {}'.format(n_cells-number_nans))
                sc_matrix[i, j] = probability

    return sc_matrix, n_cells

file_name="/home/marcnol/data/Embryo_debug_dataset/run_markus/buildsPWDmatrix/buildsPWDmatrix_3D_HiMscMatrix.npy"
# file_name="/home/marcnol/data/Embryo_debug_dataset/run_markus_raw_3D/buildsPWDmatrix/buildsPWDmatrix_3D_HiMscMatrix.npy"
# file_name='/home/marcnol/data/Embryo_debug_dataset/run_markus_raw_3D/buildsPWDmatrix/buildsPWDmatrix_3D_HiMscMatrix.npy'

shuffle = [0, 7, 1, 8, 2, 9, 10, 17, 11, 18, 12, 19, 13, 3, 14, 4, 15, 5, 16, 6]

matrixSC =np.load(file_name)
pixel_size = 0.1
matrix,_ = calculate_contact_probability_matrix(matrixSC, shuffle, pixel_size, threshold=0.25, norm="n_cells")

# matrix = np.nanmean(matrixSC,axis=2)
matrix_shuffled = shuffle_matrix(matrix, shuffle)
plt.imshow(matrix_shuffled,cmap='coolwarm',vmax=.01)