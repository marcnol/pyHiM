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

import numpy as np

import scipy.stats
import matplotlib.pyplot as plt     
import seaborn as sns

# =============================================================================
# FUNCTIONS
# =============================================================================q

def usage():
    print("-"*80)
    print("Example usage:\n")
    print("$ ls buildsPWDmatrix/*HiM*npy")
    print("buildsPWDmatrix/buildsPWDmatrix_DAPI_HiMscMatrix.npy  buildsPWDmatrix/buildsPWDmatrix_mask0_HiMscMatrix.npy")
    print("$ ls buildsPWDmatrix/*HiM*npy | compare_PWD_matrices.py\n")
    print("You should now see a file called scatter_plot.png\n")
    print("-"*80)

def parse_arguments():
    parser = argparse.ArgumentParser(
        "This script will compare two PWD matrices. Please provide input matrices by using piping.\n Example: $ ls matrix1.npy matrix2.npy | compare_PWD_matrices"
    )

    parser.add_argument("-F", "--rootFolder", help="Folder with images")

    p = {}

    args = parser.parse_args()
    if args.rootFolder:
        p["rootFolder"] = args.rootFolder
    else:
        p["rootFolder"] = "."

    p["input_files"] = list()

    if select.select([sys.stdin,], [], [], 0.0)[0]:
        p["input_files"] = [line.rstrip("\n") for line in sys.stdin]
    else:
        print("Nothing in stdin. Please provide list of localization files to process.")

    print("Input parameters\n" + "=" * 16)
    for item in p.keys():
        print("{}-->{}".format(item, p[item]))

    return p

def plot_result(x,y,p):

    plt.figure(figsize=(15,15))
    plt.rcParams.update({'font.size': 20})

    plt.scatter(x, y, label = "Pearson = {:.2f}".format(p), color = "m", 
                marker = "o", s = 30, linewidth = 5) 
      
    plt.xscale('log') 
    plt.yscale('log') 
      
    plt.title('Correlation between PWD distances')
    plt.xlabel('dataset 1') 
    plt.ylabel('dataset 2') 
    plt.legend() 
    
    axisX = np.linspace(np.min([x[x!=0],y[y!=0]]), np.max([x,y]),100)
    plt.plot(axisX, axisX, color='red', linewidth=3)        
    
    #plt.show()
    filename = 'scatter_plot.png'
    plt.savefig(filename)
    print(f"> Output image saved as : {filename}")
    
def calculates_pearson_correlation(x, y):
     
    r, p = scipy.stats.pearsonr(x, y)
    
    return r

def parses_matrix_to_vector(matrix):
    matrix_shape = matrix.shape
    matrix = np.nan_to_num(matrix)
    return matrix.reshape(matrix_shape[0]*matrix_shape[1])

def load_matrix(file):
    return np.load(file)

def main_script(p):
    print("Processing files\n" + "=" * 16)

    files = p["input_files"]
    matrices = [load_matrix(file) for file in files]

    matrices = [np.nanmean(matrix,axis=2) for matrix in matrices]

    [x,y] = [parses_matrix_to_vector(matrix) for matrix in matrices]

    r = calculates_pearson_correlation(x, y)

    print("Pearson Correlation Coefficient: ", r)
    
    plot_result(x,y, r)
    
# =============================================================================
# MAIN
# =============================================================================


def main():

    usage()

    # [parsing arguments]
    p = parse_arguments()

    # [loops over lists of datafolders]
    folder = p["rootFolder"]
    
    n_files = len(p["input_files"])
    print(f"> Number of input files: {n_files}")
    if n_files < 1:
        print(
            "Please provide matrices by piping. \n Example: \n$ ls matrix1.npy matrix2.npy | compare_PWD_matrices.\n"
        )
        sys.exit()
    elif n_files > 2:
        print("Only two matrices can be processed at once.\n")
        sys.exit()

    print("Input files: ")
    for file in p["input_files"]:
        print(f"{file}")
    print("-"*80)
    
    main_script(p)

    print("Finished execution")
    print("-"*80)


if __name__ == "__main__":
    main()
