#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 21:04:10 2023

@author: marcnol

script to compare PWD matrices from two experiments
- pearson correlation in ensemble
- map barcodes to allow for experiments with different barcode combinations
- same but single cell


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
        p["localization_files"] = [line.rstrip("\n") for line in sys.stdin]
    else:
        print("Nothing in stdin. Please provide list of localization files to process.")

    print("Input parameters\n" + "-" * 15)
    for item in p.keys():
        print("{}-->{}".format(item, p[item]))

    return p

def plot_result_sns(x,y, p):

    # Set figure size
    plt.figure(figsize = (20,20))
    
    # Plot scatterplot
    sns.scatterplot(x, y, logx = True, logy = True)
    
    # Show plot
    plt.show()

def plot_result(x,y, p):

    plt.figure(figsize=(20,20))
    plt.scatter(x, y, label = "Pearson = {}".format(p), color = "m", 
                marker = "o", s = 30, linewidth = 5) 
      
    plt.xscale('log') 
    plt.yscale('log') 
      
    plt.title('Scatterplot in log log') 
    plt.xlabel('x-axis') 
    plt.ylabel('y-axis') 
    plt.legend() 
      
    plt.show()
    
def calculates_pearson_correlation(x, y):
   
    x = [1, 2, 3, 4, 5]
    y = [2, 4, 5, 4, 5]
    
    r, p = scipy.stats.pearsonr(x, y)
    
    return r

def parses_matrix_to_vector(matrix):
    
    return np.flatten(matrix)

def load_matrix(file):
    return np.load(file)

def main_script(p):

    files = p["input_files"]

    matrices = [load_matrix(file) for file in files]

    [x,y] = [parses_matrix_to_vector(matrix) for matrix in matrices]

    r = calculates_pearson_correlation(x, y)

    print("Pearson Correlation Coefficient: ", r)
    
    plot_result(x,y, p)

# =============================================================================
# MAIN
# =============================================================================


def main():

    # [parsing arguments]
    p = parse_arguments()

    # [loops over lists of datafolders]
    folder = p["rootFolder"]

    if len(p["input_files"]) < 1:
        print(
            "Please provide matrices by piping. \n Example: \n$ ls matrix1.npy matrix2.npy | compare_PWD_matrices.\n"
        )
        sys.exit()
    elif len(p["input_files"]) > 2:
        print("Only two matrices can be processed at once.\n")
        sys.exit()

    print("Input files:\n ")
    for file in p["input_files"]:
        print(f"{file}\n")

    main_script(p)

    print("Finished execution")


if __name__ == "__main__":
    main()
