#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:07:05 2022

@author: marcnol

This script will load a trace file and analyze a number of properties such as:
    - number of barcodes detected per trace
    - number of duplicated barcodes
    - trace Rg

$ trace_analyzer.py

output:

trace_stats.csv

trace_ID, number of barcodes, number of duplications, Rg, 


"""

# =============================================================================
# IMPORTS
# =============================================================================q

import numpy as np
import os, sys
import json
from datetime import datetime
import argparse
import csv
import glob
import select

import collections
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from matrixOperations.chromatin_trace_table import ChromatinTraceTable
from imageProcessing.imageProcessing import Image

font = {'weight' : 'normal',
        'size'   : 30}

matplotlib.rc('font', **font)

# =============================================================================
# FUNCTIONS
# =============================================================================q


def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("--pipe", help="inputs Trace file list from stdin (pipe)", action = 'store_true')

    p = {}

    args = parser.parse_args()
    if args.rootFolder:
        p["rootFolder"] = args.rootFolder
    else:
        p["rootFolder"] = "."

    p["trace_files"] = list()
    if args.pipe:
        p["pipe"] = True
        if select.select([sys.stdin, ], [], [], 0.0)[0]:
            p["trace_files"] = [line.rstrip("\n") for line in sys.stdin]
        else:
            print("Nothing in stdin")
    else:
        p["pipe"] = False

    return p

def get_barcode_statistics(trace, output_filename='test.png'):
    '''
    Function that calculates the 
        - number of barcodes per trace
        - number of unique barcodes per trace
        - number of repeated barcodes per trace
    
    Parameters
    ----------
    trace : TYPE
        Trace table in ASTROPY Table format.
    output_filename : TYPE, optional
        Output figure in PNG. The default is 'test.png'.

    Returns
    -------
    None.

    '''
    trace_by_ID = trace.group_by('Trace_ID')

    trace_lengths = list()
    trace_unique_barcodes = list()
    trace_repeated_barcodes = list()
    number_unique_barcodes = list()
    number_repeated_barcodes = list()
    
    for sub_trace_table in trace_by_ID.groups:
        trace_lengths.append(len(sub_trace_table))
        
        unique_barcodes = list(set(sub_trace_table['Barcode #']))
        trace_unique_barcodes.append(unique_barcodes)
        number_unique_barcodes.append(len(unique_barcodes))
        
        repeated_barcodes = [item for item, count in collections.Counter(sub_trace_table['Barcode #']).items() if count > 1]
        trace_repeated_barcodes.append(repeated_barcodes)
        number_repeated_barcodes.append(len(repeated_barcodes))
        
    distributions = [trace_lengths,number_unique_barcodes,number_repeated_barcodes]
    axis_x_labels = ['number of barcodes','number of unique barcodes','number of repeated barcodes']
    number_plots = len(distributions)
    
    fig = plt.figure(constrained_layout=True)
    im_size = 10
    fig.set_size_inches(( im_size*number_plots,im_size))
    gs = fig.add_gridspec(1,number_plots)
    axes = [fig.add_subplot(gs[0,i]) for i in range(number_plots)]
    
    for axis, distribution,xlabel in zip(axes,distributions,axis_x_labels):
        axis.hist(distribution, alpha=.3)
        axis.set_xlabel(xlabel)
        axis.set_ylabel('counts')
        axis.set_title('n = '+str(len(distribution))+' | median = '+str(np.median(distribution)))
        
    plt.savefig(output_filename)
    

def analyze_trace(trace, trace_file):
    '''
    Launcher function that will perform different kinds of trace analyses

    Parameters
    ----------
    trace : ChromatinTraceTable Class
        Trace table, instance of the ChromatinTraceTable Class.
    trace_file : string
        file name of trace table in ecsv format.

    Returns
    -------
    None.

    '''
    
    trace_table = trace.data

    print(f"$ Number of lines in trace: {len(trace_table)}")

    output_filename = [trace_file.split(".")[0], "_trace_statistics", ".png"]
    
    get_barcode_statistics(trace_table, "".join(output_filename))


def process_traces(folder, trace_files = list()):
    '''
    Processes list of trace files and sends each to get analyzed individually

    Parameters
    ----------
    folder : TYPE
        DESCRIPTION.
    trace_files : TYPE, optional
        DESCRIPTION. The default is list().

    Returns
    -------
    None.

    '''
    trace_folder = folder.rstrip("/") + os.sep + "buildsPWDmatrix" + os.sep

    if len(trace_files)<1:
        trace_files = [x for x in glob.glob(trace_folder + "Trace*ecsv") if "uniqueBarcodes" not in x]

    # removes already labeled trace files
    trace_files = [x for x in trace_files if "labeled" not in x]

    print("\n{} trace files to process= {}".format(len(trace_files), "\n".join(map(str, trace_files))))

    if len(trace_files) > 0:
        # iterates over traces in folder
        for trace_file in trace_files:

            trace = ChromatinTraceTable()
            trace.initialize()

            # reads new trace
            trace.load(trace_file)
            
            print(f"> Analyzing traces for {trace_file}")            

            print(f"> Plotting traces for {trace_file}")
            trace.plots_traces([trace_file.split(".")[0], "_traces_XYZ", ".png"])

            analyze_trace(trace,trace_file)

            # outputfile = trace_file.rstrip(".ecsv") + "_labeled" + ".ecsv"


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    # [parsing arguments]
    p = parseArguments()
    
    # [loops over lists of datafolders]
    folder = p["rootFolder"]
    process_traces(folder, trace_files = p["trace_files"])

    print("Finished execution")
