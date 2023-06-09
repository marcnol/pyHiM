#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 13:39:06 2023

@author: marcnol

trace_plot

script to plot one or multiple traces in 3D

Takes a trace file and either:
    - ranks traces and plots a selection
    - plots a user-selected trace in .ecsv (barcode, xyz) and PDF formats. The output files contain the trace name.
    - saves output coordinates for selected traces in pdb format so they can be loaded by other means
    including https://www.rcsb.org/3d-view, pymol, or nglviewer.
    
future:
    - options for 3D plotting in 3D: colors, backbone, background, etc
    - need to find how connectivites are made in the PDB to enforce them
    
installs:
    pip install nglview, pdbparser
    
example usage:
    
ls Trace_3D_barcode_KDtree_ROI:1.ecsv | trace_plot.py --selected_trace 5b1e6f89-0362-4312-a7ed-fc55ae98a0a5

this pipes the file 'Trace_3D_barcode_KDtree_ROI:1.ecsv' into trace_plot and then selects a trace for conversion.

format for json dict:
    {"12": "C  ", "18": "C  ", "9": "P  "}
keys provide barcode names in the trace file, these should be attributed to 3 character codes
"""


# =============================================================================
# IMPORTS
# =============================================================================q

import argparse
import csv
import glob
import json
import os
import select
import sys
from datetime import datetime
from astropy.io import ascii

import numpy as np

from imageProcessing.imageProcessing import Image
from matrixOperations.chromatin_trace_table import ChromatinTraceTable
from matrixOperations.HIMmatrixOperations import  write_xyz_2_pdb
from pdbparser.pdbparser import pdbparser

# =============================================================================
# FUNCTIONS
# =============================================================================q


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("-n", "--number_traces", help="Number of traces treated")
    parser.add_argument("-N","--N_barcodes", help="minimum_number_barcodes. Default = 2")
    parser.add_argument("--selected_trace", help="Selected trace for analysis")
    parser.add_argument("--barcode_type_dict", help="Json dictionnary linking barcodes and atom types (MUST BE 3 characters long!). ")

    p = {}

    args = parser.parse_args()
    if args.rootFolder:
        p["rootFolder"] = args.rootFolder
    else:
        p["rootFolder"] = "."

    if args.N_barcodes:
        p["N_barcodes"] = int(args.N_barcodes)
    else:
        p["N_barcodes"] = 2

    if args.number_traces:
        p["number_traces"] = int(args.number_traces)
    else:
        p["number_traces"] = 2

    if args.selected_trace:
        p["selected_trace"] = args.selected_trace
    else:
        p["selected_trace"] = 'fa9f0eb5-abcc-4730-bcc7-ba1da682d776'

    if args.barcode_type_dict:
        p["barcode_type_dict"] = args.barcode_type_dict
    else:
        p["barcode_type_dict"] = 'barcode_type_dict.json'
        
    p["trace_files"] = []
    if select.select([sys.stdin,], [], [], 0.0)[0]:
        p["trace_files"] = [line.rstrip("\n") for line in sys.stdin]
    else:
        print("Nothing in stdin. Please provide the list of files to treat using piping.")

    return p

def convert_trace_to_pdb(single_trace, export=None):

    #natoms = int(lines[0])
    # skip first 2 lines and build records
    records = []
    barcodes, X, Y, Z= single_trace["Barcode #"], single_trace["x"], single_trace["y"], single_trace["z"]

    center_of_mass = np.mean(X), np.mean(Y), np.mean(Z)
    unit_conversion = 10.0 # converts from nm to Angstroms

    for barcode,x,y,z in zip(barcodes, X,Y,Z):

        records.append( { "record_name"       : 'ATOM',
                          "serial_number"     : len(records)+1,
                          "atom_name"         : str(barcode),# 'C', 
                          "location_indicator": '',
                          "residue_name"      : 'XYZ',#str(barcode),
                          "chain_identifier"  : '',
                          "sequence_number"   : len(records)+1,
                          "code_of_insertion" : '',
                          "coordinates_x"     : unit_conversion*float(x - center_of_mass[0]),
                          "coordinates_y"     : unit_conversion*float(y - center_of_mass[1]),
                          "coordinates_z"     : unit_conversion*float(z - center_of_mass[2]),
                          "occupancy"         : 1.0,
                          "temperature_factor": 0.0,
                          "segment_identifier": '',
                          "element_symbol"    : str(barcode),
                          "charge"            : '',
                          } )
    # create and return pdb
    pdb = pdbparser(filePath = None)
    pdb.records = records
    # export
    if export is not None:
        pdb.export_pdb(export)
    # return
    return pdb

def runtime(folder, N_barcodes=2, trace_files=[], selected_trace = 'fa9f0eb5-abcc-4730-bcc7-ba1da682d776',barcode_type=dict()):

    # gets trace files

    if len(trace_files) > 0:

        print(
            "\n{} trace files to process= {}".format(
                len(trace_files), "\n".join(map(str, trace_files))
            )
        )

        # iterates over traces in folder
        for trace_file in trace_files:

            trace = ChromatinTraceTable()
            trace.initialize()

            # reads new trace
            trace.load(trace_file)

            # filters trace
            trace.filter_traces_by_n(minimum_number_barcodes=N_barcodes)

            # saves output trace
            outputfile = trace_file.rstrip(".ecsv") + "_selected_traces" + ".xyz"

            # indexes traces by Trace_ID
            trace_table = trace.data
            trace_table_indexed = trace_table.group_by("Trace_ID")
            print("$ number of traces to process: {}".format(len(trace_table_indexed)))

            # iterates over traces
            for idx, single_trace in enumerate(trace_table_indexed.groups):
                trace_id = single_trace["Trace_ID"][0]
                if trace_id == selected_trace:
                    print("Converting trace ID: {}".format(trace_id))
                    
                    # sorts by barcode
                    new_trace = single_trace.copy()
                    new_trace = new_trace.group_by('Barcode #')
                    ascii.write(new_trace['Barcode #', 'x','y','z'], selected_trace+'.ecsv', overwrite=True)  
 
                    #convert_trace_to_pdb(new_trace, export=selected_trace+'_v1.pdb')
                    write_xyz_2_pdb(selected_trace+'_v2.pdb', new_trace, barcode_type)
    else:
        print("No trace file found to process!")

    return len(trace_files)


# =============================================================================
# MAIN
# =============================================================================
def loads_barcode_type(p):
    import json

    filename = p["barcode_type_dict"]
    # Opening JSON file
    f = open(filename)
      
    # returns JSON object as a dictionary
    barcode_type = json.load(f)
      
    # Closing file
    f.close()

    print("$ {} barcode dictionary loaded")
    return barcode_type

def main():
    begin_time = datetime.now()

    # [parsing arguments]
    p = parse_arguments()
    # [loops over lists of datafolders]
    folder = p["rootFolder"]
    barcode_type = loads_barcode_type(p)
    
    n_traces_processed = runtime(
        folder, N_barcodes=p["N_barcodes"], trace_files=p["trace_files"],
        selected_trace=p["selected_trace"],
        barcode_type=barcode_type,
    )

    print(f"Processed <{n_traces_processed}> trace file(s)")
    print("Finished execution")


if __name__ == "__main__":
    main()



