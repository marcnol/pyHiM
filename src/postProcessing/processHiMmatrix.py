#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:36:20 2020

@author: marcnol


This script takes JSON file with folders where datasets are
stored and processes multiple PWD matrices together.

$ processHiMmatrix.py -F root_folder

outputs

sc_matrix_collated: 3D npy matrix. PWD matrix for single cells. Axes:0-1 barcodes, Axis:2, cellID
unique_barcodes: npy array. list of unique barcodes
SClabeledCollated: npy array. binary label indicating if cell is in pattern or not. Axis:0 cellID

"""

# =============================================================================
# IMPORTS
# =============================================================================q

import argparse
import csv
import json
import os
import sys
from datetime import datetime

import numpy as np

from core.data_file import save_json
from core.pyhim_logging import write_string_to_file
from matrixOperations.HIMmatrixOperations import (
    load_sc_data,
    load_sc_data_matlab,
    plot_ensemble_3_way_contact_matrix,
    plot_ensemble_contact_probability_matrix,
    plot_inverse_pwd_matrix,
    plot_single_contact_probability_matrix,
    plot_single_pwd_matrice,
)

# Olivier
csv.field_size_limit(sys.maxsize)

# =============================================================================
# FUNCTIONS
# =============================================================================q


def joinsListArrays(ListArrays, axis=0):
    joinedArray = np.zeros(0)
    for iArray in ListArrays:
        if joinedArray.shape[0] == 0:
            joinedArray = iArray
        else:
            joinedArray = np.concatenate((joinedArray, iArray), axis=axis)
    return joinedArray


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument(
        "-P",
        "--parameters",
        help="Provide name of parameter files. folders_to_load.json assumed as default",
    )
    parser.add_argument("-A", "--label", help="Add name of label (e.g. doc)")
    parser.add_argument(
        "-W", "--action", help="Select: [all], [labeled] or [unlabeled] cells plotted "
    )
    parser.add_argument(
        "--matlab", help="Use to load matlab formatted data", action="store_true"
    )
    parser.add_argument(
        "--saveMatrix", help="Use to load matlab formatted data", action="store_true"
    )
    parser.add_argument(
        "--getStructure", help="Use to save ShEc3D PDB structure", action="store_true"
    )
    parser.add_argument("--pixelSize", help="pixelSize in um")
    parser.add_argument(
        "--HiMnormalization",
        help="Normalization of contact matrix: nonNANs (default) or n_cells",
    )
    parser.add_argument("--d3", help="Use to load 3D maps", action="store_true")

    p = {}

    args = parser.parse_args()
    if args.rootFolder:
        p["rootFolder"] = args.rootFolder
    else:
        p["rootFolder"] = "."
        # p["rootFolder"] = "/home/marcnol/grey/docPaper_fullDatasets/updatedDatasets/white_wt_docTAD_nc14"

    if args.parameters:
        p["parametersFileName"] = args.parameters
    else:
        p["parametersFileName"] = "folders_to_load.json"

    if args.label:
        p["label"] = args.label
    else:
        p["label"] = "M"

    if args.action:
        p["action"] = args.action
    else:
        p["action"] = "all"

    if args.matlab:
        p["format"] = "matlab"
    else:
        p["format"] = "pyHiM"

    if args.saveMatrix:
        p["saveMatrix"] = True
    else:
        p["saveMatrix"] = False

    if args.getStructure:
        p["getStructure"] = True
    else:
        p["getStructure"] = False

    if args.pixelSize:
        p["pixelSize"] = float(args.pixelSize)
    else:
        p["pixelSize"] = 0.1

    if args.HiMnormalization:
        p["HiMnormalization"] = args.HiMnormalization
    else:
        p["HiMnormalization"] = "nonNANs"

    if args.d3:
        p["d3"] = args.d3
    else:
        p["d3"] = False

    return p


# =============================================================================
# MAIN
# =============================================================================


def main():
    begin_time = datetime.now()

    # [parsing arguments]
    p = parse_arguments()

    # [ initialises MD file]
    now = datetime.now()
    date_time = now.strftime("%d%m%Y_%H%M%S")
    filename_root = "processHiMmatrixAnalysis_"

    # [ Lists and loads datasets from different embryos]
    filename_list_data_json = p["rootFolder"] + os.sep + p["parametersFileName"]
    print("\n-------------------------------------------------------------------")
    if os.path.exists(filename_list_data_json):
        with open(filename_list_data_json, encoding="utf-8") as json_file:
            list_data = json.load(json_file)
        print(
            "Loaded JSON file with {} datasets from {}\n".format(
                len(list_data), filename_list_data_json
            )
        )
    else:
        print("File not found: {}".format(filename_list_data_json))
        sys.exit()

    # [ creates output folder]
    p["outputFolder"] = p["rootFolder"] + os.sep + "scHiMmatrices"
    if not os.path.exists(p["outputFolder"]):
        os.mkdir(p["outputFolder"])
        print("Folder created: {}".format(p["outputFolder"]))

    # [loops over lists of datafolders]
    for dataset_name in list(list_data.keys()):
        # [loads SC matrices]
        if p["format"] == "pyHiM":
            print(">>> Loading pyHiM-formatted dataset")
            (
                sc_matrix_collated,
                unique_barcodes,
                build_pwd_matrix_collated,
                run_name,
                sc_labeled_collated,
            ) = load_sc_data(list_data, dataset_name, p)
        elif p["format"] == "matlab":
            print(">>> Loading MATLAB-formatted dataset")
            (
                sc_matrix_collated,
                unique_barcodes,
                run_name,
                sc_labeled_collated,
            ) = load_sc_data_matlab(list_data, dataset_name, p)

        markdown_filename = (
            p["rootFolder"]
            + os.sep
            + filename_root
            + "_"
            + dataset_name
            + "_Cells:"
            + p["action"]
            + "_"
            + date_time
            + ".md"
        )
        write_string_to_file(
            markdown_filename, "# Post-processing of Hi-M matrices", "w"
        )
        write_string_to_file(
            markdown_filename,
            f"""**dataset: {dataset_name}** - **Cells: {p["action"]}**""",
            "a",
        )
        p["SClabeledCollated"] = sc_labeled_collated

        if len(sc_matrix_collated) > 0:
            # [plots distance matrix for each dataset]
            write_string_to_file(
                markdown_filename, "## single dataset PWD matrices", "a"
            )
            print(
                ">>> Producing {} PWD matrices for dataset {}\n".format(
                    len(sc_matrix_collated), dataset_name
                )
            )
            plot_single_pwd_matrice(
                sc_matrix_collated,
                unique_barcodes,
                run_name,
                list_data[dataset_name],
                p,
                markdown_filename,
                dataset_name=dataset_name,
            )

            write_string_to_file(
                markdown_filename, "## single dataset inverse PWD matrices", "a"
            )
            print(
                ">>> Producing {} inverse PWD matrices for dataset {}\n".format(
                    len(sc_matrix_collated), dataset_name
                )
            )
            plot_inverse_pwd_matrix(
                sc_matrix_collated,
                unique_barcodes,
                run_name,
                list_data[dataset_name],
                p,
                markdown_filename,
                dataset_name=dataset_name,
            )

            # [Plots contact probability matrices for each dataset]
            write_string_to_file(
                markdown_filename, "## single dataset Contact Probability matrices", "a"
            )
            print(
                ">>> Producing {} contact matrices for dataset {}\n".format(
                    len(sc_matrix_collated), dataset_name
                )
            )
            plot_single_contact_probability_matrix(
                sc_matrix_collated,
                unique_barcodes,
                run_name,
                list_data[dataset_name],
                p,
                markdown_filename=markdown_filename,
                dataset_name=dataset_name,
            )

            # [combines matrices from different embryos and calculates integrated contact probability matrix]
            write_string_to_file(
                markdown_filename, "## Ensemble contact probability", "a"
            )
            print(
                ">>> Producing ensemble contact matrix for dataset {}\n".format(
                    dataset_name
                )
            )
            (
                sc_matrixCollatedEnsemble,
                common_set_unique_barcodes,
            ) = plot_ensemble_contact_probability_matrix(
                sc_matrix_collated,
                unique_barcodes,
                run_name,
                list_data[dataset_name],
                p,
                markdown_filename=markdown_filename,
                dataset_name=dataset_name,
            )

            anchors = list_data[dataset_name]["3wayContacts_anchors"]
            anchors = [a - 1 for a in anchors]  # convert to zero-based
            s_out = "Probability"  # Probability or Counts
            write_string_to_file(markdown_filename, "## Ensemble 3way contacts", "a")
            print(
                ">>> Producing ensemble 3way contact matrix for dataset {}\n".format(
                    dataset_name
                )
            )
            plot_ensemble_3_way_contact_matrix(
                sc_matrix_collated,
                unique_barcodes,
                anchors,
                s_out,
                run_name,
                list_data[dataset_name],
                p,
                markdown_filename=markdown_filename,
                dataset_name=dataset_name,
            )

            # [deletes variables before starting new iteration]
            # del sc_matrix_collated, unique_barcodes, build_pwd_matrix_collated, run_name
            print("\nDone with dataset {}".format(dataset_name))
        else:
            print("\n Could not load ANY dataset!\n")

        # [saves output files]

        # creates output_filename root
        output_filename = (
            p["outputFolder"]
            + os.sep
            + dataset_name
            + "_label:"
            + p["label"]
            + "_action:"
            + p["action"]
        )

        # saves npy arrays
        if "SCmatrixCollatedEnsemble" in locals():
            np.save(
                output_filename + "_ensembleContactProbability.npy",
                sc_matrixCollatedEnsemble,
            )

        np.save(
            output_filename + "_SCmatrixCollated.npy",
            joinsListArrays(sc_matrix_collated, axis=2),
        )
        np.save(
            output_filename + "_SClabeledCollated.npy",
            joinsListArrays(sc_labeled_collated, axis=0),
        )

        # saves lists
        if "SCmatrixCollatedEnsemble" in locals():
            with open(
                output_filename + "_uniqueBarcodes.csv",
                mode="w",
                newline="",
                encoding="utf-8",
            ) as csvfile:
                spamwriter = csv.writer(
                    csvfile, delimiter=" ", quotechar="|", quoting=csv.QUOTE_MINIMAL
                )
                spamwriter.writerow(common_set_unique_barcodes)

        p["SClabeledCollated"] = []

        save_json(p, output_filename + "_parameters.json")

        with open(
            output_filename + "_runName.csv", mode="w", newline="", encoding="utf-8"
        ) as csvfile:
            spamwriter = csv.writer(
                csvfile, delimiter=" ", quotechar="|", quoting=csv.QUOTE_MINIMAL
            )
            spamwriter.writerow(run_name)

        print("Finished execution")


if __name__ == "__main__":
    main()
