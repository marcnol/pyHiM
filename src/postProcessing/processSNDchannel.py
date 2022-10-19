#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:21:01 2020

@author: marcnol

processSNDchannel


"""


# =============================================================================
# IMPORTS
# =============================================================================q

import glob
import os

try:
    del os.environ["MPLBACKEND"]  # unset before importing matplotlib
except:
    print("No environment variable MPLBACKEND found. Continuing anyway.")
import argparse
from datetime import datetime

import numpy as np
from astropy.table import Column, Table, vstack
from astropy.visualization import simple_norm
from matplotlib import pyplot as plt
from roipoly import MultiRoi

from fileProcessing.fileManagement import (
    Folders,
    Log,
    Parameters,
    Session,
    write_string_to_file,
)
from imageProcessing.imageProcessing import Image

# =============================================================================
# FUNCTIONS
# =============================================================================q


def show_image(data_2d, normalization="simple", size=(10, 10)):
    fig = plt.figure()
    fig.set_size_inches(size)

    ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()

    if normalization == "simple":
        norm = simple_norm(data_2d, "sqrt", percent=99.9)

    fig.add_axes(ax)
    ax.imshow(data_2d, origin="lower", cmap="Greys_r", norm=norm)
    ax.set_title("2D Data")
    return fig, ax


def createsUserMask(current_param, current_log, file_name, output_filename):

    # loads image

    # displays image
    im_obj = Image(current_param, current_log)
    im_obj.data_2d = np.load(file_name).squeeze()
    im_obj.show_image(show=True, normalization="simple")
    print("Click on the button to add a new ROI")

    # Draw multiple rois
    multiroi_named = MultiRoi(roi_names=["My first ROI", "My second ROI"])

    number_rois = len(multiroi_named.rois)
    print("Number of rois drawn: {}".format(number_rois))

    masks = np.zeros((im_obj.data_2d.shape[0], im_obj.data_2d.shape[1], number_rois))

    # Display image with all rois
    fig, ax = show_image(im_obj.data_2d)

    roi_names = []
    i = 0
    for name, roi in multiroi_named.rois.items():
        roi.display_roi()
        roi.display_mean(im_obj.data_2d)
        roi_names.append(name)
        masks[:, :, i] = roi.get_mask(im_obj.data_2d)
        i += 1
    plt.legend(roi_names, bbox_to_anchor=(1.2, 1.05))
    # plt.show()

    print("Saving and closing image with rois...")
    fig.savefig(output_filename + "_segmentedSources.png")
    plt.close(fig)

    # saves result
    np.save(output_filename, masks)


def processesUserMasks(current_param, current_log, processingList):

    if current_param.param_dict["segmentedObjects"]["operation"] == "overwrite":

        # session
        session_name = "processesUserMasks"

        # processes folders and files
        data_folder = Folders(current_param.param_dict["rootFolder"])
        data_folder.set_folders()
        current_log.add_simple_text(
            "\n===================={}====================\n".format(session_name)
        )
        current_log.report("folders read: {}".format(len(data_folder.list_folders)))

        position_roi_information = current_param.param_dict["acquisition"][
            "positionROIinformation"
        ]
        numberMaskedFiles = 0

        allresultsTable = Table()
        for current_folder in data_folder.list_folders:

            files_folder = glob.glob(current_folder + os.sep + "*.tif")
            data_folder.create_folders(current_folder, current_param)
            current_log.report("-------> Processing Folder: {}".format(current_folder))

            # generates lists of files to process
            current_param.find_files_to_process(files_folder)
            current_log.report(
                "About to process {} files\n".format(
                    len(current_param.files_to_process)
                )
            )

            if len(current_param.files_to_process) > 0:
                files_to_process = [
                    file
                    for file in glob.glob(
                        data_folder.output_folders["segmentedObjects"] + os.sep + "*"
                    )
                    if "SNDmask" in file.split("_")
                ]

                if processingList["cleanAllMasks"]:

                    # [clears all SND masks]
                    numberMaskedFiles = len(files_to_process)
                    for file_name in files_to_process:
                        os.remove(file_name)
                        current_log.report(
                            "Removing SND mask: {}".format(os.path.basename(file_name)),
                            "info",
                        )

                elif (
                    processingList["addMask"] == ""
                    and not processingList["cleanAllMasks"]
                ):

                    # [assigns cells to exsting masks]
                    resultsTable = assignsSNDmask2Cells(
                        files_to_process, position_roi_information
                    )
                    allresultsTable = vstack([allresultsTable, resultsTable])

                    current_log.report("assigning masks", "info")

                elif (
                    not processingList["cleanAllMasks"]
                    and len(processingList["addMask"]) > 0
                ):

                    # [makes new set of masks]
                    for file_name in current_param.files_to_process:

                        # gets filename information
                        ROI = os.path.basename(file_name).split("_")[
                            position_roi_information
                        ]

                        # checks that SND channel was projected and aligned
                        registeredFileName = (
                            data_folder.output_folders["alignImages"]
                            + os.sep
                            + os.path.basename(file_name).split(".")[0]
                            + "_2d_registered.npy"
                        )

                        if os.path.exists(registeredFileName):
                            output_filename = (
                                data_folder.output_folders["segmentedObjects"]
                                + os.sep
                                + os.path.basename(registeredFileName).split(".")[0]
                                + "_SNDmask"
                                + "_"
                                + processingList["addMask"]
                                + ".npy"
                            )
                            createsUserMask(
                                current_param,
                                current_log,
                                registeredFileName,
                                output_filename,
                            )
                            numberMaskedFiles += 1
                            current_log.report(
                                "Segmented image for SND channel ROI#{}: {}".format(
                                    ROI, output_filename
                                ),
                                "info",
                            )
                        else:
                            current_log.report(
                                "Could not find registered image for SND channel:{}".format(
                                    registeredFileName
                                ),
                                "error",
                            )

        tableOutputFileName = (
            data_folder.output_folders["segmentedObjects"]
            + os.sep
            + "snd_assigned_cells.ecsv"
        )
        print(f"\n>Number of cells found: {len(allresultsTable)}")
        print("\n>Saving output to: {}".format(tableOutputFileName))
        if os.path.exists(tableOutputFileName):
            print("\Results appended to {}".format(tableOutputFileName))
            previousresultsTable = Table.read(
                tableOutputFileName, format="ascii.ecsv"
            )  # ascii.ecsv
            allresultsTable = vstack([previousresultsTable, allresultsTable])
        if len(allresultsTable) > 0:
            print("\nTable written")
            allresultsTable.write(
                tableOutputFileName, format="ascii.ecsv", overwrite=True
            )

        return numberMaskedFiles
    else:
        return 0


def assignsSNDmask2Cells(files_to_process, position_roi_information):
    resultsTable = Table()

    numberFilesProcessed = 0
    # print(f"\nfiles2Process: {files_to_process}")

    for file_name in files_to_process:
        print(f"\n-----> Processing: {file_name}")
        ROI = os.path.basename(file_name).split("_")[position_roi_information]

        # [checks if DAPI mask exists for the file to process]
        fileNameDAPImask = (
            os.path.dirname(file_name)
            + os.sep
            + "_".join(os.path.basename(file_name).split("_")[0:7])
            + "_ch00_Masks.npy"
        )

        print(f"\nWill search masks in: {fileNameDAPImask}")

        if os.path.exists(fileNameDAPImask) and file_name.split(".")[-1] == "npy":

            # load DAPI mask
            maskDAPI = Image()
            maskDAPI.data_2d = np.load(fileNameDAPImask).squeeze()

            # load SND mask
            print(f"\nWill attemp to match masks from: {file_name}")
            maskSND = Image()
            maskSND.data_2d = np.load(file_name, allow_pickle=False).squeeze()

            # matches cells and SND masks
            new_matrix = maskDAPI.data_2d * maskSND.data_2d
            cellsWithinMask = np.unique(new_matrix)

            if len(cellsWithinMask) > 0:

                newTable = Table()
                colMask = Column(
                    [file_name.split("_")[-1].split(".")[0]] * len(cellsWithinMask),
                    name="MaskID #",
                    dtype=str,
                )
                col_roi = Column(
                    int(ROI) * np.ones(len(cellsWithinMask)), name="ROI #", dtype=int
                )
                col_cell_id = Column(cellsWithinMask, name="CellID #", dtype=int)
                newTable.add_column(col_roi, index=0)
                newTable.add_column(col_cell_id, index=1)
                newTable.add_column(colMask, index=2)

                resultsTable = vstack([resultsTable, newTable])
                del newTable

            # outputs list of cells within (1) and outside SND mask (0)
            current_log.report(
                "Matched DAPI and SND masks: \n SND: {}\n DAPI: {}\n".format(
                    file_name, fileNameDAPImask
                ),
                "info",
            )
            numberFilesProcessed += 1

        else:
            print(
                "\nERROR: Could not find expected DAPI mask: {}\n or file does not have ({}) the expected extension (.npy)".format(
                    fileNameDAPImask, file_name.split(".")[-1]
                )
            )

    return resultsTable


# =============================================================================
# MAIN
# =============================================================================

def main():
    begin_time = datetime.now()

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("-A", "--addMask", help="Add manual segmentation")
    parser.add_argument("--cleanAllMasks", help="Clear all masks", action="store_true")

    args = parser.parse_args()

    print(
        "\n--------------------------------------------------------------------------"
    )
    processingList = {}

    if args.rootFolder:
        root_folder = args.rootFolder
    else:
        root_folder = "."
        # root_folder='/home/marcnol/data/Experiment_4/0_Embryo/alignImages/'
        # fileNameRNA = root_folder+'scan_002_DAPI_001_ROI_converted_decon_ch01_2d_registered.npy'

    if args.addMask:
        processingList["addMask"] = args.addMask
    else:
        processingList["addMask"] = ""

    processingList["cleanAllMasks"] = args.cleanAllMasks

    print("parameters> root_folder: {}".format(root_folder))
    now = datetime.now()

    current_param = Parameters(root_folder=root_folder, file_name="infoList.json")
    labels = current_param.param_dict["labels"]

    session_name = "processSNDchannel"

    # setup logs
    current_log = Log(root_folder)
    current_log.add_simple_text(
        "\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format(
            session_name
        )
    )
    current_log.report(
        "Process SND channel MD: {}".format(current_log.markdown_filename)
    )
    write_string_to_file(
        current_log.markdown_filename,
        "# Process SND channel {}".format(now.strftime("%d/%m/%Y %H:%M:%S")),
        "w",
    )  # initialises MD file

    for label in labels:

        # sets parameters
        current_param = Parameters(
            root_folder=root_folder, label=label, file_name="infoList.json"
        )
        print("**Analyzing label: {}**".format(label))

        # processes Secondary masks
        if label == "RNA":
            numberMaskedFiles = processesUserMasks(
                current_param, current_log, processingList
            )
            if numberMaskedFiles > 0:
                print("Files processed: {}".format(numberMaskedFiles), "info")
            else:
                print("Nothing processed...", "warning")
        print("\n")
        del current_param

    # exits
    print("\n===================={}====================\n".format("Normal termination"))

    print("Elapsed time: {}".format(datetime.now() - begin_time))

if __name__ == "__main__":
    main()