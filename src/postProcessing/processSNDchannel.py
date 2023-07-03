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

import argparse
import glob
import os
from datetime import datetime

import numpy as np
from astropy.table import Column, Table, vstack
from astropy.visualization import simple_norm
from matplotlib import pyplot as plt
from roipoly import MultiRoi

from core.data_manager import DataManager
from core.folder import Folders
from core.parameters import Parameters
from core.pyhim_logging import Log, write_string_to_file
from imageProcessing.imageProcessing import Image

try:
    del os.environ["MPLBACKEND"]  # unset before importing matplotlib
except:
    print("No environment variable MPLBACKEND found. Continuing anyway.")
# =============================================================================
# FUNCTIONS
# =============================================================================q


def show_image(data_2d, normalization="simple", size=(10, 10)):
    fig = plt.figure()
    fig.set_size_inches(size)

    axes = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
    axes.set_axis_off()

    if normalization == "simple":
        norm = simple_norm(data_2d, "sqrt", percent=99.9)

    fig.add_axes(axes)
    axes.imshow(data_2d, origin="lower", cmap="Greys_r", norm=norm)
    axes.set_title("2D Data")
    return fig, axes


def creates_user_mask(current_param, current_log, file_name, output_filename):
    # loads image

    # displays image
    im_obj = Image(current_param, current_log)
    im_obj.data_2d = np.load(file_name).squeeze()
    im_obj.show_image(show=True, normalization="simple")
    print("Click on the button to add a new ROI")

    # Draw multiple rois
    multiroi_named = MultiRoi(roi_names=["My first ROI", "My second ROI"])

    number_rois = len(multiroi_named.rois)
    print(f"Number of rois drawn: {number_rois}")

    masks = np.zeros((im_obj.data_2d.shape[0], im_obj.data_2d.shape[1], number_rois))

    # Display image with all rois
    fig, _ = show_image(im_obj.data_2d)

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


def processes_user_masks(current_param, current_log, processing_list):
    if current_param.param_dict["segmentedObjects"]["operation"] == "overwrite":
        # session
        session_name = "processes_user_masks"

        # processes folders and files
        data_folder = Folders(current_param.param_dict["rootFolder"])
        data_folder.set_folders()
        current_log.add_simple_text(
            f"\n===================={session_name}====================\n"
        )
        current_log.report(f"folders read: {len(data_folder.list_folders)}")

        position_roi_information = current_param.param_dict["acquisition"][
            "positionROIinformation"
        ]
        number_masked_files = 0

        all_results_table = Table()
        for current_folder in data_folder.list_folders:
            files_folder = glob.glob(current_folder + os.sep + "*.tif")
            data_folder.create_folders(current_folder, current_param)
            current_log.report(f"-------> Processing Folder: {current_folder}")

            # generates lists of files to process
            current_param.find_files_to_process(files_folder)
            file_quantity = len(current_param.files_to_process)
            current_log.report(f"About to process {file_quantity} files\n")

            if file_quantity > 0:
                files_to_process = [
                    file
                    for file in glob.glob(
                        data_folder.output_folders["segmentedObjects"] + os.sep + "*"
                    )
                    if "SNDmask" in file.split("_")
                ]

                if processing_list["cleanAllMasks"]:
                    # [clears all SND masks]
                    number_masked_files = len(files_to_process)
                    for file_name in files_to_process:
                        os.remove(file_name)
                        current_log.report(
                            f"Removing SND mask: {os.path.basename(file_name)}",
                            "info",
                        )

                elif (
                    processing_list["addMask"] == ""
                    and not processing_list["cleanAllMasks"]
                ):
                    # [assigns cells to exsting masks]
                    results_table = assigns_snd_mask2cells(
                        files_to_process, position_roi_information, current_log
                    )
                    all_results_table = vstack([all_results_table, results_table])

                    current_log.report("assigning masks", "info")

                elif (
                    not processing_list["cleanAllMasks"]
                    and len(processing_list["addMask"]) > 0
                ):
                    # [makes new set of masks]
                    for file_name in current_param.files_to_process:
                        # gets filename information
                        roi = os.path.basename(file_name).split("_")[
                            position_roi_information
                        ]

                        # checks that SND channel was projected and aligned
                        registered_filename = (
                            data_folder.output_folders["alignImages"]
                            + os.sep
                            + os.path.basename(file_name).split(".")[0]
                            + "_2d_registered.npy"
                        )

                        if os.path.exists(registered_filename):
                            output_filename = (
                                data_folder.output_folders["segmentedObjects"]
                                + os.sep
                                + os.path.basename(registered_filename).split(".")[0]
                                + "_SNDmask"
                                + "_"
                                + processing_list["addMask"]
                                + ".npy"
                            )
                            creates_user_mask(
                                current_param,
                                current_log,
                                registered_filename,
                                output_filename,
                            )
                            number_masked_files += 1
                            current_log.report(
                                f"Segmented image for SND channel ROI#{roi}: {output_filename}",
                                "info",
                            )
                        else:
                            current_log.report(
                                "Could not find registered image for SND channel: "
                                + str(registered_filename),
                                "error",
                            )

        table_output_filename = (
            data_folder.output_folders["segmentedObjects"]
            + os.sep
            + "snd_assigned_cells.ecsv"
        )
        print(f"\n>Number of cells found: {len(all_results_table)}")
        print(f"\n>Saving output to: {table_output_filename}")
        if os.path.exists(table_output_filename):
            print(f"\nResults appended to {table_output_filename}")
            previous_results_table = Table.read(
                table_output_filename, format="ascii.ecsv"
            )  # ascii.ecsv
            all_results_table = vstack([previous_results_table, all_results_table])
        if len(all_results_table) > 0:
            print("\nTable written")
            all_results_table.write(
                table_output_filename, format="ascii.ecsv", overwrite=True
            )

        return number_masked_files
    else:
        return 0


def assigns_snd_mask2cells(files_to_process, position_roi_information, current_log):
    results_table = Table()

    for file_name in files_to_process:
        print(f"\n-----> Processing: {file_name}")
        roi = os.path.basename(file_name).split("_")[position_roi_information]

        # [checks if DAPI mask exists for the file to process]
        filename_dapi_mask = (
            os.path.dirname(file_name)
            + os.sep
            + "_".join(os.path.basename(file_name).split("_")[0:7])
            + "_ch00_Masks.npy"
        )

        print(f"\nWill search masks in: {filename_dapi_mask}")

        if os.path.exists(filename_dapi_mask) and file_name.split(".")[-1] == "npy":
            # load DAPI mask
            mask_dapi = Image()
            mask_dapi.data_2d = np.load(filename_dapi_mask).squeeze()

            # load SND mask
            print(f"\nWill attemp to match masks from: {file_name}")
            mask_snd = Image()
            mask_snd.data_2d = np.load(file_name, allow_pickle=False).squeeze()

            # matches cells and SND masks
            new_matrix = mask_dapi.data_2d * mask_snd.data_2d
            cells_within_mask = np.unique(new_matrix)

            if len(cells_within_mask) > 0:
                new_table = Table()
                col_mask = Column(
                    [file_name.split("_")[-1].split(".")[0]] * len(cells_within_mask),
                    name="MaskID #",
                    dtype=str,
                )
                col_roi = Column(
                    int(roi) * np.ones(len(cells_within_mask)), name="ROI #", dtype=int
                )
                col_cell_id = Column(cells_within_mask, name="CellID #", dtype=int)
                new_table.add_column(col_roi, index=0)
                new_table.add_column(col_cell_id, index=1)
                new_table.add_column(col_mask, index=2)

                results_table = vstack([results_table, new_table])
                del new_table

            # outputs list of cells within (1) and outside SND mask (0)
            current_log.report(
                f"Matched DAPI and SND masks: \n SND: {file_name}\n DAPI: {filename_dapi_mask}\n",
                "info",
            )

        else:
            print(
                "\nERROR: Could not find expected DAPI mask: "
                + str(filename_dapi_mask)
                + "\n or file does not have ("
                + str(file_name.split(".")[-1])
                + ") the expected extension (.npy)"
            )

    return results_table


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

    print("\n-------------------------------------------------------------------")
    processing_list = {}

    if args.rootFolder:
        root_folder = args.rootFolder
    else:
        root_folder = "."
        # root_folder='/home/marcnol/data/Experiment_4/0_Embryo/alignImages/'
        # fileNameRNA = root_folder+'scan_002_DAPI_001_ROI_converted_decon_ch01_2d_registered.npy'

    if args.addMask:
        processing_list["addMask"] = args.addMask
    else:
        processing_list["addMask"] = ""

    processing_list["cleanAllMasks"] = args.cleanAllMasks

    print(f"parameters> root_folder: {root_folder}")
    now = datetime.now()

    datam = DataManager(root_folder)
    raw_dict = datam.load_user_param()
    current_param = Parameters(raw_dict, root_folder=datam.m_data_path)
    labels = current_param.param_dict["labels"]

    session_name = "processSNDchannel"

    # setup logs
    current_log = Log(root_folder)
    current_log.add_simple_text(
        f"\n^^^^^^^^^^^^^^^^^^^^^^^^^^{session_name}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n"
    )
    current_log.report(f"Process SND channel MD: {current_log.markdown_filename}")
    write_string_to_file(
        current_log.markdown_filename,
        f"""# Process SND channel {now.strftime("%d/%m/%Y %H:%M:%S")}""",
        "w",
    )  # initialises MD file

    for label in labels:
        # sets parameters
        current_param = Parameters(raw_dict, root_folder=datam.m_data_path, label=label)
        print(f"**Analyzing label: {label}**")

        # processes Secondary masks
        if label == "RNA":
            number_masked_files = processes_user_masks(
                current_param, current_log, processing_list
            )
            if number_masked_files > 0:
                print(f"Files processed: {number_masked_files}", "info")
            else:
                print("Nothing processed...", "warning")
        print("\n")
        del current_param

    # exits
    print("\n==================== Normal termination ====================\n")

    print(f"Elapsed time: {datetime.now() - begin_time}")


if __name__ == "__main__":
    main()
