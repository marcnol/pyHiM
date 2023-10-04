#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 09:22:18 2020

@author: marcnol

Sets of functions that do alignment of 2D fiducial images. It also contains
code to apply these alignments to other channels (masks/ barcodes)

For the time being alignment is purely based on optimized sub-pixed accuracy
image cross correlation

"""

# =============================================================================
# IMPORTS
# =============================================================================

import glob
import os
import sys

import cv2
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import SigmaClip
from astropy.table import Table
from dask.distributed import get_client
from numpy import linalg as LA
from photutils import Background2D, MedianBackground
from scipy.ndimage import shift as shift_image
from skimage import exposure, measure
from skimage.exposure import match_histograms
from skimage.metrics import mean_squared_error, normalized_root_mse
from skimage.metrics import structural_similarity as ssim
from skimage.registration import phase_cross_correlation
from skimage.util.shape import view_as_blocks
from tqdm import tqdm, trange

from core.dask_cluster import try_get_client
from core.data_manager import load_json, save_json
from core.folder import Folders
from core.parameters import get_dictionary_value, rt_to_filename
from core.pyhim_logging import print_log, print_session_name, write_string_to_file
from core.saving import plotting_block_alignment_results, save_image_2d_cmd
from imageProcessing.imageProcessing import (
    Image,
    image_adjust,
    reassemble_3d_image,
    scatter_3d_image,
)

np.seterr(divide="ignore", invalid="ignore")

# =============================================================================
# FUNCTIONS
# =============================================================================


def display_equalization_histograms(
    i_histogram, lower_threshold, output_filename, markdown_filename, verbose=False
):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    ax1.plot(i_histogram["Im1"][0][1], i_histogram["Im1"][0][0])
    ax2.plot(i_histogram["Im2"][0][1], i_histogram["Im2"][0][0])
    ax3.plot(i_histogram["Im1"][1][1], i_histogram["Im1"][1][0])
    ax4.plot(i_histogram["Im2"][1][1], i_histogram["Im2"][1][0])
    ax3.set_yscale("log")
    ax4.set_yscale("log")
    ax1.vlines(lower_threshold["Im1"], 0, i_histogram["Im1"][0][0].max(), colors="r")
    ax2.vlines(lower_threshold["Im2"], 0, i_histogram["Im2"][0][0].max(), colors="r")
    plt.savefig(f"{output_filename}_intensityHist.png")
    write_string_to_file(
        markdown_filename,
        f"{os.path.basename(output_filename)}\n ![]({output_filename}_intensityHist.png)\n",
        "a",
    )

    if not verbose:
        plt.close(fig)


def remove_inhomogeneous_background(im, current_param):
    sigma_clip = SigmaClip(
        sigma=current_param.param_dict["alignImages"]["background_sigma"]
    )
    bkg_estimator = MedianBackground()
    bkg = Background2D(
        im,
        (64, 64),
        filter_size=(3, 3),
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
    )

    return im - bkg.background


def save_image_differences(img_1, img_2, img_3, img_4, output_filename):
    """
    Overlays two images as R and B and saves them to output file
    """

    img_1, img_2 = img_1 / img_1.max(), img_2 / img_2.max()
    img_3, img_4 = img_3 / img_3.max(), img_4 / img_4.max()

    img_1, _, _, _, _ = image_adjust(
        img_1, lower_threshold=0.5, higher_threshold=0.9999
    )
    img_2, _, _, _, _ = image_adjust(
        img_2, lower_threshold=0.5, higher_threshold=0.9999
    )
    img_3, _, _, _, _ = image_adjust(
        img_3, lower_threshold=0.5, higher_threshold=0.9999
    )
    img_4, _, _, _, _ = image_adjust(
        img_4, lower_threshold=0.5, higher_threshold=0.9999
    )

    cmap = "seismic"

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches((60, 30))

    ax1.imshow(img_1 - img_2, cmap=cmap)
    ax1.axis("off")
    ax1.set_title("uncorrected")

    ax2.imshow(img_3 - img_4, cmap=cmap)
    ax2.axis("off")
    ax2.set_title("corrected")

    fig.savefig(output_filename)

    plt.close(fig)


def align_images_by_blocks(
    img_1,
    img_2,
    block_size,
    upsample_factor=100,
    min_number_pollsters=4,
    tolerance=0.1,
    use_cv2=False,
    shift_error_tolerance=5,
):
    block_1 = view_as_blocks(img_1, block_size)
    block_2 = view_as_blocks(img_2, block_size)

    if use_cv2:
        warp_matrix = np.eye(2, 3, dtype=np.float32)
        warp_mode = cv2.MOTION_TRANSLATION

    shift_image_norm = np.zeros((block_1.shape[0], block_1.shape[1]))
    shifted_image = np.zeros((block_1.shape[0], block_1.shape[1], 2))
    rms_image = np.zeros((block_1.shape[0], block_1.shape[1]))

    for i in trange(block_1.shape[0]):
        for j in range(block_1.shape[1]):
            if not use_cv2:
                # using Scimage registration functions
                shift, _, _ = phase_cross_correlation(
                    block_1[i, j], block_2[i, j], upsample_factor=upsample_factor
                )
                shift_image_norm[i, j] = LA.norm(shift)
                shifted_image[i, j, 0], shifted_image[i, j, 1] = shift[0], shift[1]
                img_2_aligned = shift_image(img_2, shift)
            else:
                # uses CV2 cause it is 20 times faster than Scimage
                _, warp_matrix = align_cv2(block_1[i, j], block_2[i, j], warp_mode)
                shift_image_norm[i, j] = LA.norm(warp_matrix[:, 2])
                shifted_image[i, j, 0], shifted_image[i, j, 1] = (
                    warp_matrix[:, 2][0],
                    warp_matrix[:, 2][1],
                )
                img_2_aligned = apply_correction(img_2, warp_matrix)

            rms_image[i, j] = np.sum(np.sum(np.abs(img_1 - img_2_aligned), axis=1))

    # [calculates optimal shifts by polling blocks showing the best RMS]

    # threshold = filters.threshold_otsu(rms_image)
    threshold = (1 + tolerance) * np.min(rms_image)
    mask = rms_image < threshold

    contours = measure.find_contours(rms_image, threshold)

    try:
        contour = sorted(contours, key=lambda x: len(x))[-1]
    except IndexError:
        contour = np.array([0, 0])

    # [Averages shifts and errors from regions within the tolerated blocks]
    mean_shifts = [np.mean(shifted_image[mask, 0]), np.mean(shifted_image[mask, 1])]
    std_shifts = [np.std(shifted_image[mask, 0]), np.std(shifted_image[mask, 1])]
    mean_shift_norm = np.mean(shift_image_norm[mask])
    mean_error = np.mean(rms_image[mask])
    relative_shifts = np.abs(shift_image_norm - mean_shift_norm)

    # [calculates global shift, if it is better than the polled shift, or
    # if we do not have enough pollsters to fall back to then it does a global cross correlation!]
    mean_shifts_global, _, _ = phase_cross_correlation(
        img_1, img_2, upsample_factor=100
    )
    img_2_aligned_global = shift_image(img_2, shift)
    mean_error_global = np.sum(np.sum(np.abs(img_1 - img_2_aligned_global), axis=1))

    print_log(
        f"Block alignment error: {mean_error}, global alignment error: {mean_error_global}"
    )

    if (
        np.sum(mask) < min_number_pollsters
        or mean_error_global < mean_error
        or np.max(std_shifts) > shift_error_tolerance
    ):
        mean_shifts = mean_shifts_global
        mean_error = mean_error_global
        print_log("Falling back to global registration")

    print_log(
        f"*** Global XY shifts: {mean_shifts_global[0]:.2f} px | {mean_shifts_global[1]:.2f} px"
    )
    print_log(
        f"*** Mean polled XY shifts: {mean_shifts[0]:.2f}({std_shifts[0]:.2f}) px | {mean_shifts[1]:.2f}({std_shifts[1]:.2f}) px"
    )

    return np.array(mean_shifts), mean_error, relative_shifts, rms_image, contour


def save_2_images_rgb(img_1, img_2, output_filename):
    """
    Overlays two images as R and B and saves them to output file
    """

    sz = img_1.shape
    img_1, img_2 = img_1 / img_1.max(), img_2 / img_2.max()

    img_1, _, _, _, _ = image_adjust(
        img_1, lower_threshold=0.5, higher_threshold=0.9999
    )
    img_2, _, _, _, _ = image_adjust(
        img_2, lower_threshold=0.5, higher_threshold=0.9999
    )

    fig, ax1 = plt.subplots()
    fig.set_size_inches((30, 30))

    null_image = np.zeros(sz)

    rgb = np.dstack([img_1, img_2, null_image])
    ax1.imshow(rgb)
    ax1.axis("off")

    fig.savefig(output_filename)

    plt.close(fig)


def align_2_files(file_name, img_reference, current_param, data_folder, verbose):
    """
    Uses preloaded ImReference Object and aligns it against filename

    Parameters
    ----------
    file_name : npy 2D array
        file of image to be aligned
    img_reference : Image Class
        Object type <Image> with image reference
    current_param : Parameters Class
        Running parameters
    data_folder : Folders Class
        DESCRIPTION.
    verbose : boolean
        True for display images

    Returns are returned as arguments!
    -------
    shift : float list, 2 dimensions
        offset in Y and X
    table_entry : Table Class
        results zipped in Table Class form

    """
    filename_1 = img_reference.file_name
    filename_2 = file_name

    output_filename = (
        data_folder.output_folders["alignImages"]
        + os.sep
        + os.path.basename(filename_2).split(".")[0]
    )

    # loads image
    img_2 = Image(current_param)
    img_2.load_image_2d(filename_2, data_folder.output_folders["zProject"])

    # Normalises images
    image1_uncorrected = img_reference.data_2d / img_reference.data_2d.max()
    image2_uncorrected = img_2.data_2d / img_2.data_2d.max()

    # removes inhomogeneous background
    image1_uncorrected = remove_inhomogeneous_background(
        image1_uncorrected, current_param
    )
    image2_uncorrected = remove_inhomogeneous_background(
        image2_uncorrected, current_param
    )

    lower_threshold = get_dictionary_value(
        current_param.param_dict["alignImages"], "lower_threshold", default=0.999
    )
    higher_threshold = get_dictionary_value(
        current_param.param_dict["alignImages"], "higher_threshold", default=0.9999999
    )
    align_by_block = get_dictionary_value(
        current_param.param_dict["alignImages"], "alignByBlock", default=False
    )
    tolerance = get_dictionary_value(
        current_param.param_dict["alignImages"], "tolerance", default=0.1
    )
    dict_block_size = get_dictionary_value(
        current_param.param_dict["alignImages"], "blockSize", default=256
    )

    if not align_by_block:
        # [calculates unique translation for the entire image using cross-correlation]
        (
            shift,
            error,
            diffphase,
            lower_threshold,
            i_histogram,
            _,
            _,
            _,
        ) = align_2_images_cross_correlation(
            image1_uncorrected,
            image2_uncorrected,
            lower_threshold=lower_threshold,
            higher_threshold=higher_threshold,
        )

        # displays intensity histograms
        display_equalization_histograms(
            i_histogram,
            lower_threshold,
            output_filename,
            current_param.param_dict["fileNameMD"],
            verbose,
        )

    else:
        # [calculates block translations by cross-correlation and gets overall shift by polling]

        # normalizes images
        image1_uncorrected, image2_uncorrected = (
            np.float32(image1_uncorrected),
            np.float32(image2_uncorrected),
        )

        # matches histograms
        image2_uncorrected = np.float32(
            match_histograms(image2_uncorrected, image1_uncorrected)
        )

        # calculates block shifts and polls for most favourable shift
        upsample_factor = 100
        block_size = (dict_block_size, dict_block_size)

        (
            shift,
            error,
            relative_shifts,
            rms_image,
            contour,
        ) = align_images_by_blocks(
            image1_uncorrected,
            image2_uncorrected,
            block_size,
            upsample_factor=upsample_factor,
            min_number_pollsters=4,
            tolerance=tolerance,
        )
        diffphase = 0

        plotting_block_alignment_results(
            relative_shifts,
            rms_image,
            contour,
            file_name=f"{output_filename}_block_alignments.png",
        )

        write_string_to_file(
            current_param.param_dict["fileNameMD"],
            f"{os.path.basename(output_filename)}\n ![]({output_filename}_block_alignments.png)\n",
            "a",
        )

        # saves mask of valid regions with a correction within the tolerance
        save_image_2d_cmd(rms_image, f"{output_filename}_rmsBlockMap")
        save_image_2d_cmd(relative_shifts, f"{output_filename}_errorAlignmentBlockMap")

    image2_corrected_raw = shift_image(image2_uncorrected, shift)

    image2_corrected_raw[image2_corrected_raw < 0] = 0

    error = np.sum(np.sum(np.abs(image1_uncorrected - image2_corrected_raw), axis=1))

    print_log(f"$ Detected subpixel offset (y, x): {shift} px")

    # [displays and saves results]

    # thresholds corrected images for better display and saves
    image1_uncorrected[image1_uncorrected < 0] = 0
    image2_uncorrected[image2_uncorrected < 0] = 0

    save_2_images_rgb(
        image1_uncorrected,
        image2_corrected_raw,
        f"{output_filename}_overlay_corrected.png",
    )

    save_image_differences(
        image1_uncorrected,
        image2_uncorrected,
        image1_uncorrected,
        image2_corrected_raw,
        f"{output_filename}_referenceDifference.png",
    )

    # reports image in MD file
    write_string_to_file(
        current_param.param_dict["fileNameMD"],
        f"{os.path.basename(output_filename)}\n ![]({output_filename}_overlay_corrected.png)\n ![]({output_filename}_referenceDifference.png)\n",
        "a",
    )

    # outputs results to logfile
    alignment_output = data_folder.output_files["alignImages"]
    text_to_output = f"{os.path.basename(filename_2)}\t{os.path.basename(filename_1)}\t{shift[0]:.2f}\t{shift[1]:.2f}\t{error:.2f}\t{diffphase:.2f}"
    write_string_to_file(alignment_output, text_to_output, "a")

    # creates Table entry to return
    table_entry = [
        os.path.basename(filename_2),
        os.path.basename(filename_1),
        shift[0],
        shift[1],
        error,
        diffphase,
    ]

    # saves registered fiducial image
    save_image_2d_cmd(image2_corrected_raw, f"{output_filename}_2d_registered")

    del img_2
    return shift, table_entry


def align_images_in_current_folder(
    current_folder, current_param, data_folder, current_session, file_name=None
):
    alignment_results_table = Table(
        names=(
            "aligned file",
            "reference file",
            "shift_x",
            "shift_y",
            "error",
            "diffphase",
        ),
        dtype=("S2", "S2", "f4", "f4", "f4", "f4"),
    )

    # initializes variables
    files_folder = glob.glob(current_folder + os.sep + "*.tif")
    data_folder.create_folders(current_folder, current_param)

    # generates lists of files to process for current_folder
    current_param.find_files_to_process(files_folder)
    print_log(f"> Processing Folder: {current_folder}")
    print_log(f"> About to process {len(current_param.files_to_process)} files\n")
    write_string_to_file(
        data_folder.output_files["alignImages"],
        "File1 \t File_reference \t shift_y \t shift_x \t error \t diffphase",
        "w",
    )

    # Finds and loads Reference fiducial information
    reference_barcode = current_param.param_dict["alignImages"]["referenceFiducial"]
    print_log(f"$ Reference fiducial {reference_barcode}")

    # retrieves the list of fiducial image files to be aligned
    filenames_with_ref_barcode, roi_list = rt_to_filename(
        current_param, reference_barcode
    )

    if len(filenames_with_ref_barcode) > 0:
        verbose = False
        # contains dictionary of shifts for each folder
        dict_shifts = {}
        session_name = "register_global"
        # loops over fiducials images one ROI at a time
        for filename_reference in filenames_with_ref_barcode:
            # loads reference fiducial image for this ROI
            roi = roi_list[filename_reference]
            img_reference = Image(current_param)
            img_reference.load_image_2d(
                filename_reference, data_folder.output_folders["zProject"]
            )
            print_log(f"> Loading reference Image {filename_reference}")

            # saves reference 2D image of fiducial
            if not os.path.exists(
                img_reference.get_image_filename(
                    data_folder.output_folders["alignImages"], tag="_2d_registered"
                )
            ):
                img_reference.save_image_2d(
                    data_folder.output_folders["alignImages"],
                    tag="_2d_registered",
                )

            dict_shift_roi = {}

            filenames_to_process_list = [
                x
                for x in current_param.files_to_process
                if (x != filename_reference)
                and current_param.decode_file_parts(os.path.basename(x))["roi"] == roi
            ]
            print_log(f"Found {len(filenames_to_process_list)} files in ROI: {roi}")
            print_log(
                f'[roi:cycle] {"|".join([str(current_param.decode_file_parts(os.path.basename(x))["roi"]) + ":" + str(current_param.decode_file_parts(os.path.basename(x))["cycle"]) for x in filenames_to_process_list])}'
            )

            if current_param.param_dict["parallel"]:
                # running in parallel mode
                client = get_client()
                futures = []
                labels = []

                for filename_to_process in filenames_to_process_list:
                    # excludes the reference fiducial and processes files in the same ROI
                    labels.append(os.path.basename(filename_to_process).split("_")[2])
                    futures.append(
                        client.submit(
                            align_2_files,
                            filename_to_process,
                            img_reference,
                            current_param,
                            data_folder,
                            verbose,
                        )
                    )

                print_log(f"$ Waiting for {len(futures)} results to arrive")

                results = client.gather(futures)

                print_log(f"$ Retrieving {len(results)} results from cluster")

                for result, label in zip(results, labels):
                    shift, table_entry = result
                    dict_shift_roi[label] = shift.tolist()
                    alignment_results_table.add_row(table_entry)
                    # TODO: filename_to_process var doesn't exist in this scope, why and what's happend ?
                    current_session.add(filename_to_process, session_name)

            else:
                # running in sequential mode
                n_files = len(current_param.files_to_process)

                for i_file, filename_to_process in enumerate(
                    current_param.files_to_process
                ):
                    # excludes the reference fiducial and processes files in the same ROI
                    label = os.path.basename(filename_to_process).split("_")[2]
                    roi = current_param.decode_file_parts(
                        os.path.basename(filename_to_process)
                    )["roi"]
                    print_log(f"\n$ About to process file {i_file} \\ {n_files}")

                    if filename_to_process in filename_reference:
                        print_log(
                            f"\n$ Skipping reference file: {os.path.basename(filename_to_process)} "
                        )
                    elif file_name is None or (
                        file_name is not None
                        and os.path.basename(file_name)
                        == os.path.basename(filename_to_process)
                    ):
                        # aligns files and saves results to database in dict format and to a Table
                        shift, table_entry = align_2_files(
                            filename_to_process,
                            img_reference,
                            current_param,
                            data_folder,
                            verbose,
                        )
                        dict_shift_roi[label] = shift.tolist()
                        alignment_results_table.add_row(table_entry)
                        current_session.add(filename_to_process, session_name)
            # accumulates shifst for this ROI into global dictionary
            dict_shifts[f"ROI:{roi}"] = dict_shift_roi
            del img_reference

        # saves dicShifts dictionary with shift results
        dictionary_filename = (
            os.path.splitext(data_folder.output_files["dictShifts"])[0] + ".json"
        )
        save_json(dict_shifts, dictionary_filename)
        print_log(f"$ Saved alignment dictionary to {dictionary_filename}")

    else:
        print_log(f"# Reference Barcode file does not exist: {reference_barcode}")
        raise ValueError

    return alignment_results_table


def align_images(current_param, current_session, file_name=None):
    """
    From a given parameters class it aligns all the fiducial images

    Parameters
    ----------
    current_param : Parameters class
        running parameters
    current_session : Session Class
        logs session information.

    Returns
    -------
    None.

    """
    session_name = "register_global"

    # processes folders and adds information to log files
    data_folder = Folders(current_param.param_dict["rootFolder"])
    data_folder.set_folders()
    print_session_name(session_name)
    print_log(f"folders read: {len(data_folder.list_folders)}")
    write_string_to_file(
        current_param.param_dict["fileNameMD"],
        f"""## {session_name}: {current_param.param_dict["acquisition"]["label"]}\n""",
        "a",
    )

    # loops over folders
    for current_folder in data_folder.list_folders:
        alignment_results_table = align_images_in_current_folder(
            current_folder, current_param, data_folder, current_session, file_name
        )

    # saves Table with all shifts
    alignment_results_table.write(
        data_folder.output_files["alignImages"].split(".")[0] + ".table",
        format="ascii.ecsv",
        overwrite=True,
    )

    del data_folder


def apply_registrations_to_filename(
    filename_to_process, current_param, data_folder, current_session, dict_shifts
):
    """
    Applies registration of filename_to_process

    Parameters
    ----------
    filename_to_process : string
        file to apply registration to
    current_param : Parameters class
    data_folder : data_folder class
    current_session : Session class
    dict_shifts : Dictionnary
        contains the shifts to be applied to all rois

    Returns
    -------
    None.

    """
    # gets shift from dictionary
    # ROI = os.path.basename(filename_to_process).split("_")[position_roi_information]
    roi = current_param.decode_file_parts(os.path.basename(filename_to_process))["roi"]

    label = os.path.basename(filename_to_process).split("_")[2]  # to FIX

    try:
        shift_array = dict_shifts[f"ROI:{roi}"][label]
    except KeyError:
        shift_array = None
        print_log(
            f"$ Could not find dictionary with alignment parameters for this ROI: ROI:{roi}, label: {label}"
        )

    if shift_array is not None:
        shift = np.asarray(shift_array)
        # loads 2D image and applies registration
        im_obj = Image(current_param)
        im_obj.load_image_2d(
            filename_to_process, data_folder.output_folders["zProject"]
        )
        im_obj.data_2d = shift_image(im_obj.data_2d, shift)
        print_log(f"$ Image registered using ROI:{roi}, label:{label}, shift={shift}")

        # saves registered 2D image
        im_obj.save_image_2d(
            data_folder.output_folders["alignImages"],
            tag="_2d_registered",
        )

        # session
        session_name = "register_global"

        # logs output
        current_session.add(filename_to_process, session_name)
    elif label == current_param.param_dict["alignImages"]["referenceFiducial"]:
        im_obj = Image(current_param)
        im_obj.load_image_2d(
            filename_to_process, data_folder.output_folders["zProject"]
        )
        im_obj.save_image_2d(
            data_folder.output_folders["alignImages"],
            tag="_2d_registered",
        )
        print_log(f"$ Saving image for referenceRT ROI:{roi}, label:{label}")

    else:
        print_log(
            f"# No shift found in dictionary for ROI:{roi}, label:{label}",
            status="WARN",
        )


def apply_registrations_to_current_folder(
    current_folder, current_param, data_folder, current_session, file_name=None
):
    """
    applies registrations to all files in current_folder

    Parameters
    ----------
    current_folder : TYPE
        DESCRIPTION.
    current_param : Parameters class
    data_folder : data_folder class
    current_session : Session class
    file_name : string, optional
        File to process. The default is None.

    Returns
    -------
    None.

    """

    # current_folder=data_folder.list_folders[0] # only one folder processed so far...
    files_folder = glob.glob(current_folder + os.sep + "*.tif")
    data_folder.create_folders(current_folder, current_param)
    print_log(f"> Processing Folder: {current_folder}")

    # loads dicShifts with shifts for all rois and all labels
    dict_filename = (
        os.path.splitext(data_folder.output_files["dictShifts"])[0] + ".json"
    )

    # dict_filename = data_folder.output_files["dictShifts"] + ".json"
    dict_shifts = load_json(dict_filename)
    if len(dict_shifts) == 0:
        print_log(f"# File with dictionary not found!: {dict_filename}")
    else:
        print_log(f"$ Dictionary File loaded: {dict_filename}")

    # generates lists of files to process
    current_param.find_files_to_process(files_folder)
    n_files = len(current_param.files_to_process)
    print_log(f"\n$ About to process {n_files} files\n")

    if len(current_param.files_to_process) > 0:
        # loops over files in file list
        for i, filename_to_process in enumerate(current_param.files_to_process):
            if file_name is None or (
                file_name is not None
                and os.path.basename(file_name) == os.path.basename(filename_to_process)
            ):
                print_log(f"\n$ About to process file {i} \\ {n_files}")
                apply_registrations_to_filename(
                    filename_to_process,
                    current_param,
                    data_folder,
                    current_session,
                    dict_shifts,
                )


def apply_registrations(current_param, current_session, file_name=None):
    """This function will
    - load masks, RNA and barcode 2D projected images,
    - apply registrations
    - save registered images as npy arrays
    """

    session_name = "register_global"

    # verbose=False

    # processes folders and files
    data_folder = Folders(current_param.param_dict["rootFolder"])
    data_folder.set_folders()
    print_session_name(session_name)
    print_log(f"$ folders read: {len(data_folder.list_folders)}")

    for current_folder in data_folder.list_folders:
        apply_registrations_to_current_folder(
            current_folder, current_param, data_folder, current_session, file_name
        )

    del data_folder


# =============================================================================
# IMAGE ALIGNMENT
# =============================================================================


def apply_xy_shift_3d_images(image, shift, parallel_execution=True):
    """
    Applies XY shift to a 3D stack

    Parameters
    ----------
    images : 3D numpy array
        image to process.

    Returns
    -------
    shifted 3D image.

    """
    client = try_get_client() if parallel_execution else None
    number_planes = image.shape[0]

    if client is None:
        print_log(f"> Shifting {number_planes} planes with 1 thread...")
        shift_3d = np.zeros((3))
        shift_3d[0], shift_3d[1], shift_3d[2] = 0, shift[0], shift[1]
        output = shift_image(image, shift_3d)
    else:
        print_log(
            f"> Shifting {number_planes} planes using {len(client.scheduler_info()['workers'])} workers..."
        )

        image_list_scattered = scatter_3d_image(image)

        futures = [
            client.submit(shift_image, img, shift) for img in image_list_scattered
        ]

        output = reassemble_3d_image(client, futures, image.shape)

        del futures
        del image_list_scattered

    print_log("$ Done shifting 3D image.")

    return output


def image_block_alignment_3d(images, block_size_xy=256, upsample_factor=100):
    # sanity checks
    if len(images) < 2:
        sys.exit(f"# Error, number of images must be 2, not {len(images)}")

    # - break in blocks
    num_planes = images[0].shape[0]
    block_size = (num_planes, block_size_xy, block_size_xy)

    print_log("$ Breaking images into blocks")
    blocks = [view_as_blocks(x, block_shape=block_size).squeeze() for x in images]

    block_ref = blocks[0]
    block_target = blocks[1]

    # - loop thru blocks and calculates block shift in xyz:
    shift_matrices = [np.zeros(block_ref.shape[:2]) for _ in range(3)]

    for i in trange(block_ref.shape[0]):
        for j in range(block_ref.shape[1]):
            # - cross correlate in 3D to find 3D shift
            shifts_xyz, _, _ = phase_cross_correlation(
                block_ref[i, j], block_target[i, j], upsample_factor=upsample_factor
            )
            for matrix, _shift in zip(shift_matrices, shifts_xyz):
                matrix[i, j] = _shift

    return shift_matrices, block_ref, block_target


def combine_blocks_image_by_reprojection(
    block_ref, block_target, shift_matrices=None, axis1=0
):
    """
    This routine will overlap block_ref and block_target images block by block.
    block_ref will be used as a template.
    - block_target will be first translated in ZXY using the corresponding values in shift_matrices
    to realign each block
    - then an rgb image will be created with block_ref in the red channel, and the reinterpolated
    block_target block in the green channel.
    - the Blue channel is used for the grid to improve visualization of blocks.


    Parameters
    ----------
    block_ref : npy array
        return of view_as_blocks()
    block_target : npy array
        return of view_as_blocks()
    shift_matrices : list of npy arrays
        index 0 contains Z, index 1 X and index 2 Y
    axis1 : int
        axis used for the reprojection: The default is 0.
        - 0 means an XY projection
        - 1 an ZX projection
        - 2 an ZY projection

    Returns
    -------
    output : NPY array of size im_size x im_size x 3
        rgb image.
    ssim_as_blocks = NPY array of size number_blocks x number_blocks
        Structural similarity index between ref and target blocks
    """
    number_blocks = block_ref.shape[0]
    block_sizes = list(block_ref.shape[2:])
    block_sizes.pop(axis1)
    img_sizes = [x * number_blocks for x in block_sizes]

    # gets ranges for slicing
    slice_coordinates = [
        [range(x * block_size, (x + 1) * block_size) for x in range(number_blocks)]
        for block_size in block_sizes
    ]

    # creates output images
    output = np.zeros((img_sizes[0], img_sizes[1], 3))
    ssim_as_blocks = np.zeros((number_blocks, number_blocks))
    mse_as_blocks = np.zeros((number_blocks, number_blocks))
    nrmse_as_blocks = np.zeros((number_blocks, number_blocks))

    # blank image for blue channel to show borders between blocks
    blue = np.zeros(block_sizes)
    blue[0, :], blue[:, 0], blue[:, -1], blue[-1, :] = [0.5] * 4

    # reassembles image
    # takes one plane block
    for i, i_slice in enumerate(tqdm(slice_coordinates[0])):
        for j, j_slice in enumerate(slice_coordinates[1]):
            imgs = [block_ref[i, j]]
            if shift_matrices is not None:
                shift_3d = np.array(
                    [x[i, j] for x in shift_matrices]
                )  # gets 3D shift from block decomposition
                imgs.append(
                    shift_image(block_target[i, j], shift_3d)
                )  # realigns and appends to image list
            else:
                imgs.append(
                    block_target[i, j]
                )  # appends original target with no re-alignment

            imgs = [np.sum(x, axis=axis1) for x in imgs]  # projects along axis1
            imgs = [
                exposure.rescale_intensity(x, out_range=(0, 1)) for x in imgs
            ]  # rescales intensity values
            imgs = [
                image_adjust(x, lower_threshold=0.5, higher_threshold=0.9999)[0]
                for x in imgs
            ]  # adjusts pixel intensities

            nrmse_as_blocks[i, j] = normalized_root_mse(
                imgs[0], imgs[1], normalization="euclidean"
            )
            mse_as_blocks[i, j] = mean_squared_error(imgs[0], imgs[1])
            ssim_as_blocks[i, j] = ssim(
                imgs[0], imgs[1], data_range=imgs[1].max() - imgs[1].min()
            )

            imgs.append(blue)  # appends last channel with grid

            rgb = np.dstack(imgs)  # makes block rgb image

            output[
                i_slice[0] : i_slice[-1] + 1, j_slice[0] : j_slice[-1] + 1, :
            ] = rgb  # inserts block into final rgb stack

    return output, ssim_as_blocks, mse_as_blocks, nrmse_as_blocks


def align_2_images_cross_correlation(
    image1_uncorrected,
    image2_uncorrected,
    lower_threshold=0.999,
    higher_threshold=0.9999999,
    upsample_factor=100,
):
    """
    Aligns 2 images by contrast adjust and cross correlation
    Parameters
    ----------
    img_reference : TYPE
        DESCRIPTION.
    img_2 : TYPE
        DESCRIPTION.

    Returns
    -------
    shift : TYPE
        DESCRIPTION.
    error : TYPE
        DESCRIPTION.
    diffphase : TYPE
        DESCRIPTION.
    lower_threshold : TYPE
        DESCRIPTION.
    i_histogram : TYPE
        DESCRIPTION.
    image2_corrected : TYPE
        DESCRIPTION.
    image1_adjusted : TYPE
        DESCRIPTION.
    image2_adjusted : TYPE
        DESCRIPTION.
    image2_corrected_raw : TYPE
        DESCRIPTION.

    """

    (
        image1_adjusted,
        hist1_before,
        hist1_after,
        lower_cutoff1,
        _,
    ) = image_adjust(
        image1_uncorrected,
        lower_threshold=lower_threshold,
        higher_threshold=higher_threshold,
    )
    (
        image2_adjusted,
        hist2_before,
        hist2_after,
        lower_cutoff2,
        _,
    ) = image_adjust(
        image2_uncorrected,
        lower_threshold=lower_threshold,
        higher_threshold=higher_threshold,
    )

    # zips histograms
    lower_threshold = {"Im1": lower_cutoff1, "Im2": lower_cutoff2}
    i_histogram = {
        "Im1": (hist1_before, hist1_after),
        "Im2": (hist2_before, hist2_after),
    }

    # calculates shift
    shift, error, diffphase = phase_cross_correlation(
        image1_adjusted, image2_adjusted, upsample_factor=upsample_factor
    )

    # corrects image
    # The shift corresponds to the pixel offset relative to the reference image
    image2_corrected = shift_image(image2_adjusted, shift)
    image2_corrected = exposure.rescale_intensity(image2_corrected, out_range=(0, 1))

    return (
        shift,
        error,
        diffphase,
        lower_threshold,
        i_histogram,
        image2_corrected,
        image1_adjusted,
        image2_adjusted,
    )


def align_cv2(im1, im2, warp_mode):
    # Define 2x3 or 3x3 matrices and initialize the matrix to identity
    if warp_mode == cv2.MOTION_HOMOGRAPHY:
        warp_matrix = np.eye(3, 3, dtype=np.float32)
    else:
        warp_matrix = np.eye(2, 3, dtype=np.float32)

    # Specify the number of iterations.
    number_of_iterations = 1000  # 5000

    # Specify the threshold of the increment
    # in the correlation coefficient between two iterations
    termination_eps = 1e-10

    # Define termination criteria
    criteria = (
        cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT,
        number_of_iterations,
        termination_eps,
    )

    # Run the ECC algorithm. The results are stored in warp_matrix.
    try:
        cc, warp_matrix = cv2.findTransformECC(
            im1, im2, warp_matrix, warp_mode, criteria, inputMask=None, gaussFiltSize=1
        )
    except TypeError:
        cc, warp_matrix = cv2.findTransformECC(
            im1, im2, warp_matrix, warp_mode, criteria
        )
    except cv2.error:
        cc = 0
        # print_log('Warning: find transform failed. Set warp as identity')

    return cc, warp_matrix


def apply_correction(im2, warp_matrix):
    sz = im2.shape

    return cv2.warpAffine(
        im2, warp_matrix, (sz[1], sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP
    )
