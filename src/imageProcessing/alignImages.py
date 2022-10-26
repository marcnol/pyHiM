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

# to remove in a future version
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import SigmaClip
from astropy.table import Table
from dask.distributed import get_client
from photutils import Background2D, MedianBackground
from scipy.ndimage import shift as shift_image
from skimage.exposure import match_histograms
from skimage.registration._phase_cross_correlation import _upsampled_dft

from fileProcessing.fileManagement import (
    Folders,
    get_dictionary_value,
    load_json,
    print_log,
    rt_to_filename,
    save_json,
    write_string_to_file,
)
from imageProcessing.imageProcessing import (
    Image,
    align_2_images_cross_correlation,
    align_images_by_blocks,
    plotting_block_alignment_results,
    save_2_images_rgb,
    save_image_2d_cmd,
    save_image_differences,
)

warnings.filterwarnings("ignore")
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
    ax1.vlines(
        lower_threshold["Im1"], 0, i_histogram["Im1"][0][0].max(), colors="r"
    )
    ax2.vlines(
        lower_threshold["Im2"], 0, i_histogram["Im2"][0][0].max(), colors="r"
    )
    plt.savefig(output_filename + "_intensityHist.png")
    write_string_to_file(
        markdown_filename, 
        f"{os.path.basename(output_filename)}\n ![]({output_filename}_intensityHist.png)\n", 
        "a",
    )

    if not verbose:
        plt.close(fig)


def show_cc_image(
    image1_uncorrected, image2_uncorrected, output_filename, shift, verbose=False
):
    image_product = (
        np.fft.fft2(image1_uncorrected) * np.fft.fft2(image2_uncorrected).conj()
    )
    cc_image = _upsampled_dft(image_product, 150, 100, (shift * 100) + 75).conj()
    if verbose:
        plt.figure(figsize=(60, 30))
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
        rgb_falsecolor_image = np.dstack(
            [image1_uncorrected, image2_uncorrected, np.zeros([2048, 2048])]
        )
        ax1.imshow(rgb_falsecolor_image, origin="lower", interpolation="nearest")
        ax1.set_axis_off()
        ax1.set_title("Super-imposed images")

        ax2.imshow(cc_image.real)
        ax2.set_axis_off()
        ax2.set_title("Supersampled XC sub-area")
    else:
        plt.figure(figsize=(30, 30))
        plt.imsave(output_filename + "_CC.png", cc_image)
        plt.close()


def save_image_adjusted(file_name, markdown_filename, image1):
    plt.figure(figsize=(30, 30))
    plt.imsave(file_name + "_adjusted.png", image1, cmap="hot")
    write_string_to_file(
        markdown_filename,
        f"{os.path.basename(file_name)}\n ![]({file_name}_adjusted.png)\n",
        "a",
    )
    plt.close()


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

    im1_bkg_substracted = im - bkg.background

    return im1_bkg_substracted


def align_2_files(
    file_name, img_reference, current_param, data_folder, verbose
):
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

        (shift, error, relative_shifts, rms_image, contour,) = align_images_by_blocks(
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
            file_name=output_filename + "_block_alignments.png",
        )

        write_string_to_file(
            current_param.param_dict["fileNameMD"],
            f"{os.path.basename(output_filename)}\n ![]({output_filename}_block_alignments.png)\n",
            "a",
        )

        # saves mask of valid regions with a correction within the tolerance
        save_image_2d_cmd(rms_image, output_filename + "_rmsBlockMap")
        save_image_2d_cmd(relative_shifts, output_filename + "_errorAlignmentBlockMap")

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
        output_filename + "_overlay_corrected.png",
    )

    save_image_differences(
        image1_uncorrected,
        image2_uncorrected,
        image1_uncorrected,
        image2_corrected_raw,
        output_filename + "_referenceDifference.png",
    )

    # reports image in MD file
    write_string_to_file(
        current_param.param_dict["fileNameMD"],
        f"{os.path.basename(output_filename)}\n ![]({output_filename}_overlay_corrected.png)\n ![]({output_filename}_referenceDifference.png)\n",
        "a",
    )

    # outputs results to logfile
    alignment_output = data_folder.output_files["alignImages"]
    list_to_output = "{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(
        os.path.basename(filename_2),
        os.path.basename(filename_1),
        shift[0],
        shift[1],
        error,
        diffphase,
    )
    write_string_to_file(alignment_output, list_to_output, "a")

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
    save_image_2d_cmd(image2_corrected_raw, output_filename + "_2d_registered")

    del img_2
    return shift, table_entry


def align_images_in_current_folder(
    current_folder, current_param, data_folder, current_session, file_name=None
):
    # session
    session_name = "alignImages"
    verbose = False

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
    dict_shifts = (
        {}
    )  # defaultdict(dict) # contains dictionary of shifts for each folder

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
                    data_folder.output_folders["alignImages"], tag="_2d_registered",
                )

            dict_shift_roi = {}

            filenames_to_process_list = [
                x
                for x in current_param.files_to_process
                if (x not in filename_reference)
                and current_param.decode_file_parts(os.path.basename(x))["roi"] == roi
            ]
            print_log(
                "Found {} files in ROI: {}".format(len(filenames_to_process_list), roi)
            )
            print_log(
                "[roi:cycle] {}".format(
                    "|".join(
                        [
                            str(
                                current_param.decode_file_parts(os.path.basename(x))[
                                    "roi"
                                ]
                            )
                            + ":"
                            + str(
                                current_param.decode_file_parts(os.path.basename(x))[
                                    "cycle"
                                ]
                            )
                            for x in filenames_to_process_list
                        ]
                    )
                )
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

                print_log("$ Waiting for {} results to arrive".format(len(futures)))

                results = client.gather(futures)

                print_log("$ Retrieving {} results from cluster".format(len(results)))

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
                    print_log("\n$ About to process file {} \\ {}".format(i_file, n_files))

                    if (filename_to_process not in filename_reference) and roi == roi:
                        if file_name is None or (
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
                    elif filename_to_process in filename_reference:
                        print_log(
                            "\n$ Skipping reference file: {} ".format(
                                os.path.basename(filename_to_process)
                            )
                        )
            # accumulates shifst for this ROI into global dictionary
            dict_shifts["ROI:" + roi] = dict_shift_roi
            del img_reference

        # saves dicShifts dictionary with shift results
        dictionary_filename = (
            os.path.splitext(data_folder.output_files["dictShifts"])[0] + ".json"
        )
        save_json(dictionary_filename, dict_shifts)
        print_log("$ Saved alignment dictionary to {}".format(dictionary_filename))

    else:
        print_log(
            "# Reference Barcode file does not exist: {}".format(reference_barcode)
        )
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
    session_name = "registersImages"

    # processes folders and adds information to log files
    data_folder = Folders(current_param.param_dict["rootFolder"])
    data_folder.set_folders()
    print_log("\n===================={}====================\n".format(session_name))
    print_log("folders read: {}".format(len(data_folder.list_folders)))
    write_string_to_file(
        current_param.param_dict["fileNameMD"],
        "## {}: {}\n".format(
            session_name, current_param.param_dict["acquisition"]["label"]
        ),
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
    # session
    session_name = "registersImages"

    # gets shift from dictionary
    # ROI = os.path.basename(filename_to_process).split("_")[position_roi_information]
    roi = current_param.decode_file_parts(os.path.basename(filename_to_process))["roi"]

    label = os.path.basename(filename_to_process).split("_")[2]  # to FIX

    try:
        shift_array = dict_shifts["ROI:" + roi][label]
    except KeyError:
        shift_array = None
        print_log(
            "$ Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(
                "ROI:" + roi, label
            )
        )

    if shift_array is not None:

        shift = np.asarray(shift_array)
        # loads 2D image and applies registration
        im_obj = Image(current_param)
        im_obj.load_image_2d(filename_to_process, data_folder.output_folders["zProject"])
        im_obj.data_2d = shift_image(im_obj.data_2d, shift)
        print_log(
            "$ Image registered using ROI:{}, label:{}, shift={}".format(
                roi, label, shift
            )
        )

        # saves registered 2D image
        im_obj.save_image_2d(
            data_folder.output_folders["alignImages"], tag="_2d_registered",
        )

        # logs output
        current_session.add(filename_to_process, session_name)
    elif (
        shift_array is None
        and label == current_param.param_dict["alignImages"]["referenceFiducial"]
    ):
        im_obj = Image(current_param)
        im_obj.load_image_2d(filename_to_process, data_folder.output_folders["zProject"])
        im_obj.save_image_2d(
            data_folder.output_folders["alignImages"], tag="_2d_registered",
        )
        print_log("$ Saving image for referenceRT ROI:{}, label:{}".format(roi, label))

    else:
        print_log(
            "# No shift found in dictionary for ROI:{}, label:{}".format(roi, label),
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
    print_log("> Processing Folder: {}".format(current_folder))

    # loads dicShifts with shifts for all rois and all labels
    dict_filename = (
        os.path.splitext(data_folder.output_files["dictShifts"])[0] + ".json"
    )

    # dict_filename = data_folder.output_files["dictShifts"] + ".json"
    dict_shifts = load_json(dict_filename)
    if len(dict_shifts) == 0:
        print_log("# File with dictionary not found!: {}".format(dict_filename))
    else:
        print_log("$ Dictionary File loaded: {}".format(dict_filename))

    # generates lists of files to process
    current_param.find_files_to_process(files_folder)
    n_files = len(current_param.files_to_process)
    print_log("\n$ About to process {} files\n".format(n_files))

    if len(current_param.files_to_process) > 0:
        # loops over files in file list
        for i, filename_to_process in enumerate(current_param.files_to_process):
            if file_name is None or (
                file_name is not None
                and os.path.basename(file_name) == os.path.basename(filename_to_process)
            ):
                print_log("\n$ About to process file {} \\ {}".format(i, n_files))
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

    session_name = "registersImages"

    # verbose=False

    # processes folders and files
    data_folder = Folders(current_param.param_dict["rootFolder"])
    data_folder.set_folders()
    print_log("\n===================={}====================\n".format(session_name))
    print_log("$ folders read: {}".format(len(data_folder.list_folders)))

    for current_folder in data_folder.list_folders:
        apply_registrations_to_current_folder(
            current_folder, current_param, data_folder, current_session, file_name
        )

    del data_folder
