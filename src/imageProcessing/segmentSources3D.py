#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 15:41:58 2021

@author: marcnol

Purpose: Corrects drift in 3D

Often there is drift in the z-position from cycle to cycle.

The drift correction routines take care of the corrections in XY but not in Z.

steps:
    - iterate over rois
    - load 3D file for cycle <i>
    - substract background
    - load global alignment for this cycle
    - re-align 3D image using XY alignment
    - segment objects to get labeled masks
    - Get weighted moments and gaussian fits

    - display results:
        - overlap XY, XZ, YZ image projections with moment and gaussian fits
    - output results in a Table() using same formatting as that used in segmentMasks.py

"""

# =============================================================================
# IMPORTS
# =============================================================================

import glob
import os
import uuid
from datetime import datetime

import numpy as np
from apifish.detection.spot_modeling import fit_subpixel
from astropy.table import Table, vstack
from skimage import exposure, io
from skimage.measure import regionprops

from fileProcessing.fileManagement import (
    Folders,
    get_dictionary_value,
    load_alignment_dictionary,
    print_dict,
    print_log,
    retrieve_number_rois_folder,
    try_get_client,
    write_string_to_file,
)
from imageProcessing.imageProcessing import (
    _plot_image_3d,
    _reinterpolate_focal_plane,
    _segment_3d_volumes_by_thresholding,
    _segment_3d_volumes_stardist,
    apply_xy_shift_3d_images,
    image_adjust,
    preprocess_3d_image,
    reinterpolate_z,
)

# =============================================================================
# CLASSES
# =============================================================================


class SegmentSources3D:
    def __init__(self, param, current_session, parallel=False):
        self.current_param = param
        self.current_session = current_session
        self.parallel = parallel

        self.p = {}
        self.roi_list = []
        self.number_rois = None
        self.dict_shifts = None
        self.dict_shifts_available = None
        self.roi = None
        self.filenames_to_process_list = []
        self.inner_parallel_loop = None
        self.data_folder = None
        self.current_folder = None
        self.label = None
        self.output_filename = None

        # parameters from infoList.json
        self.p["referenceBarcode"] = self.current_param.param_dict["alignImages"][
            "referenceFiducial"
        ]
        self.p["brightest"] = self.current_param.param_dict["segmentedObjects"][
            "brightest"
        ]
        self.p["blockSizeXY"] = self.current_param.param_dict["zProject"]["blockSize"]
        self.p["regExp"] = self.current_param.param_dict["acquisition"][
            "fileNameRegExp"
        ]
        self.p["zBinning"] = get_dictionary_value(
            self.current_param.param_dict["acquisition"], "zBinning", default=1
        )

        self.p["zWindow"] = int(
            self.current_param.param_dict["zProject"]["zwindows"] / self.p["zBinning"]
        )
        self.p["pixelSizeXY"] = self.current_param.param_dict["acquisition"][
            "pixelSizeXY"
        ]
        self.p["pixelSizeZ"] = self.current_param.param_dict["acquisition"][
            "pixelSizeZ"
        ]

        self.p["parallelizePlanes"] = get_dictionary_value(
            self.current_param.param_dict["acquisition"], "parallelizePlanes", default=1
        )

        # decides what segmentation method to use
        self.p["3Dmethod"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "3Dmethod",
            default="thresholding",
        )
        self.p["reducePlanes"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "reducePlanes",
            default=True,
        )

        # parameters used for 3D segmentation and deblending
        self.p["threshold_over_std"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "3D_threshold_over_std",
            default=1,
        )
        self.p["sigma"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"], "3D_sigma", default=3
        )
        box_size = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"], "3D_boxSize", default=32
        )
        self.p["boxSize"] = (box_size, box_size)
        filter_size = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "3D_filter_size",
            default=3,
        )
        self.p["filter_size"] = (filter_size, filter_size)
        self.p["area_min"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"], "3D_area_min", default=3
        )
        self.p["area_max"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "3D_area_max",
            default=1000,
        )
        self.p["nlevels"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"], "3D_nlevels", default=64
        )
        self.p["contrast"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "3D_contrast",
            default=0.001,
        )

        # parameters for stardist
        self.p["stardist_basename"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "stardist_basename",
            default="/mnt/PALM_dataserv/DATA/JB/2021/Data_single_loci/Annotated_data/data_loci_small/models/",
        )
        self.p["stardist_network"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "stardist_network3D",
            default="stardist_18032021_single_loci",
        )

        # parameters used for 3D gaussian fitting
        self.p["voxel_size_z"] = float(1000 * self.p["pixelSizeZ"] * self.p["zBinning"])
        self.p["voxel_size_yx"] = float(1000 * self.p["pixelSizeXY"])
        self.p["psf_z"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"], "3D_psf_z", default=500
        )
        self.p["psf_yx"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"], "3D_psf_yx", default=200
        )

        # range used for adjusting image levels during pre-precessing
        self.p["lower_threshold"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "3D_lower_threshold",
            default=0.9,
        )
        self.p["higher_threshold"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "3D_higher_threshold",
            default=0.9999,
        )

        # parameters used for plotting 3D image
        # sets the number of planes around the center of the image used to represent localizations in XZ and ZY
        self.p["windowDisplay"] = 10

    def plot_image_3d(self, image_3d, localizations=None, masks=None, normalize=None):
        """
        makes list with XY, XZ and ZY projections and sends for plotting

        Parameters
        ----------
        image_3d : numpy array
            image in 3D.
        localizations : list
            list of localizations to overlay onto 3D image. The default is None.

        Returns
        -------
        figure handle

        """
        window = self.p["windowDisplay"]
        fig1 = _plot_image_3d(
            image_3d,
            localizations=localizations,
            masks=masks,
            normalize=normalize,
            window=window,
        )
        return fig1

    def _segment_3d_volumes(self, image_3d_aligned):
        p = self.p

        if "stardist" in p["3Dmethod"]:
            binary, segmented_image_3d = _segment_3d_volumes_stardist(
                image_3d_aligned,
                deblend_3d=True,
                axis_norm=(0, 1, 2),
                model_dir=p["stardist_basename"],
                model_name=p["stardist_network"],
            )
        else:
            binary, segmented_image_3d = _segment_3d_volumes_by_thresholding(
                image_3d_aligned,
                threshold_over_std=p["threshold_over_std"],
                sigma=p["sigma"],
                box_size=p["boxSize"],
                filter_size=p["filter_size"],
                area_min=p["area_min"],
                area_max=p["area_max"],
                nlevels=p["nlevels"],
                contrast=p["contrast"],
                deblend_3d=True,
                parallel_execution=self.inner_parallel_loop,
            )

        return binary, segmented_image_3d

    def segment_sources_3d_file(self, filename_to_process):

        p = self.p
        # excludes the reference fiducial and processes files in the same ROI
        roi = self.current_param.decode_file_parts(
            os.path.basename(filename_to_process)
        )["roi"]
        label = str(
            self.current_param.decode_file_parts(os.path.basename(filename_to_process))[
                "cycle"
            ]
        )

        # creates Table that will hold results
        output_table = create_output_table()

        # - load  and preprocesses 3D fiducial file
        print_log("\n\n>>>Processing roi:[{}] cycle:[{}]<<<".format(roi, label))
        print_log("$ File:{}".format(os.path.basename(filename_to_process)))
        image_3d_0 = io.imread(filename_to_process).squeeze()

        # reinterpolates image in z if necessary
        image_3d_0 = reinterpolate_z(
            image_3d_0, range(0, image_3d_0.shape[0], p["zBinning"]), mode="remove"
        )

        # restricts analysis to a sub volume containing sources
        if p["reducePlanes"]:
            _, z_range, _ = _reinterpolate_focal_plane(
                image_3d_0, block_size_xy=p["blockSizeXY"], window=p["zWindow"]
            )
            z_offset = z_range[1][0]
            image_3d = image_3d_0[z_range[1], :, :].copy()
            print_log(
                "$ Focal plane found: {}, z_range = {}, image_size = {}".format(
                    z_range[0], z_range[1], image_3d.shape
                )
            )
        else:
            image_3d = image_3d_0.copy()
            z_offset = 0
            print_log("$ z_range used = 0-{}".format(image_3d.shape[0]))

        # preprocesses image by background substraction and level normalization
        if "stardist" not in p["3Dmethod"]:
            image_3d = preprocess_3d_image(
                image_3d,
                p["lower_threshold"],
                p["higher_threshold"],
                parallel_execution=self.inner_parallel_loop,
            )

        # drifts 3D stack in XY
        shift = None
        if self.dict_shifts_available and label != p["referenceBarcode"]:
            # uses existing shift calculated by align_images
            try:
                shift = self.dict_shifts["ROI:" + roi][label]
                print_log("> Applying existing XY shift...")
            except KeyError:
                shift = None

        if shift is None and label != p["referenceBarcode"]:
            raise SystemExit(
                "> Existing with ERROR: Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(
                    "ROI:" + self.roi, label
                )
            )

        # applies XY shift to 3D stack
        if label != p["referenceBarcode"]:
            # print_log("$ Applies shift = {:.2f}".format(shift))
            print_log("$ Applies shift = [{:.2f} ,{:.2f}]".format(shift[0], shift[1]))
            image_3d_aligned = apply_xy_shift_3d_images(
                image_3d, shift, parallel_execution=self.inner_parallel_loop
            )
        else:
            print_log("$ Running reference fiducial cycle: no shift applied!")
            shift = np.array([0.0, 0.0])
            image_3d_aligned = image_3d

        # segments 3D volumes
        _, segmented_image_3d = self._segment_3d_volumes(image_3d_aligned)

        # gets centroids and converts to spot int64 NPY array
        (
            spots,
            sharpness,
            roundness1,
            roundness2,
            npix,
            sky,
            peak,
            flux,
            mag,
        ) = get_mask_properties(
            segmented_image_3d,
            image_3d_aligned,
            threshold=p["threshold_over_std"],
            n_tolerance=p["brightest"],
        )

        number_sources = len(peak)
        print_log(
            "$ Number of sources detected by image segmentation: {}".format(
                number_sources
            )
        )

        if number_sources > 0:
            print_log("> Refits spots using gaussian 3D fittings...")

            print_log(" > Rescales image values after reinterpolation")
            image_3d_aligned = exposure.rescale_intensity(
                image_3d_aligned, out_range=(0, 1)
            )  # removes negative intensity levels

            # calls bigfish to get 3D sub-pixel coordinates based on 3D gaussian fitting
            # compatibility with latest version of bigfish. To be removed if stable.
            try:
                # version 0.4 commit fa0df4f
                spots_subpixel = fit_subpixel(
                    image_3d_aligned,
                    spots,
                    voxel_size_z=p["voxel_size_z"],
                    voxel_size_yx=p["voxel_size_yx"],
                    psf_z=p["psf_z"],
                    psf_yx=p["psf_yx"],
                )
            except TypeError:
                # version > 0.5
                spots_subpixel = fit_subpixel(
                    image_3d_aligned,
                    spots,
                    (
                        p["voxel_size_z"],
                        p["voxel_size_yx"],
                        p["voxel_size_yx"],
                    ),  # voxel size
                    (p["psf_z"], p["psf_yx"], p["psf_yx"]),
                )  # spot radius

            print_log(" > Updating table and saving results")
            # updates table
            for i in range(spots_subpixel.shape[0]):
                z, x, y = spots_subpixel[i, :]
                table_entry = [
                    str(uuid.uuid4()),
                    roi,
                    0,
                    int(label.split("RT")[1]),
                    i,
                    z + z_offset,
                    y,
                    x,
                    sharpness[i],
                    roundness1[i],
                    roundness2[i],
                    npix[i],
                    sky[i],
                    peak[i],
                    flux[i],
                    mag[i],
                ]
                output_table.add_row(table_entry)

            # represents image in 3D with localizations
            figures = []
            if "stardist" in p["3Dmethod"]:
                image_3d_aligned = image_adjust(
                    image_3d_aligned,
                    lower_threshold=p["lower_threshold"],
                    higher_threshold=p["higher_threshold"],
                )[0]

            figures.append(
                [
                    self.plot_image_3d(
                        image_3d_aligned,
                        masks=segmented_image_3d,
                        localizations=[spots_subpixel, spots],
                    ),
                    "_3DimageNlocalizations.png",
                ]
            )

            # saves figures
            output_filenames = [
                self.data_folder.output_folders["segmentedObjects"]
                + os.sep
                + os.path.basename(filename_to_process)
                + x[1]
                for x in figures
            ]

            for fig, file in zip(figures, output_filenames):
                fig[0].savefig(file)

        del image_3d_aligned, image_3d, image_3d_0

        return output_table

    def segment_sources_3d_in_folder(self):
        """
        Fits sources in all files in root_folder

        Returns
        -------
        None.

        """
        now = datetime.now()

        # Reads list of parameters assigned upon Class initialization
        p = self.p
        print_dict(p)

        # Finds images to process
        files_folder = glob.glob(self.current_folder + os.sep + "*.tif")
        self.current_param.find_files_to_process(files_folder)
        self.roi_list = retrieve_number_rois_folder(
            self.current_folder, p["regExp"], ext="tif"
        )
        self.number_rois = len(self.roi_list)
        print_log("\n$ Detected {} rois".format(self.number_rois))
        print_log(
            "$ Number of images to be processed: {}".format(
                len(self.current_param.files_to_process)
            )
        )

        # loads dicShifts with shifts for all rois and all labels
        self.dict_shifts, self.dict_shifts_available = load_alignment_dictionary(
            self.data_folder
        )

        # creates Table that will hold results
        output_table_global = create_output_table()
        output_tables = []

        if self.p["parallelizePlanes"]:
            client = None
        else:
            client = try_get_client()

        if self.number_rois > 0:

            # loops over rois
            for roi in self.roi_list:
                # loads reference fiducial image for this ROI
                self.roi = roi

                self.filenames_to_process_list = [
                    x
                    for x in self.current_param.files_to_process
                    if self.current_param.decode_file_parts(os.path.basename(x))["roi"]
                    == roi
                    and "RT"
                    in self.current_param.decode_file_parts(os.path.basename(x))[
                        "cycle"
                    ]
                ]

                n_files_to_process = len(self.filenames_to_process_list)
                print_log("$ Found {} files in ROI [{}]".format(n_files_to_process, roi))
                print_log(
                    "$ [roi:cycle] {}".format(
                        " | ".join(
                            [
                                str(
                                    self.current_param.decode_file_parts(
                                        os.path.basename(x)
                                    )["roi"]
                                )
                                + ":"
                                + str(
                                    self.current_param.decode_file_parts(
                                        os.path.basename(x)
                                    )["cycle"]
                                )
                                for x in self.filenames_to_process_list
                            ]
                        )
                    )
                )

                if client is None:
                    self.inner_parallel_loop = True

                    for file_index, filename_to_process in enumerate(
                        self.filenames_to_process_list
                    ):  # self.current_param.files_to_process):
                        print_log(
                            "\n\n>>>Iteration: {}/{}<<<".format(
                                file_index, n_files_to_process
                            )
                        )

                        output_tables.append(
                            self.segment_sources_3d_file(filename_to_process)
                        )
                else:
                    self.inner_parallel_loop = False
                    print_log(
                        "> Aligning {} files using {} workers...".format(
                            n_files_to_process, len(client.scheduler_info()["workers"])
                        )
                    )

                    futures = [
                        client.submit(self.segment_sources_3d_file, x)
                        for x in self.filenames_to_process_list
                    ]

                    output_tables = client.gather(futures)
                    print_log(
                        " > Retrieving {} results from cluster".format(
                            len(output_tables)
                        )
                    )

                # Merges Tables for different cycles and appends results Table to that of previous ROI
                output_table_global = vstack([output_table_global] + output_tables)

        print_log("$ segmentSources3D procesing time: {}".format(datetime.now() - now))

        # saves Table with all shifts in every iteration to avoid loosing computed data
        output_table_global.write(
            self.output_filename, format="ascii.ecsv", overwrite=True,
        )

    def segment_sources_3d(self):
        """
        runs 3D fitting routine in root_folder

        Returns
        -------
        None.

        """
        session_name = "segmentSources3D"

        # processes folders and files
        print_log("\n===================={}====================\n".format(session_name))
        self.data_folder = Folders(self.current_param.param_dict["rootFolder"])
        print_log("$ folders read: {}".format(len(self.data_folder.list_folders)))
        write_string_to_file(
            self.current_param.param_dict["fileNameMD"],
            "## {}\n".format(session_name),
            "a",
        )

        # creates output folders and filenames
        self.current_folder = self.data_folder.list_folders[0]

        self.data_folder.create_folders(self.current_folder, self.current_param)
        self.label = self.current_param.param_dict["acquisition"]["label"]
        self.output_filename = (
            self.data_folder.output_files["segmentedObjects"]
            + "_3D_"
            + self.label
            + ".dat"
        )

        print_log("> Processing Folder: {}".format(self.current_folder))

        self.segment_sources_3d_in_folder()

        self.current_session.add(self.current_folder, session_name)

        print_log("$ segmentedObjects run in {} finished".format(self.current_folder))

        return 0

##########################################
# FUNCTIONS
##########################################

def get_mask_properties(segmented_image_3d, image_3d_aligned, threshold=10, n_tolerance=1000):
    """
    get object properties from labeled image and formats
    centroids in NPY array

    Parameters
    ----------
    segmented_image_3d : NPY 3D array
        labeled 3D image.
    image_3d_aligned : NPY 3D array
        pre-processed 3D image.

    Returns
    -------
    spots : NPY int64 array
        list of spots with the format: zyx

    """

    # gets object properties
    properties = regionprops(segmented_image_3d, intensity_image=image_3d_aligned)

    if len(properties) > 0:
        # selects n_tolerance brightest spots and keeps only these for further processing
        try:
            # compatibility with scikit_image versions <= 0.18
            peak0 = [x.max_intensity for x in properties]
        except AttributeError:
            # compatibility with scikit_image versions >=0.19
            peak0 = [x.intensity_max for x in properties]

        peak_list = peak0.copy()
        peak_list.sort()

        if n_tolerance == "None":
            last2keep = len(peak_list)
        else:
            last2keep = np.min([n_tolerance,len(peak_list)])

        highest_peak_value = peak_list[-last2keep]
        selection = list(np.nonzero(peak0 > highest_peak_value)[0])

        # attributes properties using the brightests spots selected
        try:
            # compatibility with scikit_image versions <= 0.18
            peak = [properties[x].max_intensity for x in selection]
            centroids = [properties[x].weighted_centroid for x in selection]
            sharpness = [
                float(properties[x].filled_area / properties[x].bbox_area)
                for x in selection
            ]
            roundness1 = [properties[x].equivalent_diameter for x in selection]
        except AttributeError:
            # compatibility with scikit_image versions >=0.19
            peak = [properties[x].intensity_max for x in selection]
            centroids = [properties[x].centroid_weighted for x in selection]
            sharpness = [
                float(properties[x].area_filled / properties[x].area_bbox)
                for x in selection
            ]
            roundness1 = [properties[x].equivalent_diameter_area for x in selection]

        roundness2 = [properties[x].extent for x in selection]
        npix = [properties[x].area for x in selection]
        sky = [0.0 for x in selection]

        try:
            # compatibility with scikit_image versions <= 0.18
            peak = [properties[x].max_intensity for x in selection]
            flux = [
                100 * properties[x].max_intensity / threshold for x in selection
            ]  # peak intensity over t$
        except AttributeError:
            # compatibility with scikit_image versions >=0.19
            peak = [properties[x].intensity_max for x in selection]
            flux = [
                100 * properties[x].intensity_max / threshold for x in selection
            ]  # peak intensity$

        mag = [-2.5 * np.log10(x) for x in flux]  # -2.5 log10(flux)

        # converts centroids to spot coordinates for bigfish to run 3D gaussian fits
        z = [x[0] for x in centroids]
        y = [x[1] for x in centroids]
        x = [x[2] for x in centroids]

        spots = np.zeros((len(z), 3))
        spots[:, 0] = z
        spots[:, 1] = y
        spots[:, 2] = x
        spots = spots.astype("int64")

        return (
            spots,
            sharpness,
            roundness1,
            roundness2,
            npix,
            sky,
            peak,
            flux,
            mag,
        )
    else:
        # creates output lists to return
        return [], [], [], [], [], [], [], [], []

def create_output_table():
    output = Table(
        names=(
            "Buid",
            "ROI #",
            "CellID #",
            "Barcode #",
            "id",
            "zcentroid",
            "xcentroid",
            "ycentroid",
            "sharpness",
            "roundness1",
            "roundness2",
            "npix",
            "sky",
            "peak",
            "flux",
            "mag",
        ),
        dtype=(
            "S2",
            "int",
            "int",
            "int",
            "int",
            "f4",
            "f4",
            "f4",
            "f4",
            "f4",
            "f4",
            "int",
            "f4",
            "f4",
            "f4",
            "f4",
        ),
    )
    return output
