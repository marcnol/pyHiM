#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Purpose: Segments Masks in 3D using AI (stardist)

steps:
    - iterate over rois
    - load 3D file for cycle <i>
    - load global alignment for this cycle
    - re-align 3D image using XY alignment
    - segment objects to get labeled masks

    - display results:
    - output results
"""

# =============================================================================
# IMPORTS
# =============================================================================

import glob
import os
from datetime import datetime

import numpy as np
from skimage import io
from skimage.measure import regionprops

from fileProcessing.fileManagement import (
    Folders,
    get_dictionary_value,
    load_alignment_dictionary,
    print_dict,
    print_log,
    retrieve_number_rois_folder,
    write_string_to_file,
)
from imageProcessing.imageProcessing import (
    _segment_3d_masks,
    apply_xy_shift_3d_images,
    plot_raw_images_and_labels,
    reinterpolate_z,
)

# =============================================================================
# CLASSES
# =============================================================================


class SegmentMasks3D:
    def __init__(self, param, current_session, parallel=False):
        self.current_param = param
        self.current_session = current_session
        self.parallel = parallel

        self.p = {}
        self.roi_list = []
        self.number_rois = None
        self.dict_shifts = None
        self.dict_shifts_available = None
        self.filenames_to_process_list = []
        self.inner_parallel_loop = None
        self.data_folder = None
        self.current_folder = ""
        self.label = ""

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
            default="/mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks",
        ).rstrip("/")

        self.p["stardist_network"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "stardist_network",
            default="stardist_18032021_single_loci",
        )

        self.p["stardist_network3D"] = get_dictionary_value(
            self.current_param.param_dict["segmentedObjects"],
            "stardist_network3D",
            default="stardist_20210625_deconvolved",
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

    def plot_image_3d(self, image_3d, masks, normalize=False):
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
        
        fig1 = plot_raw_images_and_labels(image_3d, masks)

        return fig1

    def _segment_3d_volumes(self, image_3d_aligned):
        p = self.p

        binary, segmented_image_3d = _segment_3d_masks(
            image_3d_aligned,
            axis_norm=(0, 1, 2),
            pmin=1,
            pmax=99.8,
            model_dir=p["stardist_basename"],
            model_name=p["stardist_network3D"],
        )

        return binary, segmented_image_3d

    def segment_masks_3d_file(self, filename_to_process):

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

        # load  and preprocesses 3D fiducial file
        print_log("\n\n>>>Processing roi:[{}] cycle:[{}]<<<".format(roi, label))
        print_log("$ File:{}".format(os.path.basename(filename_to_process)))
        image_3d_0 = io.imread(filename_to_process).squeeze()

        # reinterpolates image in z if necessary
        image_3d = reinterpolate_z(
            image_3d_0, range(0, image_3d_0.shape[0], p["zBinning"]), mode="remove"
        )

        # drifts 3D stack in XY
        if self.dict_shifts_available:
            # uses existing shift calculated by align_images
            try:
                shift = self.dict_shifts["ROI:" + roi][label]
                print_log("> Applying existing XY shift...")
            except KeyError:
                shift = None
                raise SystemExit(
                    "# Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(
                        "ROI:" + roi, label
                    )
                )

        # applies XY shift to 3D stack
        if label != p["referenceBarcode"]:
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

        number_masks = np.max(segmented_image_3d)
        print_log("$ Number of masks detected: {}".format(number_masks))

        if number_masks > 0:
            output_extension = "_3Dmasks"
            npy_labeled_image_filename = (
                self.data_folder.output_folders["segmentedObjects"]
                + os.sep
                + os.path.basename(filename_to_process)
            )
            npy_labeled_image_filename = (
                npy_labeled_image_filename.split(".")[0]
                + "."
                + output_extension
                + ".npy"
            )
            print_log(
                " > Saving output labeled image: {}".format(npy_labeled_image_filename)
            )
            np.save(npy_labeled_image_filename, segmented_image_3d)

            # represents image in 3D with localizations
            print_log("> plotting outputs...")

            figures = []
            figures.append(
                [
                    self.plot_image_3d(image_3d_aligned, segmented_image_3d,),
                    output_extension + ".png",
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

    def segment_masks_3d_in_folder(self):
        """
        Segments 3D Masks in all files in root_folder

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
        print_log("$ Detected {} rois".format(self.number_rois))
        print_log(
            "$ Images to be processed: {}".format(self.current_param.files_to_process)
        )
        print_log(
            "$ Number of images to be processed: {}".format(
                len(self.current_param.files_to_process)
            )
        )

        # loads dicShifts with shifts for all rois and all labels
        self.dict_shifts, self.dict_shifts_available = load_alignment_dictionary(
            self.data_folder
        )

        if self.number_rois > 0:

            # loops over rois
            for roi in self.roi_list:

                # loads reference fiducial image for this ROI
                self.filenames_to_process_list = [
                    x
                    for x in self.current_param.files_to_process
                    if self.current_param.decode_file_parts(os.path.basename(x))["roi"]
                    == roi
                    and (
                        "DAPI"
                        in self.current_param.decode_file_parts(os.path.basename(x))[
                            "cycle"
                        ]
                        or "mask"
                        in self.current_param.decode_file_parts(os.path.basename(x))[
                            "cycle"
                        ]
                    )
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

                self.inner_parallel_loop = True
                # processes files in this ROI
                for file_index, filename_to_process in enumerate(
                    self.filenames_to_process_list
                ):
                    print_log(
                        "\n\n>>>Iteration: {}/{}<<<".format(file_index, n_files_to_process)
                    )
                    self.segment_masks_3d_file(filename_to_process)

        print_log("$ segmentMasks3D procesing time: {}".format(datetime.now() - now))

    def segment_masks_3d(self):
        """
        segments 3D masks in root_folder

        Returns
        -------
        None.

        """
        session_name = "segmentMasks3D"

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

        print_log("> Processing Folder: {}".format(self.current_folder))

        self.segment_masks_3d_in_folder()

        self.current_session.add(self.current_folder, session_name)

        print_log("$ segmentedObjects run in {} finished".format(self.current_folder))

        return 0

# =============================================================================
# FUNCTIONS
# =============================================================================

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
        peak0 = [x.max_intensity for x in properties]
        peak_list = peak0.copy()
        peak_list.sort()
        last2keep = np.min([n_tolerance, len(peak_list)])
        highest_peak_value = peak_list[-last2keep]
        selection = list(np.nonzero(peak0 > highest_peak_value)[0])

        # attributes properties using the brightests spots selected
        peak = [properties[x].max_intensity for x in selection]
        centroids = [properties[x].weighted_centroid for x in selection]
        sharpness = [
            float(properties[x].filled_area / properties[x].bbox_area)
            for x in selection
        ]
        roundness1 = [properties[x].equivalent_diameter for x in selection]
        roundness2 = [properties[x].extent for x in selection]
        npix = [properties[x].area for x in selection]
        sky = [0.0 for x in selection]
        peak = [properties[x].max_intensity for x in selection]
        flux = [
            100 * properties[x].max_intensity / threshold for x in selection
        ]  # peak intensity over the detection threshold
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
