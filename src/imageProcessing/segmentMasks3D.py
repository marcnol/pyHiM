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

import glob
import os
from datetime import datetime

import numpy as np
from skimage import io

from core.data_manager import create_folder
from core.parameters import AcquisitionParams, SegmentationParams, load_alignment_dict
from core.pyhim_logging import print_log, print_session_name, write_string_to_file
from core.saving import plot_raw_images_and_labels
from imageProcessing.alignImages import apply_xy_shift_3d_images
from imageProcessing.makeProjections import reinterpolate_z
from imageProcessing.segmentMasks import _segment_3d_masks

# =============================================================================
# CLASSES
# =============================================================================


class Mask3D:
    def __init__(
        self,
        param,
        parallel=False,
    ):
        self.current_param = param
        self.parallel = parallel

        self.dict_shifts = None
        self.dict_shifts_available = None
        self.filenames_to_process_list = []
        self.inner_parallel_loop = None

    def _segment_3d_volumes(self, image_3d_aligned, seg_params: SegmentationParams):
        if seg_params.stardist_basename is not None and os.path.exists(
            seg_params.stardist_basename
        ):
            base_dir = seg_params.stardist_basename
        else:
            base_dir = os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                os.pardir,
                "stardist_models",
            )

        if seg_params.stardist_network3D is not None and os.path.exists(
            os.path.join(base_dir, seg_params.stardist_network3D)
        ):
            model_name = seg_params.stardist_network3D
        else:
            model_name = "DAPI_3D_stardist_17032021_deconvolved"
        binary, segmented_image_3d = _segment_3d_masks(
            image_3d_aligned,
            axis_norm=(0, 1, 2),
            pmin=1,
            pmax=99.8,
            model_dir=base_dir,
            model_name=model_name,
        )

        return binary, segmented_image_3d

    def segment_masks_3d_file(
        self,
        filename_to_process,
        data_path,
        seg_params,
        roi_name,
        acq_params: AcquisitionParams,
        reference_fiducial: str,
    ):
        # excludes the reference fiducial and processes files in the same ROI
        label = str(
            self.current_param.decode_file_parts(os.path.basename(filename_to_process))[
                "cycle"
            ]
        )

        # load  and preprocesses 3D fiducial file
        print_log(f"\n\n>>>Processing roi:[{roi_name}] cycle:[{label}]<<<")
        print_log(f"$ File:{os.path.basename(filename_to_process)}")
        image_3d_0 = io.imread(filename_to_process).squeeze()

        # reinterpolates image in z if necessary
        image_3d = reinterpolate_z(
            image_3d_0,
            range(0, image_3d_0.shape[0], acq_params.zBinning),
            mode="remove",
        )

        # # drifts 3D stack in XY
        # if self.dict_shifts_available:
        # uses existing shift calculated by align_images
        try:
            shift = self.dict_shifts[f"ROI:{roi_name}"][label]
            print_log("> Applying existing XY shift...")
        except KeyError as e:
            shift = None
            raise SystemExit(
                f"# Could not find dictionary with alignment parameters for this ROI: ROI:{roi_name}, label: {label}"
            ) from e

        # applies XY shift to 3D stack
        if label != reference_fiducial:
            print_log(f"$ Applies shift = [{shift[0]:.2f} ,{shift[1]:.2f}]")
            image_3d_aligned = apply_xy_shift_3d_images(
                image_3d, shift, parallel_execution=self.inner_parallel_loop
            )
        else:
            print_log("$ Running reference fiducial cycle: no shift applied!")
            shift = np.array([0.0, 0.0])
            image_3d_aligned = image_3d

        # segments 3D volumes
        _, segmented_image_3d = self._segment_3d_volumes(image_3d_aligned, seg_params)

        number_masks = np.max(segmented_image_3d)
        print_log(f"$ Number of masks detected: {number_masks}")

        if number_masks > 0:
            output_extension = {"2D": "_Masks", "3D": "_3Dmasks"}
            npy_labeled_image_base_2d = (
                data_path
                + os.sep
                + seg_params.mask_2d_folder
                + os.sep
                + "data"
                + os.sep
            )
            npy_labeled_image_filename_2d = (
                npy_labeled_image_base_2d
                + os.path.basename(filename_to_process).split(".")[0]
                + output_extension["2D"]
                + ".npy"
            )
            npy_labeled_image_base_3d = (
                data_path
                + os.sep
                + seg_params.mask_3d_folder
                + os.sep
                + "data"
                + os.sep
                + os.path.basename(filename_to_process)
            )
            npy_labeled_image_filename_3d = (
                npy_labeled_image_base_3d.split(".")[0]
                + output_extension["3D"]
                + ".npy"
            )
            print_log("> Saving output labeled images:")
            print_log(f"\t2D:{npy_labeled_image_filename_2d}")
            print_log(f"\t3D:{npy_labeled_image_filename_3d}")

            # saves 3D image
            np.save(npy_labeled_image_filename_3d, segmented_image_3d)

            # saves 2D image
            segmented_image_2d = np.max(segmented_image_3d, axis=0)
            create_folder(npy_labeled_image_base_2d)
            np.save(npy_labeled_image_filename_2d, segmented_image_2d)

            # represents image in 3D with localizations
            print_log("> plotting outputs...")

            # figures = []
            # figures.append(
            #     [
            fig = plot_image_3d(image_3d_aligned, segmented_image_3d)
            end_file = output_extension["3D"] + ".png"
            #     ]
            # )

            # saves figures
            # output_filenames = [
            output_filename = (
                data_path
                + os.sep
                + seg_params.mask_3d_folder
                + os.sep
                + os.path.basename(filename_to_process)
                + end_file
            )
            #     + x[1]
            #     for x in figures
            # ]

            # for fig, file in zip(figures, output_filenames):
            # fig[0].savefig(file)
            fig.savefig(output_filename)

        del image_3d_aligned, image_3d, image_3d_0

    def segment_masks_3d_in_folder(
        self,
        roi_name: str,
        data_path,
        dict_shifts_path,
        seg_params,
        acq_params: AcquisitionParams,
        reference_fiducial,
    ):
        """
        Segments 3D Masks in all files in root_folder

        Returns
        -------
        None.

        """
        now = datetime.now()

        # Reads list of parameters assigned upon Class initialization

        # Finds images to process
        files_folder = glob.glob(data_path + os.sep + "*.tif")
        self.current_param.find_files_to_process(files_folder)
        print_log(f"$ Images to be processed: {self.current_param.files_to_process}")
        nb_imgs = len(self.current_param.files_to_process)
        print_log(f"$ Number of images to be processed: {nb_imgs}")

        # loads dicShifts with shifts for all rois and all labels
        self.dict_shifts, self.dict_shifts_available = load_alignment_dict(
            dict_shifts_path
        )

        # loads reference fiducial image for this ROI
        self.filenames_to_process_list = [
            x
            for x in self.current_param.files_to_process
            if self.current_param.decode_file_parts(os.path.basename(x))["roi"]
            == roi_name
            and (
                "DAPI"
                in self.current_param.decode_file_parts(os.path.basename(x))["cycle"]
                or "mask"
                in self.current_param.decode_file_parts(os.path.basename(x))["cycle"]
            )
        ]
        n_files_to_process = len(self.filenames_to_process_list)
        print_log(f"$ Found {n_files_to_process} files in ROI [{roi_name}]")
        print_log(
            "$ [roi:cycle] {}".format(
                " | ".join(
                    [
                        str(
                            self.current_param.decode_file_parts(os.path.basename(x))[
                                "roi"
                            ]
                        )
                        + ":"
                        + str(
                            self.current_param.decode_file_parts(os.path.basename(x))[
                                "cycle"
                            ]
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
            print_log(f"\n\n>>>Iteration: {file_index+1}/{n_files_to_process}<<<")
            self.segment_masks_3d_file(
                filename_to_process,
                data_path,
                seg_params,
                roi_name,
                acq_params,
                reference_fiducial,
            )

        print_log(f"$ mask_3d procesing time: {datetime.now() - now}")

    def segment_masks_3d(
        self,
        roi_name: str,
        data_path,
        dict_shifts_path,
        seg_params,
        acq_params,
        reference_fiducial,
    ):
        """
        segments 3D masks in root_folder

        Returns
        -------
        None.

        """
        session_name = "mask_3d"

        # processes folders and files

        print_session_name(session_name)
        write_string_to_file(
            self.current_param.param_dict["fileNameMD"],
            f"## {session_name}\n",
            "a",
        )

        print_log(f"> Processing Folder: {data_path}")

        self.segment_masks_3d_in_folder(
            roi_name,
            data_path,
            dict_shifts_path,
            seg_params,
            acq_params,
            reference_fiducial,
        )

        print_log(f"$ segmentedObjects run in {data_path} finished")

        return 0


# =============================================================================
# FUNCTIONS
# =============================================================================


def plot_image_3d(image_3d, masks):
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

    return plot_raw_images_and_labels(image_3d, masks)
