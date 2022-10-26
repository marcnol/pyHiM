# -*- coding: utf-8 -*-
"""
Spyder Editor


Purpose: Corrects drift in 3D

The drift correction routines in alignImages.py take care of the corrections in XY but not in Z.

Drift in the z-position from cycle to cycle is small, typically 100-200 nm, but needs to be checked
and corrected. In addition, sample deformation can lead to inhomogeneous drift that cannot be corrected by a rigid body translation.

This routine solves both issues by correcting drift in 3D by block decomposition.

steps:
    - iterate over rois
    - load 3D fiducial file for reference fiducial
    - iterate over cycles <i>
    - load 3D fiducial file for fiducial barcode <i>
    - re-align 3D fiducial image using XY alignment
    - perform block alignment in 3D by cross=correlating blocks in 3D.

    - store in database.
    - display results:
        - drift correction maps in X-Y-Z
        - corrected blocks in XY, ZX, ZY

During buildMatrix this database is loaded, if available.
    - check database exist and load it
    - correct z-coordinate of the barcode provided the correction given in the dict

"""
# =============================================================================
# IMPORTS
# =============================================================================

import glob
import os
import resource
from datetime import datetime

import matplotlib.pylab as plt
import numpy as np
from astropy.table import Table, vstack

from pympler import tracker
from skimage import io
from skimage.registration import phase_cross_correlation

from fileProcessing.fileManagement import (
    Folders,
    get_dictionary_value,
    load_alignment_dictionary,
    print_dict,
    print_log,
    rt_to_filename,
    try_get_client,
    write_string_to_file,
)
from imageProcessing.imageProcessing import (
    apply_xy_shift_3d_images,
    combine_blocks_image_by_reprojection,
    image_block_alignment_3d,
    plot_3d_shift_matrices,
    plot_4_images,
    preprocess_3d_image,
    reinterpolate_z,
)

# =============================================================================
# CLASSES
# =============================================================================


class Drift3D:
    def __init__(self, param, current_session, parallel=False):
        self.current_param = param
        self.current_session = current_session
        self.window = 3
        self.parallel = parallel
        self.p = {}
        self.image_ref_0 = None
        self.image_ref = None
        self.filenames_to_process_list = []
        self.filenames_with_ref_barcode = []
        self.roi_list = []
        self.number_rois = None
        self.dict_shifts = None
        self.dict_shifts_available = None
        self.inner_parallel_loop = None
        self.data_folder = None
        self.current_folder = None
        self.output_filename = None

        self.p["blockSizeXY"] = 128
        self.p["upsample_factor"] = 100
        self.p["lower_threshold"] = 0.9
        self.p["higher_threshold"] = 0.9999

        self.p["blockSizeXY"] = get_dictionary_value(
            self.current_param.param_dict["alignImages"], "blockSizeXY", default=128
        )
        self.p["lower_threshold"] = get_dictionary_value(
            self.current_param.param_dict["alignImages"],
            "3D_lower_threshold",
            default=0.9,
        )
        self.p["higher_threshold"] = get_dictionary_value(
            self.current_param.param_dict["alignImages"],
            "3D_higher_threshold",
            default=0.9999,
        )

        self.p["axes2Plot"] = range(3)
        self.p["referenceBarcode"] = self.current_param.param_dict["alignImages"][
            "referenceFiducial"
        ]

        if "zBinning" in self.current_param.param_dict["acquisition"]:
            self.p["zBinning"] = self.current_param.param_dict["acquisition"][
                "zBinning"
            ]
        else:
            self.p["zBinning"] = 1

        if "parallelizePlanes" in self.current_param.param_dict["acquisition"]:
            self.p["parallelizePlanes"] = self.current_param.param_dict["acquisition"][
                "parallelizePlanes"
            ]
        else:
            self.p["parallelizePlanes"] = 1

    def find_file_to_process(self, n_barcode, n_roi):
        barcode = "RT" + str(n_barcode)
        roi = str(n_roi) + "_ROI"
        channelbarcode = self.current_param.set_channel("barcode_channel", "ch01")

        files_folder = glob.glob(self.data_folder.master_folder + os.sep + "*.tif")
        image_file = [
            x for x in files_folder if roi in x and barcode in x and channelbarcode in x
        ]

        return image_file

    def align_fiducials_3d_file(self, filename_to_process):
        """
        Aligns <filename_to_process> fiducial against reference

        Returns
        -------
        None.

        """

        p = self.p
        alignment_results_table = create_output_table()

        # excludes the reference fiducial and processes files in the same ROI
        roi = self.current_param.decode_file_parts(
            os.path.basename(filename_to_process)
        )["roi"]
        label = str(
            self.current_param.decode_file_parts(os.path.basename(filename_to_process))[
                "cycle"
            ]
        )
        inner_parallel_loop = self.inner_parallel_loop
        image_ref = self.image_ref
        image_ref_0 = self.image_ref_0
        dict_shifts_available = self.dict_shifts_available
        dict_shifts = self.dict_shifts
        output_folder = self.data_folder.output_folders["alignImages"]

        return _align_fiducials_3d_file(
            filename_to_process,
            alignment_results_table,
            p,
            roi,
            label,
            inner_parallel_loop,
            image_ref,
            image_ref_0,
            dict_shifts,
            dict_shifts_available,
            output_folder,
        )

    def load_reference_fiducial(self, filename_reference):
        """
        Loads Reference fiducial image and reports on the number of cycles to process for this reference

        Returns
        -------
        None.

        """

        self.p["fileNameReference"] = filename_reference
        self.p["ROI"] = self.roi_list[filename_reference]
        print_log("Loading reference 3D image: {}".format(filename_reference))

        self.image_ref_0, self.image_ref = load_n_preprocess_image(
            filename_reference,
            self.p["zBinning"],
            self.p["lower_threshold"],
            self.p["higher_threshold"],
            parallel_execution=False,
        )

        self.image_ref_0 = np.sum(
            self.image_ref_0, axis=0
        )  # replaces 3D with a 2D projection

        self.filenames_to_process_list = [
            x
            for x in self.current_param.files_to_process
            if (x not in filename_reference)
            and self.current_param.decode_file_parts(os.path.basename(x))["roi"]
            == self.p["ROI"]
        ]

        print_log(
            "$ Found {} files in ROI: {}".format(
                len(self.filenames_to_process_list), self.p["ROI"]
            )
        )
        print_log(
            "$ [roi:cycle] {}".format(
                "|".join(
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

    def load_dict_shifts(self):
        """
        Lods dictionary of XY shifts

        Returns
        -------
        None.

        """
        print_log("\nReference barcode: {}".format(self.p["referenceBarcode"]))
        self.filenames_with_ref_barcode, self.roi_list = rt_to_filename(
            self.current_param, self.p["referenceBarcode"]
        )

        self.number_rois = len(self.roi_list)
        print_log("\nDetected {} rois".format(self.number_rois))

        # loads dicShifts with shifts for all rois and all labels
        self.dict_shifts, self.dict_shifts_available = load_alignment_dictionary(
            self.data_folder
        )

    def align_fiducials_3d_in_folder(self):
        """
        Refits all the barcode files found in root_folder

        Returns
        -------
        None.

        """
        now = datetime.now()
        print_dict(self.p)
        tra = tracker.SummaryTracker()

        # gets files to process
        files_folder = glob.glob(self.current_folder + os.sep + "*.tif")
        self.current_param.find_files_to_process(files_folder)

        # loads dictinary of shifts
        self.load_dict_shifts()

        # creates Table that will hold results
        alignment_results_table_global = create_output_table()
        alignment_results_tables = []

        if self.p["parallelizePlanes"]:
            client = None
        else:
            client = try_get_client()

        if self.number_rois > 0:

            # loops over rois
            for filename_reference in self.filenames_with_ref_barcode:

                # loads reference fiducial image for this ROI
                self.load_reference_fiducial(filename_reference)
                number_files = len(self.filenames_to_process_list)
                tra.print_diff()

                # for file_index, filename_to_process in enumerate(self.current_param.files_to_process):
                if client is None:
                    self.inner_parallel_loop = True
                    for file_index, filename_to_process in enumerate(
                        self.filenames_to_process_list
                    ):

                        print_log(
                            "\n\n>>>Iteration: {}/{}<<<".format(
                                file_index, number_files
                            )
                        )

                        alignment_results_tables.append(
                            self.align_fiducials_3d_file(filename_to_process)
                        )

                        tra.print_diff()

                else:
                    self.inner_parallel_loop = False
                    print_log(
                        "> Aligning {} files using {} workers...".format(
                            number_files, len(client.scheduler_info()["workers"])
                        )
                    )

                    futures = [
                        client.submit(self.align_fiducials_3d_file, x)
                        for x in self.filenames_to_process_list
                    ]

                    alignment_results_tables = client.gather(futures)
                    print_log(
                        " > Retrieving {} results from cluster".format(
                            len(alignment_results_tables)
                        )
                    )

                    # del futures

                # Merges Tables for different cycles and appends results Table to that of previous ROI
                alignment_results_table_global = vstack(
                    [alignment_results_table_global] + alignment_results_tables
                )

        # saves Table with all shifts
        output_filename = (
            self.data_folder.output_files["alignImages"].split(".")[0] + "_block3D.dat"
        )
        alignment_results_table_global.write(
            output_filename, format="ascii.ecsv", overwrite=True,
        )

        print_log("$ alignImages3D procesing time: {}".format(datetime.now() - now))
        print_log(f"$ alignImages3D output Table saved in: {output_filename}")

    def align_fiducials_3d(self):
        """
        runs refitting routine in root_folder

        Returns
        -------
        None.

        """
        session_name = "alignImages3D"

        # processes folders and files
        print_log("\n===================={}====================\n".format(session_name))
        self.data_folder = Folders(self.current_param.param_dict["rootFolder"])
        print_log("folders read: {}".format(len(self.data_folder.list_folders)))
        write_string_to_file(
            self.current_param.param_dict["fileNameMD"],
            "## {}\n".format(session_name),
            "a",
        )

        # creates output folders and filenames
        self.current_folder = self.data_folder.list_folders[0]

        self.data_folder.create_folders(self.current_folder, self.current_param)
        self.output_filename = self.data_folder.output_files["alignImages"]

        print_log("-------> Processing Folder: {}".format(self.current_folder))
        # self.current_log.parallel = self.parallel

        self.align_fiducials_3d_in_folder()

        self.current_session.add(self.current_folder, session_name)

        print_log("HiM matrix in {} processed".format(self.current_folder))

        return 0


# =============================================================================
#   FUNCTIONS
# =============================================================================


def load_n_preprocess_image(
    filename_to_process,
    z_binning,
    lower_threshold,
    higher_threshold,
    parallel_execution=True,
):

    print_log("$ File:{}".format(os.path.basename(filename_to_process)))

    image_3d_0 = io.imread(filename_to_process).squeeze()

    # reinterpolates image in z if necessary
    image_3d_0 = reinterpolate_z(
        image_3d_0, range(0, image_3d_0.shape[0], z_binning), mode="remove"
    )

    image_3d = preprocess_3d_image(
        image_3d_0,
        lower_threshold,
        higher_threshold,
        parallel_execution=parallel_execution,
    )

    return image_3d_0, image_3d

def _align_fiducials_3d_file(
    filename_to_process,
    alignment_results_table,
    p,
    roi,
    label,
    inner_parallel_loop,
    image_ref,
    image_ref_0,
    dict_shifts,
    dict_shifts_available,
    output_folder,
):

    # - load  and preprocesses 3D fiducial file
    print_log("\n\n>>>Processing roi:[{}] cycle:[{}]<<<".format(roi, label))
    image_3d_0, image_3d = load_n_preprocess_image(
        filename_to_process,
        p["zBinning"],
        p["lower_threshold"],
        p["higher_threshold"],
        parallel_execution=inner_parallel_loop,
    )

    # shows original images and background substracted
    image_3d_0 = np.sum(image_3d_0, axis=0)  # replaces by a 2D projection
    images = [image_ref, image_3d]
    images_2d = [np.sum(x, axis=0) for x in images]
    fig1 = plot_4_images(
        [image_ref_0, image_3d_0] + images_2d,
        titles=["reference", "cycle <i>", "processed reference", "processed cycle <i>"],
    )

    del image_3d_0

    # drifts 3D stack in XY
    # ---------------------
    if dict_shifts_available:
        # uses existing shift calculated by align_images
        try:
            shift = dict_shifts["ROI:" + roi][label]
            print_log("> Applying existing XY shift...")
        except KeyError:
            shift = None
            print_log(
                "Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(
                    "ROI:" + roi, label
                ),
                status="WARN",
            )
    if not dict_shifts_available or shift is None:
        # if dictionary of shift or key for this cycle was not found, then it will recalculate XY shift
        images_2d = [np.sum(x, axis=0) for x in images]

        print_log("> Calculating XY shift...")
        shift, _, _ = phase_cross_correlation(
            images_2d[0], images_2d[1], upsample_factor=p["upsample_factor"]
        )

    # applies XY shift to 3D stack
    # ----------------------------
    print_log("$ shifts XY = {}".format(shift))

    # reinterpolate second file in XY using dictionnary to get rough alignment
    images.append(
        apply_xy_shift_3d_images(
            image_3d, shift, parallel_execution=inner_parallel_loop
        )
    )

    del images[1], image_3d  # removes unshifted image to save memory

    # 3D image alignment by block
    # ---------------------------
    print_log("> Block-aligning images in 3D...")
    shift_matrices, block_ref, block_target = image_block_alignment_3d(
        images, block_size_xy=p["blockSizeXY"], upsample_factor=p["upsample_factor"]
    )
    del images  # deletes image list to save memory

    # [plots shift matrices]
    fig2 = plot_3d_shift_matrices(shift_matrices, fontsize=8)

    # combines blocks into a single matrix for display instead of plotting a matrix of subplots each with a block
    outputs = []
    for axis in p["axes2Plot"]:
        outputs.append(
            combine_blocks_image_by_reprojection(
                block_ref, block_target, shift_matrices=shift_matrices, axis1=axis
            )
        )

    ssim_matrices = [x[1] for x in outputs]
    mse_matrices = [x[2] for x in outputs]
    nrmse_matrices = [x[3] for x in outputs]

    fig3 = plt.figure(constrained_layout=False)
    fig3.set_size_inches((20 * 2, 20))
    grid_spec = fig3.add_gridspec(2, 2)
    ax = [
        fig3.add_subplot(grid_spec[:, 0]),
        fig3.add_subplot(grid_spec[0, 1]),
        fig3.add_subplot(grid_spec[1, 1]),
    ]

    titles = ["Z-projection", "X-projection", "Y-projection"]

    for axis, output, i in zip(ax, outputs, p["axes2Plot"]):
        axis.imshow(output[0])
        axis.set_title(titles[i])

    fig3.tight_layout()

    fig4 = plot_3d_shift_matrices(
        ssim_matrices, fontsize=6, log=False, valfmt="{x:.2f}"
    )
    fig4.suptitle("SSIM block matrices")

    fig5 = plot_3d_shift_matrices(mse_matrices, fontsize=6, log=False, valfmt="{x:.2f}")
    fig5.suptitle("mean square root block matrices")

    fig6 = plot_3d_shift_matrices(
        nrmse_matrices, fontsize=6, log=False, valfmt="{x:.2f}"
    )
    fig6.suptitle("normalized root mean square block matrices")

    # saves figures
    # -------------
    fig_titles = [
        "_bkgSubstracted.png",
        "_shiftMatrices.png",
        "_3Dalignments.png",
        "_SSIMblocks.png",
        "_MSEblocks.png",
        "_NRMSEblocks.png",
    ]
    output_filenames = [
        output_folder + os.sep + os.path.basename(filename_to_process) + x
        for x in fig_titles
    ]

    figs = [fig1, fig2, fig3, fig4, fig5, fig6]
    for fig, file in zip(figs, output_filenames):
        fig.savefig(file)
        plt.close(fig)

    # Saves results
    # -------------
    # dict with shift_matrix and NRMSEmatrix: https://en.wikipedia.org/wiki/Root-mean-square_deviation
    # These matrices can be used to apply and assess zxy corrections for any pixel in the 3D image
    # reference file,aligned file,ROI,label,block_i,block_j,shift_z,shift_x,shift_y,quality_xy,quality_zy,quality_zx
    num_blocks, block_xy = block_ref.shape[0], block_ref.shape[-1]
    for i in range(num_blocks):
        for j in range(num_blocks):
            table_entry = [
                os.path.basename(p["fileNameReference"]),
                os.path.basename(filename_to_process),
                int(block_xy),
                int(roi),
                label,
                i,
                j,
                shift_matrices[0][i, j],
                shift_matrices[1][i, j],
                shift_matrices[2][i, j],
                nrmse_matrices[0][i, j],
                nrmse_matrices[1][i, j],
                nrmse_matrices[2][i, j],
            ]
            alignment_results_table.add_row(table_entry)

    print_log("Erasing {} variables\n".format(len(dir()) - 1))
    for var in dir():
        if var != "alignment_results_table":
            del var

    print_log("Variables still alive: {}".format(dir()))

    return alignment_results_table


def create_output_table():
    return Table(
        names=(
            "reference file",
            "aligned file",
            "blockXY",
            "ROI #",
            "label",
            "block_i",
            "block_j",
            "shift_z",
            "shift_x",
            "shift_y",
            "quality_xy",
            "quality_zy",
            "quality_zx",
        ),
        dtype=(
            "S2",
            "S2",
            "int",
            "int",
            "S2",
            "int",
            "int",
            "f4",
            "f4",
            "f4",
            "f4",
            "f4",
            "f4",
        ),
    )
