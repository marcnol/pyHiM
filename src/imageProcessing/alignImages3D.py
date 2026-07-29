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
from datetime import datetime

import gc
import cupy as cp
import tifffile as tiff

import matplotlib.pylab as plt
import numpy as np
from astropy.table import Table, vstack
from skimage import io
from skimage import exposure
from numpy.typing import ArrayLike
from scipy.ndimage import zoom

#import warpfield 3D registartions 
from warpfield.warp import warp_volume
from warpfield import Recipe, register_volumes
from warpfield.register import WarpMap

# from skimage.registration import phase_cross_correlation
from core.dask_cluster import try_get_client
from core.parameters import RegistrationParams, load_alignment_dict, print_dict
from core.pyhim_logging import print_log, print_session_name
from core.saving import plot_3d_shift_matrices, plot_4_images
from imageProcessing.alignImages import (
    apply_xy_shift_3d_images,
    combine_blocks_image_by_reprojection,
    image_block_alignment_3d,
)
from imageProcessing.imageProcessing import preprocess_3d_image
from imageProcessing.makeProjections import reinterpolate_z

from pyHiM_tools import BothImgRbgFile


# =============================================================================
# CLASSES
# =============================================================================


class Drift3D:
    def __init__(self, param, registration_params, parallel=False):
        self.current_param = param
        self.reg_params = registration_params
        self.window = 3
        self.parallel = parallel
        self.p = {}
        self.image_ref_0 = None
        self.image_ref = None
        self.filenames_to_process_list = []
        self.filenames_with_ref_barcode = ""
        self.dict_shifts = None
        self.dict_shifts_available = None
        self.inner_parallel_loop = None

    def align_fiducials_3d_file(
        self,
        filename_to_process,
        data_path,
        params: RegistrationParams,
        roi_name,
        cycle_name,
        z_binning,
    ):
        """
        Aligns <filename_to_process> fiducial against reference

        Returns
        -------
        None.

        """

        p = self.p
        alignment_results_table = create_output_table()

        # excludes the reference fiducial and processes files in the same ROI
        inner_parallel_loop = self.inner_parallel_loop
        image_ref = self.image_ref
        image_ref_0 = self.image_ref_0
        dict_shifts_available = self.dict_shifts_available
        dict_shifts = self.dict_shifts
        output_folder = data_path + os.sep + params.register_local_folder

        return _align_fiducials_3d_file(
            filename_to_process,
            alignment_results_table,
            p,
            roi_name,
            cycle_name,
            inner_parallel_loop,
            image_ref,
            image_ref_0,
            dict_shifts,
            dict_shifts_available,
            output_folder,
            params,
            z_binning,
        )

    def load_reference_fiducial(
        self, filename_reference, z_binning, lower_threshold_3d, higher_threshold_3d
    ):
        """
        Loads Reference fiducial image

        Returns
        -------
        None.

        """

        self.p["fileNameReference"] = filename_reference
        print_log(f"Loading reference 3D image: {filename_reference}")

        self.image_ref_0, self.image_ref = load_n_preprocess_image(
            filename_reference,
            z_binning,
            lower_threshold_3d,
            higher_threshold_3d,
            parallel_execution=False,
        )

        self.image_ref_0 = np.sum(
            self.image_ref_0, axis=0
        )  # replaces 3D with a 2D projection

        print_log(f"$ Found {len(self.filenames_to_process_list)} files.")

    def load_dict_shifts(self, dict_shifts_path):
        """
        Lods dictionary of XY shifts

        Returns
        -------
        None.

        """
        print_log(f"""\nReference barcode: {self.reg_params.referenceFiducial}""")

        for file in self.current_param.files_to_process:
            if self.reg_params.referenceFiducial in file.split("_"):
                self.filenames_with_ref_barcode = file

        # loads dicShifts with shifts for all rois and all labels
        self.dict_shifts, self.dict_shifts_available = load_alignment_dict(
            dict_shifts_path
        )

    def align_fiducials_3d_in_folder(
        self,
        data_path,
        dict_shifts_path,
        params: RegistrationParams,
        roi_name,
        z_binning,
        single_file_to_process=None,
    ):
        """
        Refits all the barcode files found in root_folder

        Returns
        -------
        None.

        """
        now = datetime.now()
        print_dict(self.p)

        # gets files to process
        files_folder = glob.glob(data_path + os.sep + "*.tif")
        self.current_param.find_files_to_process(files_folder)

        # loads dictionary of shifts
        self.load_dict_shifts(dict_shifts_path)

        # creates Table that will hold results
        alignment_results_table_global = create_output_table()
        alignment_results_tables = []

        client = try_get_client()

        # loads reference fiducial image for this ROI
        self.load_reference_fiducial(
            self.filenames_with_ref_barcode,
            z_binning,
            params._3D_lower_threshold,
            params._3D_higher_threshold,
        )
        self.filenames_to_process_list = [
            x
            for x in self.current_param.files_to_process
            if (x != self.filenames_with_ref_barcode)
        ]
        if single_file_to_process:
            self.filenames_to_process_list = [
                x
                for x in self.filenames_to_process_list
                if os.path.basename(x) == os.path.basename(single_file_to_process)
            ]
            if not self.filenames_to_process_list:
                raise SystemExit(
                    f"Requested file '{single_file_to_process}' was not found among the files to process."
                )
        number_files = len(self.filenames_to_process_list)

        if client is None:
            self.inner_parallel_loop = True
            for file_index, filename_to_process in enumerate(
                self.filenames_to_process_list
            ):
                print_log(f"\n\n>>>Iteration: {file_index+1}/{number_files}<<<")

                alignment_results_tables.append(
                    self.align_fiducials_3d_file(
                        filename_to_process,
                        data_path,
                        params,
                        roi_name,
                        find_cycle(self.current_param, filename_to_process),
                        z_binning,
                    )
                )

        else:
            self.inner_parallel_loop = False
            nb_workers = len(client.scheduler_info()["workers"])
            print_log(f"> Aligning {number_files} files using {nb_workers} workers...")

            futures = [
                client.submit(
                    self.align_fiducials_3d_file,
                    x,
                    data_path,
                    params,
                    roi_name,
                    find_cycle(self.current_param, x),
                    z_binning,
                )
                for x in self.filenames_to_process_list
            ]

            alignment_results_tables = client.gather(futures)
            print_log(
                f"> Retrieving {len(alignment_results_tables)} results from cluster"
            )

            # del futures

        output_prefix = params.outputFile
        if single_file_to_process:
            output_prefix = (
                output_prefix
                + "_"
                + os.path.splitext(os.path.basename(single_file_to_process))[0]
            )

        if single_file_to_process and alignment_results_tables:
            alignment_results_table_global = alignment_results_tables[0]
        else:
            alignment_results_table_global = vstack(
                [alignment_results_table_global] + alignment_results_tables
            )

        path_name = (
            data_path
            + os.sep
            + params.register_local_folder
            + os.sep
            + "data"
            + os.sep
            + output_prefix
        )
        local_shifts_path = path_name + "_block3D.ecsv"

        alignment_results_table_global.write(
            local_shifts_path,
            format="ascii.ecsv",
            overwrite=True,
        )

        print_log(f"$ register_local processing time: {datetime.now() - now}")
        print_log(f"$ register_local output Table saved in: {local_shifts_path}")

        return local_shifts_path

    def align_fiducials_3d(
        self,
        data_path,
        params: RegistrationParams,
        dict_shifts_path,
        roi_name,
        z_binning,
        single_file_to_process=None,
    ):
        """
        runs refitting routine in root_folder

        Returns
        -------
        None.

        """
        session_name = "register_local"

        # processes folders and files
        print_session_name(session_name)

        print_log(f"-------> Processing Folder: {data_path}")
        # self.current_log.parallel = self.parallel

        local_shifts_path = self.align_fiducials_3d_in_folder(
            data_path,
            dict_shifts_path,
            params,
            roi_name,
            z_binning,
            single_file_to_process=single_file_to_process,
        )

        print_log(f"HiM matrix in {data_path} processed")

        return local_shifts_path

    def WarpfieldRegistration(self,
        data_path,
        reference,
        tomove,
    ):
        """
        runs warpfield registration for one image 

        Returns
        -------
        None.

        """
        session_name = ""
        zbin = 2
        xybin = 2
        gpu = 0
        # processes folders and files
        print_session_name(session_name)

        print_log(f"-------> Processing Folder: {data_path}")
        # self.current_log.parallel = self.parallel
        
        moving = tiff.imread(self)
        moving_dtype = moving.dtype
        original_shape=moving.shape
        reference = tiff.imread(reference)
        tomove_dtype = None
        tomove_image = None
        
        if tomove is not None:
            tomove_image = tiff.imread(tomove)
            tomove_dtype = tomove_image.dtype
    
        # binning
        if xybin > 1 or zbin > 1 :
            moving = zoom(moving, (1.0 / zbin, 1.0 / xybin, 1.0 / xybin), order=1)
            reference = zoom(reference, (1.0 / zbin,  1.0 / xybin,  1.0 / xybin), order=1)
            if tomove is not None :
                tomove_image = zoom(tomove_image, ( 1.0 / zbin, 1.0 / xybin, 1.0 / xybin), order=1)
                
        # save warpfield
        base = os.path.splitext(os.path.basename(moving_path))[0]
        h5_path = os.path.join(output, f"{base}_warp_map.h5")
        moving_registered, warp_field, tomove_registered = compute_warpfield(
            reference,
            moving,
            tomove_image,
            h5_path,
            gpu_id=gpu )
        
        # RGB overlay
        os.makedirs(output, exist_ok=True)
        overlay = BothImgRbgFile(reference.max(axis=0), moving.max(axis=0), tag='reference_original')
        overlay.save(output, f"{base}_registered")
        overlay = BothImgRbgFile(reference.max(axis=0), moving_registered.max(axis=0), tag='reference_aligned')
        overlay.save(output,  f"{base}_registered")
        
        # Plot the intensity and direction of the deformation field at the center z-plane
        z_plane = warp_field.shape[1] // 2 # (3,z,x,y)
        plot_deformation_intensity_xyz(warp_field, z_plane, f"{base}")
        plot_deformation_direction(warp_field, z_plane, f"{base}")
        
        # Upsample back to the original shape if binning was applied
        if zbin > 1 or xybin > 1 :
            zoom_factors = [original_shape[0] / moving_registered.shape[0],  # Z upsampling
                            original_shape[1] / moving_registered.shape[1],  # Y upsampling
                            original_shape[2] / moving_registered.shape[2]]  # X upsampling
            if tomove_registered is not None :
                tomove_registered = zoom(tomove_registered, zoom_factors, order=1)
                
        # Restore tomove image dtype
        if tomove_registered is not None and tomove_dtype is not None:
            if np.issubdtype(tomove_dtype, np.integer):
                info = np.iinfo(tomove_dtype)
                tomove_registered = np.clip(tomove_registered,info.min, info.max).astype(tomove_dtype)
            else:
                tomove_registered = tomove_registered.astype(tomove_dtype) 
                
       return tomove_registered 

# =============================================================================
#   FUNCTIONS
# =============================================================================


def find_cycle(param, filename):
    return str(param.decode_file_parts(os.path.basename(filename))["cycle"])


def load_n_preprocess_image(
    filename_to_process,
    z_binning,
    lower_threshold,
    higher_threshold,
    parallel_execution=True,
):
    print_log(f"$ File:{os.path.basename(filename_to_process)}")

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


def _format_xy_alignment_axis(axis):
    axis.set_xlabel("X pixel")
    axis.set_ylabel("Y pixel")
    axis.tick_params(axis="both", which="both", labelsize=8)


def _format_slice_alignment_axis(
    axis,
    image,
    slice_positions,
    number_z_planes,
    slice_axis_label,
    horizontal_axis_label,
):
    axis.set_xlabel(f"{horizontal_axis_label} pixel")
    axis.set_ylabel("Z slice (per montage row)")
    if slice_positions is None or len(slice_positions) == 0:
        return

    slice_height = number_z_planes
    separator_height = 1
    row_stride = slice_height + separator_height
    tick_positions = [
        i * row_stride + (slice_height - 1) / 2 for i in range(len(slice_positions))
    ]
    z_last = number_z_planes - 1
    axis.set_yticks(tick_positions)
    axis.set_yticklabels(
        [f"{slice_axis_label}={position}\nZ=0-{z_last}" for position in slice_positions]
    )
    axis.tick_params(axis="both", which="both", labelsize=8)

    sampled_positions = ", ".join(str(position) for position in slice_positions)
    axis.text(
        0.99,
        0.02,
        f"sampled {slice_axis_label}: {sampled_positions}",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none"},
    )
    axis.set_ylim(-0.5, image.shape[0] - 0.5)


def _align_fiducials_3d_file(
    filename_to_process,
    alignment_results_table,
    p,
    roi,
    cycle_name,
    inner_parallel_loop,
    image_ref,
    image_ref_0,
    dict_shifts,
    dict_shifts_available,
    output_folder,
    params: RegistrationParams,
    z_binning,
):
    # - load  and preprocesses 3D fiducial file
    print_log(f"\n\n>>>Processing roi:[{roi}] cycle:[{cycle_name}]<<<")
    image_3d_0, image_3d = load_n_preprocess_image(
        filename_to_process,
        z_binning,
        params._3D_lower_threshold,
        params._3D_higher_threshold,
        parallel_execution=inner_parallel_loop,
    )

    # shows original images and background subtracted
    image_3d_0 = np.sum(image_3d_0, axis=0)  # replaces by a 2D projection
    images = [image_ref, image_3d]
    images_2d = [np.sum(x, axis=0) for x in images]
    fig1 = plot_4_images(
        [image_ref_0, image_3d_0] + images_2d,
        titles=["reference", "cycle <i>", "processed reference", "processed cycle <i>"],
    )

    del image_3d_0

    # gets shift values from dictionary
    # ---------------------------------
    if dict_shifts_available:
        # uses existing shift calculated by align_images
        try:
            shift = dict_shifts["ROI:" + roi][cycle_name]
        except KeyError:
            shift = None
            print_log(
                f"Could not find dictionary with alignment parameters for this ROI: ROI:{roi}, cycle: {cycle_name}",
                status="WARN",
            )
    if not dict_shifts_available or shift is None:
        # if dictionary of shift or key for this cycle was not found, then it will exit

        raise SystemExit(
            f"> Existing with ERROR: Could not find shift value \
                for this ROI: {roi} and cycle: {cycle_name}"
        )

    # applies XY shift to 3D stack
    # ----------------------------
    print_log(f"$ shift values that will be applied = {shift}")

    # reinterpolate second file in XY or XYZ using dictionary to get rough alignment
    images.append(
        apply_xy_shift_3d_images(
            image_3d, shift, parallel_execution=inner_parallel_loop
        )
    )

    del images[1], image_3d  # removes unshifted image to save memory

    # Refines 3D image alignment by block decomposition
    # -------------------------------------------------
    print_log("> Block-aligning images in 3D...")
    shift_matrices, block_ref, block_target = image_block_alignment_3d(
        images, block_size_xy=params.blockSizeXY, upsample_factor=params.upsample_factor
    )
    del images  # deletes image list to save memory

    # [plots shift matrices]
    fig2 = plot_3d_shift_matrices(shift_matrices, fontsize=8)

    # combines blocks into a single matrix for display instead of plotting a matrix
    # of subplots each with a block
    number_blocks_y, number_blocks_x = block_ref.shape[:2]
    output_xy = combine_blocks_image_by_reprojection(
        block_ref, block_target, shift_matrices=shift_matrices, axis1=0
    )
    output_xz = combine_blocks_image_by_reprojection(
        block_ref,
        block_target,
        shift_matrices=shift_matrices,
        axis1=1,
        number_slices=number_blocks_y,
        return_slice_positions=True,
    )
    output_yz = combine_blocks_image_by_reprojection(
        block_ref,
        block_target,
        shift_matrices=shift_matrices,
        axis1=2,
        number_slices=number_blocks_x,
        return_slice_positions=True,
    )
    outputs = [output_xy, output_xz[:4], output_yz[:4]]
    slice_positions = [None, output_xz[4], output_yz[4]]

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

    titles = ["XY Z-projection", "XZ slices across Y", "YZ slices across X"]

    for axis, output, i in zip(ax, outputs, range(3)):
        if i == 0:
            axis.imshow(output[0])
            _format_xy_alignment_axis(axis)
        else:
            axis.imshow(output[0], origin="lower", aspect="auto")
            _format_slice_alignment_axis(
                axis,
                output[0],
                slice_positions[i],
                block_ref.shape[2],
                slice_axis_label="Y" if i == 1 else "X",
                horizontal_axis_label="X" if i == 1 else "Y",
            )
        axis.set_title(titles[i])

    fig3.tight_layout()

    fig5 = plot_3d_shift_matrices(mse_matrices, fontsize=6, log=False, valfmt="{x:.2f}")
    fig5.suptitle("mean square root block matrices")

    # saves figures
    # -------------
    fig_titles = [
        "_bkgSubstracted.png",
        "_shiftMatrices.png",
        "_3Dalignments.png",
        "_MSEblocks.png",
    ]
    output_filenames = [
        output_folder + os.sep + os.path.basename(filename_to_process) + x
        for x in fig_titles
    ]

    figs = [fig1, fig2, fig3, fig5]
    for fig, file in zip(figs, output_filenames):
        fig.savefig(file)
        plt.close(fig)

    # Saves results
    # -------------
    # dict with shift_matrix and NRMSEmatrix: https://en.wikipedia.org/wiki/Root-mean-square_deviation
    # These matrices can be used to apply and assess zxy corrections for any pixel in the 3D image
    # reference file,aligned file,ROI,cycle,block_i,block_j,shift_z,shift_x,shift_y,quality_xy,quality_zy,quality_zx
    num_blocks, block_xy = block_ref.shape[0], block_ref.shape[-1]
    for i in range(num_blocks):
        for j in range(num_blocks):
            table_entry = [
                os.path.basename(p["fileNameReference"]),
                os.path.basename(filename_to_process),
                int(block_xy),
                int(roi),
                cycle_name,
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

    # Erasing variables, TODO: check if it's necessary
    for var in dir():
        if var != "alignment_results_table":
            del var

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

def compute_warpfield(
    img_ref: ArrayLike, 
    img_trg: ArrayLike,
    tomove_image: ArrayLike,
    h5_path: str,
    gpu_id: int = 0
) -> tuple[ArrayLike, ArrayLike, ArrayLike, ArrayLike, ArrayLike | None] :
    """
    Compute the warpfield to warp a target image to a reference image. Applies warp_map to tomoveimage.
    """

    cp.cuda.Device(gpu_id).use()

    recipe = ( Recipe() )  # initialized with a translation level, followed by an affine registration level
    recipe.pre_filter.clip_thresh = 0  # clip DC background, if present
    
    recipe.pre_filter.soft_edge = [4,33,33]

    # affine level properties
    recipe.levels[-1].repeats = 0

    # add non-rigid registration levels:
    recipe.add_level(block_size=[15, 33, 33]) # adjust block_size to make blocks roughly isotropic in real space. 
    recipe.levels[-1].block_stride = 0.75
    recipe.levels[-1].smooth.sigmas = [1.0, 1.0, 1.0] 
    recipe.levels[-1].smooth.long_range_ratio = 0.1
    recipe.levels[-1].repeats = 2
        
    recipe.add_level(block_size=[5,11,11])
    recipe.levels[-1].block_stride = 0.70
    recipe.levels[-1].smooth.sigmas = [0.5,0.5,0.5] 
    recipe.levels[-1].smooth.long_range_ratio = 0.05 # Long range ratio for double gaussian kernel.
    recipe.levels[-1].repeats = 2

    #register moving volume
    warped_image, warp_map, _ = register_volumes(
        ref=img_ref,
        vol=img_trg,
        recipe=recipe )

    # save warpfield as h5 
    warp_map.to_h5(
        h5_path,
        group="warp_map",
        compression="gzip",
        overwrite=True,
        )
    # save as np
    warped_image = cp.asnumpy(warped_image).astype(np.float32)
    warp_field = cp.asnumpy(warp_map.warp_field).astype(np.float32)
    block_size = cp.asnumpy(warp_map.block_size).astype(np.float32)
    block_stride = cp.asnumpy(warp_map.block_stride).astype(np.float32)

    tomove_registered = None
    
    # apply warp to other channel.s
    if tomove_image is not None:
        offset = -(block_size/ block_stride/ 2)
        tomove_registered_cp = warp_volume(
            cp.asarray(tomove_image, dtype=cp.float32),
            cp.asarray(warp_field),
            cp.asarray(block_stride),
            cp.asarray(offset, dtype=cp.float32)
            )
        
        tomove_registered = cp.asnumpy(tomove_registered_cp).astype(np.float32)
        del tomove_registered_cp 

    
    del warp_map
    gc.collect()
    cp.cuda.Stream.null.synchronize()
    cp.get_default_memory_pool().free_all_blocks()
    cp.get_default_pinned_memory_pool().free_all_blocks()

    return (warped_image, warp_field, tomove_registered)

class BothImgRbgFile:
    def __init__(self, image1, image2, tag='', title=''):
        self.image1 = image1
        self.image2 = image2
        self.tag = tag
        if title is None:
            self.title = tag  # gets title from tag
        else:
            self.title = title  # New attribute to hold the title

    def save(self, folder_path, basename):
        self.folder_path = folder_path
        self.basename = f"{basename}_{self.tag}_overlay"
        self.path_name = os.path.join(self.folder_path, self.basename + ".png")
        
        # Normalize images and rescale intensity
        img_1 = self.image1 / self.image1.max()
        img_2 = self.image2 / self.image2.max()
        img_1 = exposure.rescale_intensity(img_1, out_range=(0, 1))
        img_2 = exposure.rescale_intensity(img_2, out_range=(0, 1))
        
        # Create the figure and axis
        fig, ax1 = plt.subplots()
        fig.set_size_inches((30, 30))
        
        # Create RGB overlay image
        null_image = np.zeros(img_1.shape)
        rgb = np.dstack([img_1, img_2, null_image])
        
        # Display the image and set the title
        ax1.imshow(rgb)
        ax1.axis("off")
        ax1.set_title(self.title)  # Set the title of the figure
        
        # Save the figure
        fig.savefig(self.path_name)
        plt.close(fig)

def compute_intensity_np(displacement_field, z_plane):
    dx = displacement_field[0, z_plane, :, :]
    dy = displacement_field[1, z_plane, :, :]
    dz = displacement_field[2, z_plane, :, :]

    intensity = np.sqrt(dx**2 + dy**2 + dz**2)
    return [intensity, dx, dy, dz]

def plot_deformation_intensity(displacement_field, z_plane, output_prefix):
    intensity,_,_,_ = compute_intensity_np(displacement_field, z_plane)
    plt.figure(figsize=(10, 8))
    plt.imshow(intensity, cmap='Reds')
    plt.colorbar(label='Vector Field Intensity')
    plt.title(f'Intensity of Vector Field at Z-plane {z_plane}')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.savefig(f"{output_prefix}_intensity_z{z_plane}.png")
    plt.close()

def plot_deformation_intensity_xyz(displacement_field, z_plane, output_prefix):
    data = compute_intensity_np(displacement_field, z_plane)
    titles = ["magnitude", "dx", "dy", "dz"]

    fig, axes = plt.subplots(2, 2)
    fig.set_size_inches((10, 10))
    ax = axes.ravel()

    for axis, img, title in zip(ax, data, titles):
        vmin, vmax = 0, np.max(img)
        cmap="YlOrRd"
        im = axis.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
        axis.set_title(title)    
        axis.set_xlabel('X-axis')
        axis.set_ylabel('Y-axis')
        cbar1 = fig.colorbar(im, ax=axis, shrink=0.5)
        cbar1.set_label('pixels')

    fig.tight_layout()
    fig.suptitle(f'Intensity of Vector Field at Z-plane {z_plane}')

    fig.savefig(f"{output_prefix}_DF_intensity_z{z_plane}.png")


def compute_direction_np(warp, z_plane):
    dx = warp[0, z_plane, :, :]
    dy = warp[1, z_plane, :, :]
    
    # Compute the direction in the XY plane (arctangent of dy/dx)
    direction = np.arctan2(dy, dx)
    
    # Normalize the direction to the range [0, 1] for color mapping
    norm = plt.Normalize(-np.pi, np.pi)
    direction_normalized = norm(direction)
    
    return direction_normalized

def plot_deformation_direction(displacement_field, z_plane, output_prefix):
    direction = compute_direction_np(displacement_field, z_plane)
    plt.figure(figsize=(10, 8))
    plt.imshow(direction, cmap='twilight', alpha=0.9, norm=plt.Normalize(-np.pi, np.pi))
    cbar = plt.colorbar(ticks=[-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
    cbar.ax.set_yticklabels(['-π', '-π/2', '0', 'π/2', 'π'])
    cbar.set_label('Vector Field Direction (radians)')
    plt.title(f'Direction of Vector Field at Z-plane {z_plane}')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.savefig(f"{output_prefix}_DF_direction_z{z_plane}.png")
    plt.close()


