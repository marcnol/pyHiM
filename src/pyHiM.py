#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Main file of pyHiM, include the top-level mechanism."""

__version__ = "0.8.8"

import os
import sys
from datetime import datetime

import apifish
import dask.distributed

import core.function_caller as fc
from core.data_manager import DataManager
from core.parameters import Parameters
from core.pyhim_logging import Logger, print_analyzing_label, print_log
from core.run_args import RunArgs

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"


def main(command_line_arguments=None):
    """Main function of pyHiM

    Parameters
    ----------
    command_line_arguments : List[str], optional
        Used to give inputs for the runtime when you call this function like a module.
        For example, to test the pyHiM run from tests folder.
        By default None.
    """
    begin_time = datetime.now()

    run_args = RunArgs(command_line_arguments)

    logger = Logger(
        run_args.data_path,
        parallel=run_args.parallel,
        session_name="HiM_analysis",
        init_msg=run_args.args_to_str(),
    )

    datam = DataManager(
        run_args.data_path,
        logger.md_filename,
        stardist_basename=run_args.stardist_basename,
        param_file=run_args.params_path,
    )

    pipe = fc.Pipeline(datam, run_args.cmd_list, run_args.parallel, logger)
    pipe.lauch_dask_scheduler(threads_requested=run_args.thread_nbr, maximum_load=0.8)

    pipe.run()

    # Separate no-refactor routines
    #      ||
    #      ||
    #      ||
    #      ||
    #      ||
    #      ||
    #      ||
    #      ||
    # No-refactored
    # \            /
    #  \          /
    #   \        /
    #    \      /
    #     \    /
    #      \  /
    #       \/

    print_log("\n\n\n")
    raw_dict = datam.load_user_param()
    global_param = Parameters(raw_dict, root_folder=datam.m_data_path)
    labels = datam.get_processable_labels()
    print_log(f"$ Labels to process: {labels}")
    for label in labels:
        # sets parameters with old way (temporary during pyHiM restructuration)
        current_param = Parameters(
            raw_dict,
            root_folder=datam.m_data_path,
            label=label,
            stardist_basename=datam.m_stardist_basename,
        )

        print_analyzing_label(
            f"Analyzing label: {current_param.param_dict['acquisition']['label']}"
        )

        current_param.param_dict["parallel"] = pipe.parallel
        current_param.param_dict["fileNameMD"] = logger.md_filename

        # [applies registration to DAPI and barcodes]
        if "register_global" in pipe.cmds:
            registration_params = datam.labelled_params[label].registration
            pipe.apply_registrations(
                current_param, label, datam.m_data_path, registration_params
            )

        # [aligns fiducials in 3D]
        if "register_local" in pipe.cmds:
            registration_params = datam.labelled_params[label].registration
            pipe.align_images_3d(
                current_param,
                label,
                datam.m_data_path,
                registration_params,
                datam.dict_shifts_path,
            )

        # [segments DAPI and sources in 2D]
        if "mask_2d" in pipe.cmds or "localize_2d" in pipe.cmds:
            segmentation_params = datam.labelled_params[label].segmentation
            pipe.segment_masks(
                current_param,
                label,
                datam.m_data_path,
                segmentation_params,
                datam.align_folder,
            )

        # [segments masks in 3D]
        if "mask_3d" in pipe.cmds and (label in ("DAPI", "mask")):
            segmentation_params = datam.labelled_params[label].segmentation
            pipe.segment_masks_3d(
                current_param,
                label,
                datam.processed_roi,
                datam.m_data_path,
                segmentation_params,
                datam.dict_shifts_path,
            )

        # [segments sources in 3D]
        if "localize_3d" in pipe.cmds and label == "barcode":
            segmentation_params = datam.labelled_params[label].segmentation
            pipe.segment_sources_3d(
                current_param,
                label,
                datam.processed_roi,
                datam.m_data_path,
                segmentation_params,
                datam.dict_shifts_path,
            )

        # [filters barcode localization table]
        if "filter_localizations" in pipe.cmds and label == "barcode":
            segmentation_params = datam.labelled_params[label].segmentation
            fc.filter_localizations(
                current_param, label, datam.m_data_path, segmentation_params
            )

        # [registers barcode localization table]
        if "register_localizations" in pipe.cmds and label == "barcode":
            segmentation_params = datam.labelled_params[label].segmentation
            fc.register_localizations(
                current_param,
                label,
                datam.m_data_path,
                datam.local_shifts_path,
                segmentation_params,
            )

        # [build traces]
        if "build_traces" in pipe.cmds and label == "barcode":
            segmentation_params = datam.labelled_params[label].segmentation
            matrix_params = datam.labelled_params[label].matrix
            fc.build_traces(
                current_param,
                label,
                datam.m_data_path,
                segmentation_params,
                matrix_params,
            )

        # [builds matrices]
        if "build_matrix" in pipe.cmds and label == "barcode":
            matrix_params = datam.labelled_params[label].matrix
            fc.build_matrix(current_param, label, datam.m_data_path, matrix_params)

        print_log("\n")
        del current_param

    # exits
    if pipe.parallel:
        pipe.m_dask.cluster.close()
        pipe.m_dask.client.close()

    del pipe

    print_log("\n==================== Normal termination ====================\n")
    print_log(f"Elapsed time: {datetime.now() - begin_time}")


def check_version_compatibily():
    if apifish.__version__ < "0.6.4dev":
        sys.exit("ERROR: Please update apifish (git checkout development && git pull)")
    if dask.distributed.__version__ < "2023.4.1":
        sys.exit(
            "ERROR: dask[distributed] version: deprecated. \nPlease update dask[distributed] \
                (pip install -U distributed)"
        )


if __name__ == "__main__":
    print(f"[VERSION] pyHiM {__version__}")
    check_version_compatibily()
    main()
