#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Main file of pyHiM, include the top-level mechanism."""

__version__ = "0.7.2"

import os
import sys

# to remove in a future version
import warnings
from datetime import datetime

import apifish

import core.function_caller as fc
from core.parameters import Parameters
from core.pyhim_logging import print_log
from core.run_args import RunArgs

warnings.filterwarnings("ignore")

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

    him = fc.HiMFunctionCaller(run_args, session_name="HiM_analysis")
    him.initialize()

    him.lauch_dask_scheduler(threads_requested=run_args.thread_nbr, maximum_load=0.8)
    global_param = Parameters(root_folder=run_args.data_path, file_name="infoList.json")
    labels = global_param.param_dict["labels"]

    print_log(f"$ Started logging to: {him.log_file}")
    print_log(f"$ labels to process: {labels}\n")

    for label in labels:
        # sets parameters
        current_param = Parameters(
            root_folder=run_args.data_path,
            label=label,
            file_name="infoList.json",
            stardist_basename=run_args.stardist_basename,
        )

        print_log(
            "--------------------------------------------------------------------------"
        )
        print_log(
            f">                  Analyzing label: {current_param.param_dict['acquisition']['label']}           "
        )
        print_log(
            "--------------------------------------------------------------------------"
        )

        current_param.param_dict["parallel"] = him.parallel
        current_param.param_dict["fileNameMD"] = him.markdown_filename

        # [projects 3D images in 2d]
        if "makeProjections" in run_args.cmd_list:
            him.make_projections(current_param)

        # [registers fiducials using a barcode as reference]
        if "alignImages" in run_args.cmd_list:
            him.align_images(current_param, label)

        # [applies registration to DAPI and barcodes]
        if "appliesRegistrations" in run_args.cmd_list:
            him.apply_registrations(current_param, label)

        # [aligns fiducials in 3D]
        if "alignImages3D" in run_args.cmd_list:
            him.align_images_3d(current_param, label)

        # [segments DAPI and sources in 2D]
        if "segmentMasks" in run_args.cmd_list:
            him.segment_masks(current_param, label)

        # [segments masks in 3D]
        if "segmentMasks3D" in run_args.cmd_list:
            him.segment_masks_3d(current_param, label)

        # [segments sources in 3D]
        if "segmentSources3D" in run_args.cmd_list:
            him.segment_sources_3d(current_param, label)

        # [filters barcode localization table]
        if "filter_localizations" in run_args.cmd_list:
            fc.filter_localizations(current_param, label)

        # [registers barcode localization table]
        if "register_localizations" in run_args.cmd_list:
            fc.register_localizations(current_param, label)

        # [build traces]
        if "build_traces" in run_args.cmd_list:
            fc.build_traces(current_param, label)

        # [builds matrices]
        if "build_matrix" in run_args.cmd_list:
            fc.build_matrix(current_param, label)

        # [builds PWD matrix for all folders with images]
        if "buildHiMmatrix" in run_args.cmd_list:
            him.process_pwd_matrices(current_param, label)

        print("\n")
        del current_param

    # exits
    him.current_session.save()
    print_log("\n==================== Normal termination ====================\n")

    if run_args.parallel:
        him.cluster.close()
        him.client.close()

    del him

    print_log(f"Elapsed time: {datetime.now() - begin_time}")


if __name__ == "__main__":
    if apifish.__version__ < "0.6.4dev":
        sys.exit("ERROR: Please update apifish (git checkout development && git pull)")
    else:
        main()
