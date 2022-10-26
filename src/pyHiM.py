#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 09:11:01 2020

@author: marcnol

This file contains routines to process Hi-M datasets

"""
# =============================================================================
# IMPORTS
# =============================================================================

import os
# to remove in a future version
import warnings
from datetime import datetime

from fileProcessing.fileManagement import Parameters, print_log
import fileProcessing.functionCaller as fc

warnings.filterwarnings("ignore")

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"


# =============================================================================
# MAIN
# =============================================================================

def main():
    begin_time = datetime.now()

    run_parameters = fc.him_parse_arguments()

    him = fc.HiMFunctionCaller(run_parameters, session_name="HiM_analysis")
    him.initialize()

    him.lauch_dask_scheduler(
        threads_requested=run_parameters["threads"], maximum_load=0.8
    )
    global_param = Parameters(
        root_folder=run_parameters["rootFolder"], file_name="infoList.json"
    )
    labels = global_param.param_dict["labels"]

    print_log(f"$ Started logging to: {him.log_file}")
    print_log(f"$ labels to process: {labels}\n")

    for label in labels:

        # sets parameters
        current_param = Parameters(
            root_folder=run_parameters["rootFolder"],
            label=label,
            file_name="infoList.json",
        )

        print_log("--------------------------------------------------------------------------")
        print_log(f">                  Analyzing label: {current_param.param_dict['acquisition']['label']}           ")
        print_log("--------------------------------------------------------------------------")

        current_param.param_dict["parallel"] = him.parallel
        current_param.param_dict["fileNameMD"] = him.markdown_filename

        # [projects 3D images in 2d]
        if "makeProjections" in run_parameters["cmd"]:
            him.make_projections(current_param)

        # [registers fiducials using a barcode as reference]
        if "alignImages" in run_parameters["cmd"]:
            him.align_images(current_param, label)

        # [applies registration to DAPI and barcodes]
        if "appliesRegistrations" in run_parameters["cmd"]:
            him.apply_registrations(current_param, label)

        # [aligns fiducials in 3D]
        if "alignImages3D" in run_parameters["cmd"]:
            him.align_images_3d(current_param, label)

        # [segments DAPI and sources in 2D]
        if "segmentMasks" in run_parameters["cmd"]:
            him.segment_masks(current_param, label)

        # [segments masks in 3D]
        if "segmentMasks3D" in run_parameters["cmd"]:
            him.segment_masks_3d(current_param, label)

        # [segments sources in 3D]
        if "segmentSources3D" in run_parameters["cmd"]:
            him.segment_sources_3d(current_param, label)

        # [filters barcode localization table]
        if "filter_localizations" in run_parameters["cmd"]:
            fc.filter_localizations(current_param, label)

        # [registers barcode localization table]
        if "register_localizations" in run_parameters["cmd"]:
            fc.register_localizations(current_param, label)

        # [build traces]
        if "build_traces" in run_parameters["cmd"]:
            fc.build_traces(current_param, label)

        # [builds matrices]
        if "build_matrix" in run_parameters["cmd"]:
            fc.build_matrix(current_param, label)

        # [builds PWD matrix for all folders with images]
        if "buildHiMmatrix" in run_parameters["cmd"]:
            him.process_pwd_matrices(current_param, label)

        print("\n")
        del current_param

    # exits
    him.current_session.save()
    print_log("\n==================== Normal termination ====================\n")

    if run_parameters["parallel"]:
        him.cluster.close()
        him.client.close()

    del him

    print_log(f"Elapsed time: {datetime.now() - begin_time}")


if __name__ == "__main__":
    main()