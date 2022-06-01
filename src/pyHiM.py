#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 09:11:01 2020

@author: marcnol

This file contains routines to process Hi-M datasets

"""
# =============================================================================
# IMPORTS
# =============================================================================q

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


def old_main(run_parameters):

    him = fc.HiMFunctionCaller(run_parameters, session_name="HiM_analysis")
    him.initialize()

    him.lauch_dask_scheduler(
        threads_requested=run_parameters["threads"], maximum_load=0.8
    )
    global_param = Parameters(
        root_folder=run_parameters["root_folder"], file_name="infoList.json"
    )
    labels = global_param.param_dict["labels"]

    print_log("$ Started logging to: {}".format(him.log_file))
    print_log("$ labels to process: {}\n".format(labels))

    for label in labels:

        # sets parameters
        current_param = Parameters(
            root_folder=run_parameters["root_folder"],
            label=label,
            file_name="infoList.json",
        )

        print_log(
            "--------------------------------------------------------------------------"
        )
        print_log(
            ">                  Analyzing label: {}           ".format(
                current_param.param_dict["acquisition"]["label"]
            )
        )
        print_log(
            "--------------------------------------------------------------------------"
        )

        current_param.param_dict["parallel"] = him.parallel
        current_param.param_dict["markdown_filename"] = him.markdown_filename

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
    print_log(
        "\n===================={}====================\n".format("Normal termination")
    )

    if run_parameters["parallel"]:
        him.cluster.close()
        him.client.close()

    del him


def main(command_line_arguments=None):
    """main function to run pyHiM

    Parameters
    ----------
    command_line_arguments : List, optional
        Used for test functions, by default None
    """

    # run_args = RunArgs(command_line_arguments)

    # pipeline = Pipeline(run_args.cmd_list)

    # data = InputData(run_args.data_path)
    # data.check_consistency(pipeline)

    # # TODO
    # parameters = InputParameters(run_args.param_path)
    # parameters.check_consistency(pipeline, data)

    # pipeline.build_output_folders(run_args.output_path)
    # pipeline.fill(data, parameters)
    # pipeline.run()


if __name__ == "__main__":
    begin_time = datetime.now()

    run_parameters = fc.him_parse_arguments()

    if "new" in run_parameters["cmd"]:
        main()
    else:
        old_main(run_parameters)

    print_log("Elapsed time: {}".format(datetime.now() - begin_time))
