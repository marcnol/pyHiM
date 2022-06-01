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

from datetime import datetime

# to remove in a future version
import warnings
warnings.filterwarnings("ignore")

from fileProcessing.fileManagement import Parameters, print_log
from fileProcessing.functionCaller import HiMFunctionCaller, him_parse_arguments

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    run_parameters=him_parse_arguments()
    run_parameters['root_folder'] = "/home/marcnol/data/Embryo_debug_dataset/test_dataset"
    run_parameters['cmd'] = ['alignImages3D']

    him = HiMFunctionCaller(run_parameters, session_name="HiM_analysis")
    him.initialize()
    current_session, current_log = him.current_session, him.current_log

    him.lauch_dask_scheduler(threads_requested=run_parameters["threads"], maximum_load=0.8)
    current_param = Parameters(root_folder=run_parameters["rootFolder"], file_name='infoList.json')
    labels = current_param.param_dict['labels']

    print_log('$ Started logging to: {}'.format(him.log_file))
    print_log("$ labels to process: {}\n".format(labels))

    for label in labels:#range(len(him.labels_to_process)):

        # sets parameters
        current_param = Parameters(root_folder=run_parameters["rootFolder"], label=label, file_name='infoList.json')

        print_log("--------------------------------------------------------------------------")
        print_log(">                  Analyzing label: {}           ".format(current_param.param_dict["acquisition"]["label"]))
        print_log("--------------------------------------------------------------------------")

        current_param.param_dict['parallel']=him.parallel
        current_param.param_dict['markdown_filename']=him.markdown_filename

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

        # [segments sources in 3D]
        if "segmentSources3D" in run_parameters["cmd"]:
            him.segment_sources_3d(current_param, label)

        # [builds PWD matrix for all folders with images]
        if "buildHiMmatrix" in run_parameters["cmd"]:
            him.process_pwd_matrices(current_param, label)

        print("\n")
        del current_param

    # exits
    him.current_session.save()
    print_log("\n===================={}====================\n".format("Normal termination"))

    if run_parameters["parallel"]:
        him.cluster.close()
        him.client.close()

    del him

    print_log("Elapsed time: {}".format(datetime.now() - begin_time))
