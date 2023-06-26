#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 15:26:00 2020

@author: marcnol
"""


import argparse
import logging
import os
from datetime import datetime

from core.dask_cluster import DaskCluster
from core.parameters import load_json, print_dict, save_json
from core.pyhim_logging import Log, print_log, write_string_to_file
from imageProcessing.alignImages import align_images, apply_registrations
from imageProcessing.alignImages3D import Drift3D
from imageProcessing.makeProjections import make_projections
from imageProcessing.segmentMasks import segment_masks
from imageProcessing.segmentMasks3D import SegmentMasks3D
from imageProcessing.segmentSources3D import SegmentSources3D
from matrixOperations.alignBarcodesMasks import process_pwd_matrices
from matrixOperations.build_matrix import BuildMatrix
from matrixOperations.build_traces import BuildTraces
from matrixOperations.filter_localizations import FilterLocalizations
from matrixOperations.register_localizations import RegisterLocalizations


class Session:
    def __init__(self, root_folder, name="dummy"):
        now = datetime.now()
        session_root_name = now.strftime("%d%m%Y_%H%M%S")
        self.file_name = root_folder + os.sep + "Session_" + session_root_name + ".json"
        self.name = name
        self.data = {}

    # loads existing session
    def load(self):
        self.data = load_json(self.file_name)
        print(f"Session information read: {self.file_name}")

    # saves session to file
    def save(self):
        save_json(self.file_name, self.data)
        print(f"> Saved json session file to {self.file_name}")

    # add new task to session
    def add(self, key, value):
        if key not in self.data:
            self.data[key] = value
        else:
            self.data[key] = [self.data[key], value]

    def clear_data(self):
        self.data = {}


class HiMFunctionCaller:
    def __init__(self, run_parameters, session_name="HiM_analysis"):
        self.run_parameters = run_parameters
        self.root_folder = run_parameters["rootFolder"]
        self.parallel = run_parameters["parallel"]
        self.session_name = session_name

        self.current_log = Log(root_folder=self.root_folder, parallel=self.parallel)

        self.current_session = Session(self.root_folder, self.session_name)

        self.log_file = ""
        self.markdown_filename = ""
        self.client = None
        self.cluster = None

    def initialize(self):
        print_log(
            "\n--------------------------------------------------------------------------"
        )

        print_log(f"$ root_folder: {self.root_folder}")

        begin_time = datetime.now()

        #####################
        # setup markdown file
        #####################
        print_log(
            f"\n======================{self.session_name}======================\n"
        )
        now = datetime.now()
        date_time = now.strftime("%d%m%Y_%H%M%S")

        filename_root = "HiM_analysis"
        self.log_file = self.root_folder + os.sep + filename_root + date_time + ".log"
        self.markdown_filename = self.log_file.split(".")[0] + ".md"

        print_log(f"$ Hi-M analysis will be written tos: {self.markdown_filename}")
        write_string_to_file(
            self.markdown_filename,
            f"""# Hi-M analysis {begin_time.strftime("%Y/%m/%d %H:%M:%S")}""",
            "w",
        )  # initialises MD file

        ##############
        # setupLogger
        ##############

        # creates output formats for terminal and log file
        formatter1 = logging.Formatter("%(asctime)s: %(levelname)s: %(message)s")
        formatter2 = logging.Formatter("%(message)s")

        # clears up any existing logger
        logger = logging.getLogger()
        logger.handlers = []
        for hdlr in logger.handlers[:]:
            if isinstance(hdlr, logging.FileHandler):
                logger.removeHandler(hdlr)

        # initializes handlers for terminal and file
        filehandler = logging.FileHandler(self.log_file, "w")
        ch = logging.StreamHandler()

        filehandler.setLevel(logging.INFO)
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)

        logger.addHandler(ch)
        logger.addHandler(filehandler)

        filehandler.setFormatter(formatter1)
        ch.setFormatter(formatter2)

    def lauch_dask_scheduler(self, threads_requested=25, maximum_load=0.8):
        if self.parallel:
            print_log(f"$ Requested {threads_requested} threads")

            dask_cluster_instance = DaskCluster(
                threads_requested, maximum_load=maximum_load
            )

            dask_cluster_instance.create_distributed_client()
            self.client = dask_cluster_instance.client
            self.cluster = dask_cluster_instance.cluster

    def make_projections(self, current_param):
        if not self.run_parameters["parallel"]:
            make_projections(current_param, self.current_session)
        else:
            result = self.client.submit(
                make_projections, current_param, self.current_session
            )
            _ = self.client.gather(result)

    def align_images(self, current_param, label):
        if (
            label == "fiducial"
            and current_param.param_dict["acquisition"]["label"] == "fiducial"
        ):
            print_log(f"> Making image registrations for label: {label}")
            if not self.parallel:
                align_images(current_param, self.current_session)
            else:
                result = self.client.submit(
                    align_images, current_param, self.current_session
                )
                _ = self.client.gather(result)

    def align_images_3d(self, current_param, label):
        if (
            label == "fiducial"
            and "block3D" in current_param.param_dict["alignImages"]["localAlignment"]
        ):
            print_log(f"> Making 3D image registrations label: {label}")
            _drift_3d = Drift3D(
                current_param, self.current_session, parallel=self.parallel
            )
            _drift_3d.align_fiducials_3d()

    def apply_registrations(self, current_param, label):
        if (
            label != "fiducial"
            and current_param.param_dict["acquisition"]["label"] != "fiducial"
        ):
            print_log(f"> Applying image registrations for label: {label}")

            if not self.parallel:
                apply_registrations(current_param, self.current_session)
            else:
                result = self.client.submit(
                    apply_registrations, current_param, self.current_session
                )
                _ = self.client.gather(result)

    def segment_masks(self, current_param, label):
        if "segmentedObjects" in current_param.param_dict.keys():
            operation = current_param.param_dict["segmentedObjects"]["operation"]
        else:
            operation = [""]

        if (
            label != "RNA"
            and current_param.param_dict["acquisition"]["label"] != "RNA"
            and "2D" in operation
        ):
            if not self.parallel:
                segment_masks(current_param, self.current_session)
            else:
                result = self.client.submit(
                    segment_masks, current_param, self.current_session
                )
                _ = self.client.gather(result)

    def segment_masks_3d(self, current_param, label):
        if (label in ("DAPI", "mask")) and "3D" in current_param.param_dict[
            "segmentedObjects"
        ]["operation"]:
            print_log(f"Making 3D image segmentations for label: {label}")
            print_log(f">>>>>>Label in functionCaller:{label}")

            _segment_sources_3d = SegmentMasks3D(
                current_param, self.current_session, parallel=self.parallel
            )
            _segment_sources_3d.segment_masks_3d()

    def segment_sources_3d(self, current_param, label):
        if (
            label == "barcode"
            and "3D" in current_param.param_dict["segmentedObjects"]["operation"]
        ):
            print_log(f"Making 3D image segmentations for label: {label}")
            print_log(f">>>>>>Label in functionCaller:{label}")

            _segment_sources_3d = SegmentSources3D(
                current_param, self.current_session, parallel=self.parallel
            )
            _segment_sources_3d.segment_sources_3d()

    def process_pwd_matrices(self, current_param, label):
        if label in ("DAPI", "mask"):
            if not self.parallel:
                process_pwd_matrices(current_param, self.current_session)
            else:
                result = self.client.submit(
                    process_pwd_matrices, current_param, self.current_session
                )
                _ = self.client.gather(result)


# =============================================================================
# FUNCTIONS
# =============================================================================


# filters barcode localization table
def filter_localizations(current_param, label):
    if label == "barcode":
        filter_localizations_instance = FilterLocalizations(current_param)
        filter_localizations_instance.filter_folder()


# filters barcode localization table
def register_localizations(current_param, label):
    if label == "barcode":
        register_localizations_instance = RegisterLocalizations(current_param)
        register_localizations_instance.register()


# build traces
def build_traces(current_param, label):
    if label == "barcode":
        build_traces_instance = BuildTraces(current_param)
        build_traces_instance.run()


# build matrices
def build_matrix(current_param, label):
    if label == "barcode":
        build_matrix_instance = BuildMatrix(current_param)
        build_matrix_instance.run()
