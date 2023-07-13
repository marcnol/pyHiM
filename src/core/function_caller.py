#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for high level function calling
"""


import logging
import os
from datetime import datetime

from core.dask_cluster import DaskCluster
from core.pyhim_logging import (
    Log,
    Logger,
    Session,
    print_dashes,
    print_log,
    write_string_to_file,
)
from imageProcessing.alignImages import align_images, apply_registrations
from imageProcessing.alignImages3D import Drift3D
from imageProcessing.makeProjections import Project, make_projections
from imageProcessing.segmentMasks import segment_masks
from imageProcessing.segmentMasks3D import SegmentMasks3D
from imageProcessing.segmentSources3D import SegmentSources3D
from matrixOperations.alignBarcodesMasks import process_pwd_matrices
from matrixOperations.build_matrix import BuildMatrix
from matrixOperations.build_traces import BuildTraces
from matrixOperations.filter_localizations import FilterLocalizations
from matrixOperations.register_localizations import RegisterLocalizations


class Pipeline:
    """Class for high level function calling"""

    def __init__(self, data_m, cmd_list, global_param, is_parallel, logger):
        self.m_data_m = data_m
        self.cmds = cmd_list
        self.params = global_param
        self.parallel = is_parallel
        self.m_logger = logger
        self.m_dask = None
        self.tempo_var = {"makeProjections": Project}
        self.features = []
        self.init_features()

    def init_features(self):
        for command in self.cmds:
            if command in self.tempo_var:
                self.features.append(self.tempo_var.get(command)(self.params))

    def manage_parallel_option(self, feature, *args, **kwargs):
        if not self.parallel:
            feature(*args, **kwargs)
        else:
            result = self.m_dask.client.submit(feature, *args, **kwargs)
            _ = self.m_dask.client.gather(result)

    def lauch_dask_scheduler(self, threads_requested=25, maximum_load=0.8):
        if self.parallel:
            print_log(f"$ Requested {threads_requested} threads")

            self.m_dask = DaskCluster(threads_requested, maximum_load=maximum_load)

            self.m_dask.create_distributed_client()
            # self.client = dask_cluster_instance.client
            # self.cluster = dask_cluster_instance.cluster

    def find_files_to_process(self):
        pass

    def make_projections(self, current_param):
        self.manage_parallel_option(
            make_projections, current_param, self.m_logger.m_session
        )

    def align_images(self, current_param, label):
        if (
            label == "fiducial"
            and current_param.param_dict["acquisition"]["label"] == "fiducial"
        ):
            print_log(f"> Making image registrations for label: {label}")
            self.manage_parallel_option(
                align_images, current_param, self.m_logger.m_session
            )

    def align_images_3d(self, current_param, label):
        if (
            label == "fiducial"
            and "block3D" in current_param.param_dict["alignImages"]["localAlignment"]
        ):
            print_log(f"> Making 3D image registrations label: {label}")
            _drift_3d = Drift3D(
                current_param, self.m_logger.m_session, parallel=self.parallel
            )
            _drift_3d.align_fiducials_3d()

    def apply_registrations(self, current_param, label):
        if (
            label != "fiducial"
            and current_param.param_dict["acquisition"]["label"] != "fiducial"
        ):
            print_log(f"> Applying image registrations for label: {label}")
            self.manage_parallel_option(
                apply_registrations, current_param, self.m_logger.m_session
            )

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
            self.manage_parallel_option(
                segment_masks, current_param, self.m_logger.m_session
            )

    def segment_masks_3d(self, current_param, label):
        if (label in ("DAPI", "mask")) and "3D" in current_param.param_dict[
            "segmentedObjects"
        ]["operation"]:
            print_log(f"Making 3D image segmentations for label: {label}")
            print_log(f">>>>>>Label in functionCaller:{label}")

            _segment_sources_3d = SegmentMasks3D(
                current_param, self.m_logger.m_session, parallel=self.parallel
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
                current_param, self.m_logger.m_session, parallel=self.parallel
            )
            _segment_sources_3d.segment_sources_3d()

    def process_pwd_matrices(self, current_param, label):
        if label in ("DAPI", "mask"):
            self.manage_parallel_option(
                process_pwd_matrices, current_param, self.m_logger.m_session
            )

    def run(self):
        for feat in self.features:
            (required_data, required_ref, required_table) = feat.get_required_inputs()
            # reference = self.m_data_m.load_reference(required_ref)
            # table = self.m_data_m.load_table(required_table)
            files_to_process = self.m_data_m.get_inputs(required_data)
            self.m_data_m.create_folder(feat.out_folder)
            if self.parallel:
                pass
                # C'est KC
                # # dask
                # client = get_client()
                # threads = [
                #     client.submit(feat.run, f2p.load(), reference, table)
                #     for f2p in files_to_process
                # ]
                # print_log(f"$ Waiting for {len(threads)} threads to complete ")
                # for _, _ in enumerate(threads):
                #     wait(threads)
            else:
                for f2p in files_to_process:
                    data = f2p.load()
                    results = feat.run(data, f2p.m_label)
                    # results = feat.run(data, reference, table)
                    self.m_data_m.save_data(
                        results, feat.out_tags, feat.out_folder, f2p
                    )


# =============================================================================
# FUNCTIONS
# =============================================================================


def filter_localizations(current_param, label):
    """Filters barcode localization table

    Parameters
    ----------
    current_param : Parameters
        _description_
    label : str
        Only 'barcode' are accepted
    """
    if label == "barcode":
        filter_localizations_instance = FilterLocalizations(current_param)
        filter_localizations_instance.filter_folder()


def register_localizations(current_param, label):
    """Registers barcode localization table

    Parameters
    ----------
    current_param : Parameters
        _description_
    label : str
        Only 'barcode' are accepted
    """
    if label == "barcode":
        register_localizations_instance = RegisterLocalizations(current_param)
        register_localizations_instance.register()


def build_traces(current_param, label):
    """Build traces

    Parameters
    ----------
    current_param : Parameters
        _description_
    label : str
        Only 'barcode' are accepted
    """
    if label == "barcode":
        build_traces_instance = BuildTraces(current_param)
        build_traces_instance.run()


def build_matrix(current_param, label):
    """Build matrices

    Parameters
    ----------
    current_param : Parameters
        _description_
    label : str
        Only 'barcode' are accepted
    """
    if label == "barcode":
        build_matrix_instance = BuildMatrix(current_param)
        build_matrix_instance.run()
