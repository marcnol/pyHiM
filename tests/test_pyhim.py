#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check the non regression of modules called inside the jupyter notebook on pyHiM documentation"""

import os
import shutil
import tempfile

from core.data_manager import extract_files

# sys.path.append("..")
from pyHiM import main
from tests.testing_tools.comparison import (
    compare_ecsv_files,
    compare_line_by_line,
    compare_npy_files,
)

# build a temporary directory
tmp_dir = tempfile.TemporaryDirectory()
tmp_resources = os.path.join(tmp_dir.name, "resources")
shutil.copytree("pyhim-small-dataset/resources", tmp_resources)
tmp_small_inputs = os.path.join(tmp_resources, "small_dataset/IN")
tmp_traces_inputs = os.path.join(tmp_resources, "traces_dataset/IN")


def test_make_projections():
    """Check 'project'"""
    main(["-F", tmp_small_inputs, "-C", "project"])
    tmp_z_project = os.path.join(tmp_small_inputs, "zProject/data")
    out_z_project = (
        "pyhim-small-dataset/resources/small_dataset/OUT/makeProjections/data"
    )
    out_files = extract_files(out_z_project)
    assert len(out_files) > 0
    for _, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_z_project, filename)
        out_file = os.path.join(out_z_project, filename)
        assert compare_npy_files(tmp_file, out_file)


def test_register_global():
    """Check register_global"""
    main(["-F", tmp_small_inputs, "-C", "register_global"])
    tmp_align_images = os.path.join(tmp_small_inputs, "alignImages/data")
    out_align_images = (
        "pyhim-small-dataset/resources/small_dataset/OUT/alignImages/data"
    )
    out_apply_register = (
        "pyhim-small-dataset/resources/small_dataset/OUT/appliesRegistrations/data"
    )
    out_files = extract_files(out_align_images) + extract_files(out_apply_register)
    tmp_files = extract_files(tmp_align_images)
    assert len(out_files) > 0
    assert len(out_files) == len(tmp_files)
    for file_path, short_filename, extension in out_files:
        if extension:
            filename = short_filename + "." + extension
        else:
            filename = short_filename
        tmp_file = os.path.join(tmp_align_images, filename)
        if extension == "npy":
            assert compare_npy_files(tmp_file, file_path)
        else:
            assert compare_line_by_line(tmp_file, file_path, shuffled_lines=True)


def test_align_images_3d():
    """Check register_local"""
    main(["-F", tmp_small_inputs, "-C", "register_local"])
    tmp_align_images = os.path.join(tmp_small_inputs, "alignImages/data")
    out_align_images = (
        "pyhim-small-dataset/resources/small_dataset/OUT/alignImages3D/data"
    )
    out_files = extract_files(out_align_images)
    assert len(out_files) > 0
    for _, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_align_images, filename)
        out_file = os.path.join(out_align_images, filename)
        assert compare_line_by_line(tmp_file, out_file, shuffled_lines=True)


def test_segment_masks_3d():
    """Check mask_3d"""
    main(["-F", tmp_small_inputs, "-C", "mask_3d"])
    tmp_segmented_objects = os.path.join(tmp_small_inputs, "segmentedObjects/data")
    out_segmented_objects = (
        "pyhim-small-dataset/resources/small_dataset/OUT/segmentMasks3D/data"
    )
    out_files = extract_files(out_segmented_objects)
    assert len(out_files) > 0
    for _, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_segmented_objects, filename)
        out_file = os.path.join(out_segmented_objects, filename)
        assert compare_npy_files(tmp_file, out_file)


# TODO: Find a way to test this module
# def test_segment_sources_3d():
#     """Check localize_3d"""
#     main(["-F", tmp_small_inputs, "-C", "localize_3d"])
#     tmp_segmented_objects = os.path.join(tmp_small_inputs, "segmentedObjects")
#     out_segmented_objects = "pyhim-small-dataset/resources/small_dataset/OUT/segmentSources3D/"
#     out_files = extract_files(out_segmented_objects)
#     assert len(out_files) > 0
#     for _, short_filename, extension in out_files:
#         filename = short_filename + "." + extension
#         tmp_file = os.path.join(tmp_segmented_objects, filename)
#         out_file = os.path.join(out_segmented_objects, filename)
#         assert compare_ecsv_files(tmp_file, out_file,columns_to_remove= ["Buid","id"], shuffled_lines=True)


def test_build_traces():
    """Check build_traces"""
    main(["-F", tmp_traces_inputs, "-C", "build_traces"])
    tmp_builds_pwd_matrix = os.path.join(tmp_traces_inputs, "buildsPWDmatrix/data")
    out_builds_pwd_matrix = (
        "pyhim-small-dataset/resources/traces_dataset/OUT/build_traces/data"
    )
    out_files = extract_files(out_builds_pwd_matrix)
    assert len(out_files) > 0
    for _, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_builds_pwd_matrix, filename)
        out_file = os.path.join(out_builds_pwd_matrix, filename)
        assert compare_ecsv_files(tmp_file, out_file, columns_to_remove=["Trace_ID"])


def test_build_matrix():
    """Check build_matrix"""
    main(["-F", tmp_traces_inputs, "-C", "build_matrix"])
    tmp_builds_pwd_matrix = os.path.join(tmp_traces_inputs, "buildsPWDmatrix/data")
    out_builds_pwd_matrix = (
        "pyhim-small-dataset/resources/traces_dataset/OUT/build_matrix/data"
    )
    out_files = extract_files(out_builds_pwd_matrix)
    assert len(out_files) > 0
    for _, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_builds_pwd_matrix, filename)
        out_file = os.path.join(out_builds_pwd_matrix, filename)
        if extension == "npy":
            assert compare_npy_files(tmp_file, out_file, shuffled_plans=True)
        elif extension == "ecsv":
            assert compare_line_by_line(tmp_file, out_file)
