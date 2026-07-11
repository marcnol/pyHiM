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
    compare_dat_file_structure,
    compare_ecsv_files,
    compare_line_by_line,
    compare_mask_files,
    compare_npy_files,
    compare_npy_shape,
)
from tests.testing_tools.output import assert_reference_outputs_exist

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
    tmp_align_images = os.path.join(tmp_small_inputs, "alignImages")
    out_align_images = (
        "pyhim-small-dataset/resources/small_dataset/OUT/alignImages/data"
    )
    out_apply_register = (
        "pyhim-small-dataset/resources/small_dataset/OUT/appliesRegistrations/data"
    )

    def compare(tmp_file, out_file):
        extension = out_file.rsplit(".", 1)[-1] if "." in out_file else None
        if extension == "npy":
            assert compare_npy_shape(tmp_file, out_file)
        elif extension == "json":
            assert os.path.getsize(tmp_file) > 0
        else:
            assert compare_line_by_line(tmp_file, out_file, shuffled_lines=True)

    assert_reference_outputs_exist(
        tmp_align_images,
        out_align_images,
        compare,
        aliases={"shifts.json": "register_global.json"},
    )
    assert_reference_outputs_exist(tmp_align_images, out_apply_register, compare)


def test_align_images_3d():
    """Check register_local"""
    main(["-F", tmp_small_inputs, "-C", "register_global"])
    main(["-F", tmp_small_inputs, "-C", "register_local"])
    tmp_align_images = os.path.join(tmp_small_inputs, "alignImages/data")
    out_align_images = (
        "pyhim-small-dataset/resources/small_dataset/OUT/alignImages3D/data"
    )
    out_files = extract_files(out_align_images)
    assert len(out_files) > 0
    compared = 0
    for _, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_align_images, filename)
        if not os.path.exists(tmp_file) and filename == "shifts_block3D.dat":
            tmp_file = os.path.join(tmp_align_images, "register_global_block3D.dat")
        if not os.path.exists(tmp_file):
            continue
        out_file = os.path.join(out_align_images, filename)
        assert compare_dat_file_structure(tmp_file, out_file)
        compared += 1


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
        assert compare_mask_files(tmp_file, out_file)


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
