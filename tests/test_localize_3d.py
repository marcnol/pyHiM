#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the non regression of Localize3D feature"""

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
    image_pixel_differences,
)

# Build a temporary directory
tmp_dir = tempfile.TemporaryDirectory()
# Define a "localize_3d" directory inside the temp dir
tmp_localize_3d_in = os.path.join(tmp_dir.name, "localize_3d")
# Copy the modes & IN/OUT structure for localize_3d inside the "localize_3d" temp dir
shutil.copytree("pyhim-small-dataset/localize_3d/IN", tmp_localize_3d_in)


def template_test_localize_3d(mode: str):
    """Check Localize3D feature with all possibilities"""
    inputs = os.path.join(tmp_localize_3d_in, mode)
    main(["-F", inputs, "-C", "localize_3d"])
    generated_align_images = os.path.join(inputs, "localize_3d")
    reference_outputs = f"pyhim-small-dataset/localize_3d/OUT/{mode}/segmentedObjects/"
    generated_files = extract_files(generated_align_images)
    reference_files = extract_files(reference_outputs)
    assert len(generated_files) == len(reference_files)
    for filepath, short_filename, extension in generated_files:
        if "data" in filepath.split(os.sep):
            filename = f"data{os.sep}{short_filename}.{extension}"
        else:
            filename = f"{short_filename}.{extension}"
        tmp_file = os.path.join(generated_align_images, filename)
        out_file = os.path.join(reference_outputs, filename)
        assert os.path.exists(out_file)
        if extension == "npy":
            assert compare_npy_files(tmp_file, out_file)
        elif extension == "png":
            assert image_pixel_differences(tmp_file, out_file)
        elif extension == "json":
            assert compare_line_by_line(tmp_file, out_file)
        elif extension == "dat":
            assert compare_line_by_line(
                tmp_file,
                out_file,
                line_start=len("e5c550de-d381-4b63-99d3-736ca7e549d9"),
                shuffled_lines=True,
            )
        elif extension == "table":
            assert compare_ecsv_files(tmp_file, out_file)
        else:
            raise ValueError(f"Extension file UNRECOGNIZED: {filepath}")


def test_full_stack_stardist():
    template_test_localize_3d("full_stack_stardist")


def test_reduce_planes_stardist():
    template_test_localize_3d("reduce_planes_stardist")


def test_reduce_planes_thresholding():
    template_test_localize_3d("reduce_planes_thresholding")
