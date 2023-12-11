#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the non regression of Mask2D feature"""

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
# Define a "mask_2d" directory inside the temp dir
tmp_mask_2d_in = os.path.join(tmp_dir.name, "mask_2d")
# Copy the modes & IN/OUT structure for mask_2d inside the "mask_2d" temp dir
shutil.copytree("pyhim-small-dataset/mask_2d/IN", tmp_mask_2d_in)


def template_test_mask_2d(mode: str):
    """Check Mask2D feature with all possibilities"""
    inputs = os.path.join(tmp_mask_2d_in, mode)
    main(["-F", inputs, "-C", "mask_2d"])
    generated_align_images = os.path.join(inputs, "mask_2d")
    reference_outputs = f"pyhim-small-dataset/mask_2d/OUT/{mode}/segmentedObjects/"
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
        elif extension == "table" or extension == "dat":
            assert compare_ecsv_files(tmp_file, out_file)
        else:
            raise ValueError(f"Extension file UNRECOGNIZED: {filepath}")


def test_inhomogeneous_only():
    template_test_mask_2d("inhomogeneous_only")


def test_inhomogeneous_tesselation():
    template_test_mask_2d("inhomogeneous_tesselation")


def test_stardist_only():
    template_test_mask_2d("stardist_only")


def test_stardist_tesselation():
    template_test_mask_2d("stardist_tesselation")
