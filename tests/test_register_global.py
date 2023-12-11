#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the non regression of RegisterGlobal feature
"""

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
# Define a "register_global" directory inside the temp dir
tmp_register_global_in = os.path.join(tmp_dir.name, "register_global")
# Copy the modes & IN/OUT structure for register_global inside the "register_global" temp dir
shutil.copytree("pyhim-small-dataset/register_global/IN", tmp_register_global_in)


def template_test_register_global(mode: str):
    """Check RegisterGlobal feature with all possibilities"""
    inputs = os.path.join(tmp_register_global_in, mode)
    main(["-F", inputs, "-C", "register_global"])
    generated_align_images = os.path.join(inputs, "register_global")
    reference_outputs = f"pyhim-small-dataset/register_global/OUT/{mode}/alignImages/"
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
        elif extension == "table":
            assert compare_ecsv_files(tmp_file, out_file)
        else:
            raise ValueError(f"Extension file UNRECOGNIZED: {filepath}")


def test_global_alignement():
    template_test_register_global("global")


def test_align_by_block():
    template_test_register_global("block")
