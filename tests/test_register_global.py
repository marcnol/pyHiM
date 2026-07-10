#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the non regression of RegisterGlobal feature
"""

import os
import shutil
import tempfile

# sys.path.append("..")
from pyHiM import main
from tests.testing_tools.comparison import (
    compare_ecsv_files,
    compare_line_by_line,
    compare_npy_files,
    image_pixel_differences,
)
from tests.testing_tools.output import assert_reference_outputs_exist

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

    aliases = {"data/shifts.json": "data/register_global.json"}

    def compare(tmp_file, out_file):
        extension = out_file.rsplit(".", 1)[-1] if "." in out_file else None
        if extension == "npy":
            assert compare_npy_files(tmp_file, out_file)
        elif extension == "png":
            assert image_pixel_differences(tmp_file, out_file)
        elif extension == "json":
            assert compare_line_by_line(tmp_file, out_file)
        elif extension == "table":
            assert compare_ecsv_files(tmp_file, out_file)
        else:
            raise ValueError(f"Extension file UNRECOGNIZED: {out_file}")

    assert_reference_outputs_exist(
        generated_align_images, reference_outputs, compare, aliases=aliases
    )


def test_global_alignement():
    template_test_register_global("global")


def test_align_by_block():
    template_test_register_global("block")
