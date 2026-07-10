#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the non regression of RegisterLocal feature"""

import os
import shutil
import tempfile

# sys.path.append("..")
from pyHiM import main
from tests.testing_tools.comparison import (
    compare_dat_file_structure,
    compare_ecsv_files,
    compare_line_by_line,
    compare_npy_files,
    image_pixel_differences,
)
from tests.testing_tools.output import assert_reference_outputs_exist

# Build a temporary directory
tmp_dir = tempfile.TemporaryDirectory()
# Define a "register_local" directory inside the temp dir
tmp_register_local_in = os.path.join(tmp_dir.name, "register_local")
# Copy the modes & IN/OUT structure for register_local inside the "register_local" temp dir
shutil.copytree("pyhim-small-dataset/register_local/IN", tmp_register_local_in)


def template_test_register_local(mode: str):
    """Check RegisterLocal feature with all possibilities"""
    inputs = os.path.join(tmp_register_local_in, mode)
    main(["-F", inputs, "-C", "register_local"])
    generated_register_local = os.path.join(inputs, "alignImages")
    reference_outputs = f"pyhim-small-dataset/register_local/OUT/{mode}/alignImages/"

    aliases = {"data/shifts_block3D.dat": "data/register_global_block3D.dat"}

    def compare(tmp_file, out_file):
        extension = out_file.rsplit(".", 1)[-1] if "." in out_file else None
        if extension == "npy":
            assert compare_npy_files(tmp_file, out_file)
        elif extension == "png":
            assert image_pixel_differences(tmp_file, out_file)
        elif extension == "json":
            assert compare_line_by_line(tmp_file, out_file)
        elif extension == "table":
            assert compare_ecsv_files(tmp_file, out_file, shuffled_lines=True)
        elif extension == "dat":
            assert compare_dat_file_structure(tmp_file, out_file)
        else:
            raise ValueError(f"Extension file UNRECOGNIZED: {out_file}")

    assert_reference_outputs_exist(
        generated_register_local, reference_outputs, compare, aliases=aliases
    )


def test_with_global_done():
    template_test_register_local("with_global")


def test_without_register_global():
    template_test_register_local("alone")
