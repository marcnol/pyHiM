#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Helpers for checking generated files against reference outputs."""

import os
from typing import Callable, Optional

from core.data_manager import extract_files


def relative_output_name(
    filepath: str, short_filename: str, extension: Optional[str]
) -> str:
    """Return the relative output name used by the regression fixtures."""
    filename = short_filename if extension is None else f"{short_filename}.{extension}"
    if "data" in filepath.split(os.sep):
        return os.path.join("data", filename)
    return filename


def reference_file_names(reference_outputs: str) -> list[str]:
    """List reference fixture names relative to the reference output root."""
    return [
        relative_output_name(filepath, short_filename, extension)
        for filepath, short_filename, extension in extract_files(reference_outputs)
    ]


def assert_reference_outputs_exist(
    generated_root: str,
    reference_outputs: str,
    compare: Callable[[str, str], None],
    aliases: Optional[dict[str, str]] = None,
):
    """Assert every reference output exists in generated outputs and compare it.

    The pipeline may create additional diagnostic files depending on parameters or
    dependency versions. Regression tests should therefore verify the canonical
    reference artifacts rather than fail only because extra files were emitted.
    """
    names = reference_file_names(reference_outputs)
    assert len(names) > 0
    for reference_name in names:
        generated_name = (
            aliases.get(reference_name, reference_name) if aliases else reference_name
        )
        tmp_file = os.path.join(generated_root, generated_name)
        out_file = os.path.join(reference_outputs, reference_name)
        assert os.path.exists(tmp_file), f"Missing generated output: {generated_name}"
        compare(tmp_file, out_file)
