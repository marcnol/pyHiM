#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest

from core import run_args as ra


def test_check_consistency():
    # Test SystemExit of data_path
    with pytest.raises(SystemExit, match=r".*input/path/does/not/exist.*"):
        run_args = ra.RunArgs(["-F", "input/path/does/not/exist"])
        run_args._check_consistency()

    # Test SystemExit of cmd_list
    with pytest.raises(SystemExit, match=r"'ratafia' .*"):
        run_args = ra.RunArgs(["-C", "ratafia"])
        run_args._check_consistency()

    # Test SystemExit of threads
    with pytest.raises(SystemExit, match=r".*Number of threads.*"):
        run_args = ra.RunArgs(["-T", "-1"])
        run_args._check_consistency()

    # Test SystemExit of stardist_basename
    with pytest.raises(SystemExit, match=r".*in/path/stradist/basename.*"):
        run_args = ra.RunArgs(["-S", "in/path/stradist/basename"])
        run_args._check_consistency()
