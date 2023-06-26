#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pyHiM argument parser module"""

import argparse
import os

from core.parameters import print_dict
from core.pyhim_logging import print_log


def set_of_commands():
    """Available commands for pyHiM

    Returns
    -------
    frozenset
        Set of available commands
    """
    return frozenset(
        {
            "makeProjections",
            "appliesRegistrations",
            "alignImages",
            "alignImages3D",
            "segmentMasks",
            "segmentMasks3D",
            "segmentSources3D",
            "filter_localizations",
            "register_localizations",
            "build_traces",
            "build_matrix",
            "buildHiMmatrix",  # DEPRECATED
        }
    )


def default_2d_commands():
    """Default commands for 2D pipeline

    Returns
    -------
    frozenset
        Set of 2D commands
    """
    return frozenset(
        {
            "makeProjections",
            "alignImages",
            "appliesRegistrations",
            "segmentMasks",
            "filter_localizations",
            "build_traces",
            "build_matrix",
        }
    )


def default_3d_commands():
    """Default commands for 3D pipeline

    Returns
    -------
    frozenset
        Set of 3D commands
    """
    return frozenset(
        {
            "makeProjections",
            "alignImages",
            "alignImages3D",
            "segmentMasks3D",
            "segmentSources3D",
            "filter_localizations",
            "register_localizations",
            "build_traces",
            "build_matrix",
        }
    )


def him_parse_arguments(command_line_arguments):
    """Parse run arguments of pyhiM runtime

    Parameters
    ----------
    command_line_arguments : List[str]
        Used to give inputs for the runtime when you call this function like a module.
        For example, to test the pyHiM run from tests folder.

    Returns
    -------
    ArgumentParser.args
        An accessor of run arguments

    Raises
    ------
    SystemExit
        Command not found
    """
    parser = argparse.ArgumentParser()
    run_parameters = {}
    available_commands = set_of_commands()

    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument(
        "-S",
        "--stardist_basename",
        help="Replace all stardist_basename from infoList.json",
    )
    parser.add_argument(
        "-C",
        "--cmd",
        help="Comma-separated list of routines to run (without space !): \
                        makeProjections alignImages appliesRegistrations alignImages3D \
                        segmentMasks segmentMasks3D segmentSources3D buildHiMmatrix (DEPRECATED) \
                        filter_localizations register_localizations build_traces build_matrix",
    )
    parser.add_argument(
        "-T",
        "--threads",
        help="Number of threads to run in parallel mode. \
            If none, then it will run with one thread.",
    )
    args = parser.parse_args(command_line_arguments)

    print_log(
        "\n--------------------------------------------------------------------------"
    )

    if args.rootFolder:
        run_parameters["rootFolder"] = args.rootFolder
    else:
        # pylint: disable-next=consider-iterating-dictionary
        if "docker" in os.environ.keys():
            run_parameters["rootFolder"] = "/data"
            print_log(
                f"\n\n$ Running in docker, him_data: {run_parameters['rootFolder']}"
            )
        else:
            print_log("\n\n# him_data: NOT FOUND")
            run_parameters["rootFolder"] = os.getcwd()

    if args.stardist_basename:
        run_parameters["stardist_basename"] = args.stardist_basename
    else:
        run_parameters["stardist_basename"] = None

    if args.threads:
        run_parameters["threads"] = int(args.threads)
        run_parameters["parallel"] = True
    else:
        run_parameters["threads"] = 1
        run_parameters["parallel"] = False

    default_commands = default_3d_commands()
    if args.cmd:
        cmd_dict = {
            "2D": default_2d_commands,
            "3D": default_3d_commands,
        }
        run_parameters["cmd"] = cmd_dict.get(args.cmd, lambda: args.cmd.split(","))()
    else:
        run_parameters["cmd"] = default_commands

    for cmd in run_parameters["cmd"]:
        if cmd not in available_commands:
            print_log(
                f"\n\n# ERROR: {cmd} not found in list of available commands:\n\
                {available_commands}\n",
                status="WARN",
            )
            raise SystemExit

    print_dict(run_parameters)

    return run_parameters
