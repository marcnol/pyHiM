#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pyHiM argument parser module"""


import os
from argparse import ArgumentParser

from core.pyhim_logging import print_log


def _parse_run_args(command_line_arguments):
    """Parse run arguments of pyHiM runtime

    Parameters
    ----------
    command_line_arguments : List[str]
        Used to give inputs for the runtime when you call this function like a module.
        For example, to test the pyHiM run from tests folder.

    Returns
    -------
    ArgumentParser.args
        An accessor of run arguments
    """
    parser = ArgumentParser()

    # Data path
    parser.add_argument(
        "-F",
        "--rootFolder",
        type=str,
        default=os.getcwd(),
        help="Folder path with input images",
    )
    parser.add_argument(
        "-C",
        "--cmd",
        help="Comma-separated list of routines to run (without space !): \
                        makeProjections alignImages appliesRegistrations alignImages3D \
                        segmentMasks segmentMasks3D segmentSources3D buildHiMmatrix (DEPRECATED) \
                        filter_localizations register_localizations build_traces build_matrix",
    )
    # Number of threads
    parser.add_argument(
        "-T",
        "--threads",
        type=int,
        default=1,
        help="Thread number to run with parallel mode. \
            If none or 1, then it will run with sequential mode.",
    )
    parser.add_argument(
        "-S",
        "--stardist_basename",
        type=str,
        default=None,
        help="Replace all stardist_basename from infoList.json",
    )

    return parser.parse_args(command_line_arguments)


class RunArgs:
    """Store and check run arguments"""

    def __init__(self, command_line_arguments):
        print_log("\n-----------------------------------------------------------------")
        parsed_args = _parse_run_args(command_line_arguments)
        self.data_path = parsed_args.rootFolder
        self._is_docker()
        self.cmd_list = self.parse_cmd(parsed_args.cmd)
        self.thread_nbr = parsed_args.threads
        self.parallel = self.thread_nbr > 1
        self.stardist_basename = parsed_args.stardist_basename
        self._check_consistency()
        self.print_loaded_args()

    def _is_docker(self):
        """Change the data path is run inside docker"""
        if "docker" in os.environ:
            self.data_path = "/data"
            print_log(f"\n\n$ Running in docker, him_data: {self.data_path}")

    def _check_consistency(self):
        if not os.path.isdir(self.data_path):
            raise SystemExit(f"Input data path ({self.data_path}) NOT FOUND.")
        available_commands = self.get_available_commands()
        low_available_cmds = [a_c.lower() for a_c in available_commands]
        for cmd in self.cmd_list:
            if cmd.lower() not in low_available_cmds:
                print_log(
                    f"\n\n# ERROR: {cmd} isn't an available command like:\n\
                    {available_commands}\n",
                    status="WARN",
                )
                raise SystemExit(
                    f"'{cmd}' isn't an available command like:\n\
                    {available_commands}"
                )

        if not isinstance(self.thread_nbr, int) or self.thread_nbr < 1:
            raise SystemExit(f"Number of threads ({self.thread_nbr}): INVALID.")

        if self.stardist_basename and not os.path.isdir(self.stardist_basename):
            raise SystemExit(f"Stardist basename ({self.stardist_basename}) NOT FOUND.")

    def print_loaded_args(self):
        """Print parameters in your shell terminal

        Parameters
        ----------
        dictionary : dict
            Parameters dictionary
        """

        def print_tab_spacer(name: str, val: str):
            spacer = "\t" * (3 - int(len(name) / 8))
            print("\t" + name + spacer + val)

        print("\n$ Loaded arguments:")
        print_tab_spacer("rootFolder", str(self.data_path))
        print_tab_spacer("stardist_basename", str(self.stardist_basename))
        print_tab_spacer("threads", str(self.thread_nbr))
        print_tab_spacer("parallel", str(self.parallel))
        print_tab_spacer("cmd", str(self.cmd_list))

    @classmethod
    def parse_cmd(cls, cmd):
        """Parse the input command list give by the user as a string (comma-separated)

        Parameters
        ----------
        cmd : str
            A comma-separated human-list of commands

        Returns
        -------
        List[str]
            A Python list of commands
        """
        default_commands = cls.get_3d_commands()
        if cmd:
            if cmd == "2D":
                cmd_list = cls.get_2d_commands()
            elif cmd == "3D":
                cmd_list = cls.get_3d_commands()
            else:
                cmd_list = cmd.split(",")
        else:
            cmd_list = default_commands
        return cmd_list

    @staticmethod
    def get_available_commands():
        """Available commands for pyHiM

        Returns
        -------
        frozenset
            Set of available commands
        """
        return (
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
        )

    @staticmethod
    def get_2d_commands():
        """Default commands for 2D pipeline

        Returns
        -------
        frozenset
            Set of 2D commands
        """
        return (
            "makeProjections",
            "alignImages",
            "appliesRegistrations",
            "segmentMasks",
            "filter_localizations",
            "build_traces",
            "build_matrix",
        )

    @staticmethod
    def get_3d_commands():
        """Default commands for 3D pipeline

        Returns
        -------
        frozenset
            Set of 3D commands
        """
        return (
            "makeProjections",
            "alignImages",
            "alignImages3D",
            "segmentMasks3D",
            "segmentSources3D",
            "filter_localizations",
            "register_localizations",
            "build_traces",
            "build_matrix",
        )
