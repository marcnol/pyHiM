#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and functions for pyHiM logging
"""

import json
import logging
import os
from datetime import datetime

from core.data_manager import save_json, write_string_to_file


class Logger:
    def __init__(self, root_folder, parallel=False, session_name=""):
        self.m_root_folder = root_folder
        self.m_log = Log(root_folder=root_folder, parallel=parallel)
        self.m_session = Session(root_folder, session_name)
        self.log_file = ""
        self.md_filename = ""
        self.setup_md_file()
        self.setup_logger()

    def setup_md_file(self):
        print("\n-----------------------------------------------------------------")
        print(f"$ root_folder: {self.m_root_folder}")
        begin_time = datetime.now()
        print_session_name(self.m_session.name)
        now = datetime.now()
        date_time = now.strftime("%d%m%Y_%H%M%S")

        self.log_file = (
            self.m_root_folder + os.sep + "HiM_analysis" + date_time + ".log"
        )
        self.md_filename = self.log_file.split(".")[0] + ".md"

        print(f"$ Hi-M analysis will be written to: {self.md_filename}")
        # initialises MD file
        write_string_to_file(
            self.md_filename,
            f"""# Hi-M analysis {begin_time.strftime("%Y/%m/%d %H:%M:%S")}""",
            "w",
        )

    def setup_logger(self):
        # creates output formats for terminal and log file
        formatter1 = logging.Formatter("%(asctime)s: %(levelname)s: %(message)s")
        formatter2 = logging.Formatter("%(message)s")

        # clears up any existing logger
        logger = logging.getLogger()
        logger.handlers = []
        for hdlr in logger.handlers[:]:
            if isinstance(hdlr, logging.FileHandler):
                logger.removeHandler(hdlr)

        # initializes handlers for terminal and file
        filehandler = logging.FileHandler(self.log_file, "w")
        stream_handler = logging.StreamHandler()

        logger.setLevel(logging.INFO)
        filehandler.setLevel(logging.INFO)
        stream_handler.setLevel(logging.INFO)

        logger.addHandler(stream_handler)
        logger.addHandler(filehandler)

        filehandler.setFormatter(formatter1)
        stream_handler.setFormatter(formatter2)


class Session:
    """Used to log planned tasks on this run"""

    def __init__(self, root_folder, name="dummy"):
        now = datetime.now()
        session_root_name = now.strftime("%d%m%Y_%H%M%S")
        self.file_name = root_folder + os.sep + "Session_" + session_root_name + ".json"
        self.name = name
        self.data = {}

    def save(self):
        """Saves session to file"""
        save_json(self.data, self.file_name)
        print(f"> Saved json session file to {self.file_name}")

    def add(self, key, value):
        """Add new task to session

        Parameters
        ----------
        key : str
            Task name
        value : str
            Task description
        """
        if key not in self.data:
            self.data[key] = value
        else:
            self.data[key] = [self.data[key], value]


class Log:
    """Hand made logging system to correspond for pyhim log needs"""

    def __init__(self, root_folder="./", parallel=False):
        now = datetime.now()
        date_time = now.strftime("%d%m%Y_%H%M%S")
        self.file_name = root_folder + os.sep + "logfile" + date_time + ".log"
        self.md_filename = self.file_name.split(".")[0] + ".md"
        self.parallel = parallel

    # saves to logfile, no display to cmd line
    def save(self, text=""):
        """saves to logfile, no display to cmd line

        Parameters
        ----------
        text : str, optional
            Text to write, by default ""
        """
        write_string_to_file(self.file_name, str(text), "a")

    def report(self, text, status="info"):
        """this function will output to cmd line and save in logfile

        Parameters
        ----------
        text : str
            text to write
        status : str, optional
            message status, by default "info"
        """
        if not self.parallel or status.lower() == "error":
            print(str(text))
            self.save("\n" + text)
        else:
            self.save("\n" + text)

    def add_simple_text(self, title):
        """Add simple text inside log file and print it to cmd line terminal

        Parameters
        ----------
        title : str
            Text to write
        """
        print(title)
        write_string_to_file(self.file_name, title, "a")


def print_log(message, status="INFO"):
    """
    Shows message to terminal and logs it to file

    Parameters
    ----------
    message : str
        message.
    status : str, optional
        either DEBUG, INFO or WARN. The default is 'INFO'.

    Returns
    -------
    None.

    """
    # print(message)

    if status == "INFO":
        logging.info(message)
    elif status == "WARN":
        logging.warning(message)
    elif status == "DEBUG":
        logging.debug(message)


def print_session_name(name: str):
    print_log(f"\n===================={name}====================\n")


def print_dashes():
    print_log("-----------------------------------------------------------------")
