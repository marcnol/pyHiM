#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and functions for pyHiM logging
"""

import logging
import os
from datetime import datetime


class Log:
    """Hand made logging system to correspond for pyhim log needs"""

    def __init__(self, root_folder="./", parallel=False):
        now = datetime.now()
        date_time = now.strftime("%d%m%Y_%H%M%S")
        self.file_name = root_folder + os.sep + "logfile" + date_time + ".log"
        self.markdown_filename = self.file_name.split(".")[0] + ".md"
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


def write_string_to_file(file_name, text_to_output, attribute="a"):
    """write string to the log file

    Parameters
    ----------
    file_name : str
        log file
    text_to_output : str
        text to write in log file
    attribute : str, optional
        Open file mode option, by default "a"
    """
    with open(file_name, mode=attribute, encoding="utf-8") as file_handle:
        file_handle.write(f"{text_to_output}\n")
