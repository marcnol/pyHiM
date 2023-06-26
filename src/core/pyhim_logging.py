#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and functions for pyHiM logging
"""

import json
import logging
import os
from datetime import datetime


class Session:
    """Used to log planned tasks on this run"""

    def __init__(self, root_folder, name="dummy"):
        now = datetime.now()
        session_root_name = now.strftime("%d%m%Y_%H%M%S")
        self.file_name = root_folder + os.sep + "Session_" + session_root_name + ".json"
        self.name = name
        self.data = {}

    def load(self):
        """Loads existing session"""
        self.data = load_json(self.file_name)
        print(f"Session information read: {self.file_name}")

    def save(self):
        """Saves session to file"""
        save_json(self.file_name, self.data)
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

    def clear_data(self):
        """Reset data attribute with an empty dict"""
        self.data = {}


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


def save_json(file_name, data):
    """Save a python dict as a JSON file

    Parameters
    ----------
    file_name : str
        Output JSON file name
    data : dict
        Data to save
    """
    with open(file_name, mode="w", encoding="utf-8") as json_f:
        json.dump(data, json_f, ensure_ascii=False, sort_keys=True, indent=4)


def load_json(file_name):
    """Load a JSON file like a python dict

    Parameters
    ----------
    file_name : str
        JSON file name

    Returns
    -------
    dict
        Python dict
    """
    if os.path.exists(file_name):
        with open(file_name, encoding="utf-8") as json_file:
            data = json.load(json_file)
    else:
        data = {}
    return data
