#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data manager module

Manage writing, reading and checking data.
"""

import os


def extract_files(root: str):
    """Extract recursively file informations of all files into a given directory.
    Note: filename is the name without extension

    Parameters
    ----------
    root : str
        The name of root directory

    Returns
    -------
    List[Tuple(str,str,str)]
        List of file informations: (filepath, filename, extension)
    """
    files = []
    # Iterate into dirname and each subdirectories dirpath, dirnames, filenames
    for dirpath, dirnames, filenames in os.walk(root):
        for filename in filenames:
            split_filename = filename.split(".")
            extension = split_filename.pop()
            short_filename = ".".join(split_filename)
            filepath = os.path.join(dirpath, filename)
            files.append((filepath, short_filename, extension))

        if len(dirnames) > 0:
            print("[INFO] Inside:")
            print(dirpath)
            print("[INFO] Subdirectories detected:")
            print(dirnames)

    return files


class DataManager:
    def __init__(self, data_path: str, stardist_basename: str):
        self.m_data_path = self.__set_data_path(data_path)
        self.m_stardist_basename = str(stardist_basename)
        self.all_files = extract_files(self.m_data_path)
        self.user_parameter = None
        self.data_images = []
        self.data_tables = []
        self.dispatch_files()

    @staticmethod
    def __set_data_path(data_path):
        if data_path:
            return str(data_path)
        return os.getcwd()

    def dispatch_files(self):
        img_ext = ["tif", "tiff"]
        table_ext = ["csv", "ecsv", "dat"]
        for path, name, ext in self.all_files:
            if ext in img_ext:
                self.data_images.append((path, name, ext))
            elif ext in table_ext:
                self.data_tables.append((path, name, ext))
            elif ext == "json":
                self.user_parameter = (path, name, ext)
            else:
                print(
                    f"Unrecognized data file: {name}.{ext}\n>>>Inside this folder: {path}"
                )
