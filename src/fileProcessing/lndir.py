#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 23:56:27 2020

@author: marcnol

make links of files in a directory to a second directory

In the command line, run as 

Example:

$ lndir.py "/home/marcnol/Repositories/pyHiM/*py" ~/Downloads/test

Make sure that the first argument has quotation marks if you use wildcards!

"""

import glob
import os
import sys

from fileManagement import write_string_to_file

# =============================================================================
# MAIN
# =============================================================================

def main():

    if len(sys.argv) < 3:
        raise SystemExit("Not enough arguments")

    file_list_string = sys.argv[1]
    dest_folder = sys.argv[2]

    print("file_list = {} | destDir = {}".format(file_list_string, dest_folder))

    file_list = list(glob.glob(file_list_string))

    if len(file_list) > 0:
        file_name = os.path.dirname(file_list[0]) + os.sep + "lndir.log"

        for file in file_list:

            new_file = dest_folder + os.sep + os.path.basename(file)
            print("{}-->{}".format(file, new_file))

            command = "ln -s " + file + " " + new_file
            os.system(command)

            write_string_to_file(file_name, command, attribute="a")

        print(
            "Linked {} files form {} to {}".format(
                len(file_list), os.path.dirname(file_list[0]), dest_folder
            )
        )

    else:
        print("File List is empty :(")

if __name__ == "__main__":
    main()