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
import os
import glob
import argparse
from fileManagement import Parameters, folders, writeString2File
import shutil
import sys


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    if len(sys.argv) < 3:
        raise SystemExit("Not enough arguments")

    fileListString = sys.argv[1]
    destFolder = sys.argv[2]

    print("fileList = {} | destDir = {}".format(fileListString, destFolder))

    fileList = [x for x in glob.glob(fileListString)]

    if len(fileList) > 0:
        fileName = os.path.dirname(fileList[0]) + os.sep + "lndir.log"

        for file in fileList:

            newFile = destFolder + os.sep + os.path.basename(file)
            print("{}-->{}".format(file, newFile))

            command = "ln -s " + file + " " + newFile
            os.system(command)

            writeString2File(fileName, command, attribute="a")

        print("Linked {} files form {} to {}".format(len(fileList), os.path.dirname(fileList[0]), destFolder))

    else:
        print("File List is empty :(")
