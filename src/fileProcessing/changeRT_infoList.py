#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:49:07 2020

@author: marcnol

Usage: changeRT_infolist.py first_RT_label new_RT_label.
Example changeRT_infolist.py RT33 RT95
"""

import os
import sys
import re

# import argparse


labels2Process = [
    {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
    {"label": "barcode", "parameterFile": "infoList_barcode.json"},
    {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
    {"label": "RNA", "parameterFile": "infoList_RNA.json"},
]

pattern = r"(.*)\"(?P<RT>RT\d{1,3})\""


nArgs = len(sys.argv) - 1  # sys.argv[0] is base name
print("Total arguments passed: {}".format(nArgs))

if nArgs != 1:
    print("Wrong number of arguments!\n")
    print("Please specify the label of the new RT, as follows:\n")
    print("$ changeRT_infolist.py RT95\n")
    sys.exit(-1)
else:
    newRT = sys.argv[1]


for ilabel in range(len(labels2Process)):
    label = labels2Process[ilabel]["label"]
    labelParameterFile = labels2Process[ilabel]["parameterFile"]
    print("**Modifying label {}: {}**".format(label, labelParameterFile))

    with open(labelParameterFile, "r") as f:
        oldRT = ""
        for line in f.readlines():
            match = re.search(pattern, line)
            if match:
                # print('Match found:', match.group("RT"))
                oldRT = match.group("RT")
                break

    if oldRT == "":
        print("Could not find a matching oldRT in file {}.".format(labelParameterFile))
        print("Aborting.")
        sys.exit(-1)

    command2Run1 = "sed -i '" + "s+" + oldRT + "+" + newRT + "+g' " + labelParameterFile
    print("Command: {}".format(command2Run1))

    returnValue = os.system(command2Run1)
    print("Changing {} to {}. returnValue {}.".format(oldRT, newRT, returnValue))
