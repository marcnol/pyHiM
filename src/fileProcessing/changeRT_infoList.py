#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:49:07 2020

@author: marcnol

Usage: changeRT_infolist.py first_RT_label new_RT_label.
Example changeRT_infolist.py RT33 RT95
"""

import os
import re
import sys

# import argparse


labels_to_process = [
    {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
    {"label": "barcode", "parameterFile": "infoList_barcode.json"},
    {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
    {"label": "RNA", "parameterFile": "infoList_RNA.json"},
]

PATTERN = r"(.*)\"(?P<RT>RT\d{1,3})\""


N_ARGS = len(sys.argv) - 1  # sys.argv[0] is base name
print("Total arguments passed: {}".format(N_ARGS))

if N_ARGS != 1:
    print("Wrong number of arguments!\n")
    print("Please specify the label of the new RT, as follows:\n")
    print("$ changeRT_infolist.py RT95\n")
    sys.exit(-1)
else:
    new_rt = sys.argv[1]


for label_to_process in labels_to_process:
    LABEL = label_to_process["label"]
    LABEL_PARAMETER_FILE = label_to_process["parameterFile"]
    print("**Modifying label {}: {}**".format(LABEL, LABEL_PARAMETER_FILE))

    with open(LABEL_PARAMETER_FILE, mode="r", encoding="utf-8") as f:
        # pylint: disable-next=invalid-name
        old_rt = ""
        for line in f.readlines():
            match = re.search(PATTERN, line)
            if match:
                old_rt = match.group("RT")
                break

    if old_rt == "":
        print(
            "Could not find a matching old_rt in file {}.".format(LABEL_PARAMETER_FILE)
        )
        print("Aborting.")
        sys.exit(-1)

    command_to_run_1 = (
        "sed -i '" + "s+" + old_rt + "+" + new_rt + "+g' " + LABEL_PARAMETER_FILE
    )
    print("Command: {}".format(command_to_run_1))

    return_value = os.system(command_to_run_1)
    print("Changing {} to {}. return_value {}.".format(old_rt, new_rt, return_value))
