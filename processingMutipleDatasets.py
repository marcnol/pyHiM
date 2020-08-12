#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 14:01:28 2020

@author: marcnol


This script is written to run several datasets at once.
It will take as arguments:
    - name of the root directory containing the 'Embryo_X' files
    - the numbers of embryos to be analysed,e.g. 1 2 3 will analyse Embryo_1, Embryo_2 and Embryo_3
    
Example usage:
    
$ processingMultipleDatasets.py /mnt/grey/DATA/rawData_2020/Experiment_4/deconvolved_DAPI/ 0 1 2 4 5


"""


import os
import sys

nArgs = len(sys.argv)
print("Total arguments passed: {}".format(nArgs))
EmbryoTag = "Embryo_"
if nArgs > 2:
    rootDir = sys.argv[1]
    print("parameters> rootFolder: {}".format(rootDir))

    for i in range(2, nArgs):
        print("Processing Embryo #{}".format(sys.argv[i]))
        command2Run1 = "nice -19 processingPipeline.py -F " + rootDir + EmbryoTag + sys.argv[i]
        os.system(command2Run1)
        command2Run2 = "zipHiM_run.py -F " + rootDir + EmbryoTag + str(i)
        os.system(command2Run2)
        print("Commands: {}\n{}".format(command2Run1, command2Run2))


else:
    print("not enough arguments.")
    print(
        "Example Usage: \n\n $ processingMultipleDatasets.py /mnt/grey/DATA/rawData_2020/Experiment_4/deconvolved_DAPI/ 0 1 2 4 5"
    )
    print(
        "\n analyses Embryo_0, Embryo_1, Embryo_2, Embryo_4, Embryo_5 in folder /mnt/grey/DATA/rawData_2020/Experiment_4/deconvolved_DAPI/\n"
    )
