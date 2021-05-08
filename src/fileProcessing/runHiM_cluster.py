#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 19:28:21 2021

@author: marcnol

Lauches slurm srun job

In the command line, run as
$ runHiM_cluster.py


"""
import os, sys
import glob
import argparse

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-D", "--dataset", help="dataset: folder in ~/scratch")
    parser.add_argument("-N", "--nCPU", help="Number of CPUs/Task")
    parser.add_argument("--memPerCPU", help="Memory required per allocated CPU in Mb")
    parser.add_argument("-C", "--cmd", help="Comma-separated list of routines to run (order matters !): makeProjections alignImages \
                        appliesRegistrations alignImages3D segmentMasks \
                        segmentSources3D refitBarcodes3D \
                        localDriftCorrection projectBarcodes buildHiMmatrix")
    parser.add_argument("-R", "--run", help="Deletes folders, MD files, LOG files", action="store_true")
    parser.add_argument("-F", "--dataFolder", help="Folder with data. Default: ~/scratch")
    parser.add_argument("-S", "--singleDataset", help="Folder for single Dataset.")

    args = parser.parse_args()

    runParameters = dict()
    runParameters["HOME"] = os.environ["HOME"]

    if args.dataset:
        runParameters["dataset"] = args.dataset
    else:
        runParameters["dataset"] = None

    if args.nCPU:
        runParameters["nCPU"] = int(args.nCPU)
    else:
        runParameters["nCPU"] = 1

    if args.memPerCPU:
        runParameters["memPerCPU"] = args.memPerCPU
    else:
        runParameters["memPerCPU"] = None

    if args.cmd:
        runParameters["cmd"] = args.cmd
    else:
        runParameters["cmd"] = None

    if args.run:
        runParameters["run"] = args.run
    else:
        runParameters["run"] = False

    if args.dataFolder:
        runParameters["dataFolder"] = args.dataFolder
    else:
        runParameters["dataFolder"] = runParameters["HOME"] + os.sep + "scratch"

    if args.singleDataset:
        runParameters["singleDataset"] = args.singleDataset
    else:
        runParameters["singleDataset"] = None

    runParameters["partition"] = "tests"
    runParameters["account"] = "episcope"

    print("Parameters loaded: {}\n".format(runParameters))

    if runParameters["dataset"] is None:
        print("ERROR: No dataset provided!")
        raise SystemExit

    rootFolder = runParameters["dataFolder"] + os.sep + runParameters["dataset"]

    if runParameters["singleDataset"] is None:
        folders = glob.glob(rootFolder + os.sep + "*")
        folders0 = [x for x in folders if os.path.isdir(x)] # keeps only folders
        folders = [x for x in folders0 if os.path.exists(x+os.sep+"infoList.json")]
    else:
        folders0 = folders = [runParameters["singleDataset"]]
        runParameters["dataset"] = os.path.basename(runParameters["singleDataset"])

    print("*"*50)
    print("$ Dataset: {}".format(runParameters["dataset"]))
    print("$ Folder: {}".format(rootFolder))
    print("$ Number of CPUs: {}".format(runParameters["nCPU"]))
    print("$ Command: {}".format(runParameters["cmd"]))
    print("*"*50)

    print("\n\n$ Found {} folders in {}".format(len(folders0),rootFolder))
    print("$ Of these, {} contained an infoList.json file and will be processed".format(len(folders)))
    print("Folders to process: {}".format(folders))
    print("$ Scheduling {} jobs...".format(len(folders)))
    print("-"*50)

    if runParameters["memPerCPU"] is None:
        memPerCPU = ""
    else:
        memPerCPU = " --mem_per_cpu=" + runParameters["memPerCPU"]

    for folder in folders:

        outputFile = runParameters["HOME"] + os.sep + "logs" + os.sep + runParameters["dataset"] + "_" + os.path.basename(folder) + ".log"
        print("Folder to run: {}".format(folder))
        print("Output logfile: {}".format(outputFile))

        if runParameters["cmd"] is None:
            CMD = ""
            jobName = os.path.basename(folder)+"_completePipeline"
        else:
            CMD = " -C " + runParameters["cmd"]
            jobName = os.path.basename(folder) + "_" + runParameters["cmd"]

        SRUN = (
            "srun --account="
            + runParameters["account"]
            + " --partition="
            + runParameters["partition"]
            + " --job-name="
            + jobName
            + " --cpus-per-task "
            + str(runParameters["nCPU"])
            + memPerCPU
            + " --mail-user=marcnol@gmail.com pyHiM.py -F "
            + folder
            + CMD
            + " > "
            + outputFile
            + " &"
        )

        print("Command to run: {}".format(SRUN))
        print("-"*50)

        if runParameters["run"]:
            os.system(SRUN)

    # NEED TO NOW INVOKE TO RUN
