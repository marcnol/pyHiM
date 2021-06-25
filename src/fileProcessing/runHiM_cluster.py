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
def readArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-D", "--dataset", help="dataset: name of folder with data within dataFolder")
    parser.add_argument("-F", "--dataFolder", help="Folder with data. Default: ~/scratch")
    parser.add_argument("-S", "--singleDataset", help="Folder for single Dataset.")
    parser.add_argument("-A", "--account", help="Provide your account name. Default: episcope.")
    parser.add_argument("-P", "--partition", help="Provide partition name. Default: tests")
    parser.add_argument("-N", "--nCPU", help="Number of CPUs/Task")
    parser.add_argument("--memPerCPU", help="Memory required per allocated CPU in Mb")
    parser.add_argument("--nodelist", help="Specific host names to include in job allocation.")
    parser.add_argument("-T1","--nTasksNode", help="Number of tasks per node.")
    parser.add_argument("-T2","--nTasksCPU", help="Number of tasks per CPU")
    parser.add_argument("-C", "--cmd", help="Comma-separated list of routines to run (order matters !): makeProjections alignImages \
                        appliesRegistrations alignImages3D segmentMasks \
                        segmentSources3D refitBarcodes3D \
                        localDriftCorrection projectBarcodes buildHiMmatrix")
    parser.add_argument("-R", "--srun", help="Runs using srun", action="store_true")
    parser.add_argument("--xrun", help="Runs using bash", action="store_true")

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
        runParameters["nCPU"] = None

    if args.memPerCPU:
        runParameters["memPerCPU"] = args.memPerCPU
    else:
        runParameters["memPerCPU"] = None

    if args.cmd:
        runParameters["cmd"] = args.cmd
    else:
        runParameters["cmd"] = None

    if args.xrun:
        runParameters["xrun"] = args.xrun
    else:
        runParameters["xrun"] = False

    if args.srun:
        runParameters["srun"] = args.srun
    else:
        runParameters["srun"] = False

    if args.dataFolder:
        runParameters["dataFolder"] = args.dataFolder
    else:
        runParameters["dataFolder"] = runParameters["HOME"] + os.sep + "scratch"

    if args.singleDataset:
        runParameters["singleDataset"] = args.singleDataset
    else:
        runParameters["singleDataset"] = None

    if args.account:
        runParameters["account"] = args.account
    else:
        runParameters["account"] = "episcope"

    if args.partition:
        runParameters["partition"] = args.partition
    else:
        runParameters["partition"] = "defq"

    if args.nodelist:
        runParameters["nodelist"] = args.nodelist
    else:
        runParameters["nodelist"] = None

    if args.nTasksNode:
        runParameters["nTasksNode"] = args.nTasksNode
    else:
        runParameters["nTasksNode"] = None

    if args.nTasksCPU:
        runParameters["nTasksCPU"] = args.nTasksCPU
    else:
        runParameters["nTasksCPU"] = None

    print("Parameters loaded: {}\n".format(runParameters))

    return runParameters

if __name__ == "__main__":

    runParameters = readArguments()

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
        # runParameters["dataset"] = os.path.basename(runParameters["singleDataset"])

    folders.sort()

    print("*"*50)
    print("$ Dataset: {}".format(runParameters["dataset"]))
    print("$ Folder: {}".format(rootFolder))
    print("$ Number of CPUs: {}".format(runParameters["nCPU"]))
    print("$ Command: {}".format(runParameters["cmd"]))
    print("$ Account: {}".format(runParameters["account"]))
    print("$ Partition: {}".format(runParameters["partition"]))

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

    if runParameters["nodelist"] is None:
        nodelist = ""
    else:
        nodelist = " --nodelist=" + runParameters["nodelist"]

    if runParameters["nCPU"] is None:
        CPUsPerTask = ""
    else:
        CPUsPerTask = " --cpus-per-task " + str(runParameters["nCPU"])
        
    if runParameters["nTasksCPU"] is None:
        nTasksCPU = ""
    else:
        nTasksCPU = " --ntasks-per-core=" + runParameters["nTasksCPU"]

    if runParameters["nTasksNode"] is None:
        nTasksNode = ""
    else:
        nTasksNode = " --ntasks-per-node=" + runParameters["nTasksNode"]

    if runParameters["cmd"] is None:
        cmdName=""
        CMD = ""
        jobNameExt = "_completePipeline"
    else:
        cmdName=runParameters["cmd"]
        CMD = " -C " + cmdName
        jobNameExt = "_" + cmdName

    for folder in folders:

        outputFile = runParameters["HOME"] + os.sep + "logs" + os.sep + runParameters["dataset"] + "_" + os.path.basename(folder) + "_" + cmdName +".log"
        jobName = os.path.basename(folder)+jobNameExt

        print("Folder to run: {}".format(folder))
        print("Output logfile: {}".format(outputFile))

        pyHiM = (
            "pyHiM.py -F "
            + folder
            + CMD
            + " > "
            + outputFile
            + " &"
            )

        SRUN = (
            "srun --account="
            + runParameters["account"]
            + " --partition="
            + runParameters["partition"]
            + " --job-name="
            + jobName
            + CPUsPerTask
            + nodelist
            + nTasksCPU
            + nTasksNode
            + memPerCPU
            + " --mail-user=marcnol@gmail.com "
            + pyHiM
        )

        if runParameters["xrun"]:
            os.system(pyHiM)
        elif runParameters["srun"]:
            os.system(SRUN)

        print("Command to run: {}".format(SRUN))
        print("-"*50)

