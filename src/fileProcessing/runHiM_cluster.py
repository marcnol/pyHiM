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
    parser.add_argument("--threads", help="Number of threads for parallel mode. None: sequential execution")
    parser.add_argument("--srun", help="Runs using srun", action="store_true")
    parser.add_argument("--xrun", help="Runs using bash", action="store_true")
    parser.add_argument("--sbatch", help="Runs using sbatch", action="store_true")

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

    if args.sbatch:
        runParameters["sbatch"] = args.sbatch
    else:
        runParameters["sbatch"] = False

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

    if args.threads:
        runParameters["threads"] = args.threads
    else:
        runParameters["threads"] = None

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

    if runParameters["threads"] is None:
        threads = ""
    else:
        threads = " --threads " + runParameters["threads"]

    if runParameters["cmd"] is None:
        cmdName=""
        CMD = ""
        jobNameExt = "_completePipeline"
    else:
        cmdName=runParameters["cmd"]
        CMD = " -C " + cmdName
        jobNameExt = "_" + cmdName


    if runParameters["sbatch"]:
        BATCH_file = ["#!/bin/bash"]
        SBATCH_header = [[
                "#!/bin/bash",\
                "#SBATCH "+memPerCPU,\
                "#SBATCH "+CPUsPerTask,\
                "#SBATCH "+nTasksCPU,\
                "#SBATCH --account="+runParameters["account"],\
                "#SBATCH --partition="+runParameters["partition"],\
                "#SBATCH --mail-user=marcnol@gmail.com "]]
        SBATCH_header.append(["",\
                "source /trinity/shared/apps/local/Python/Anaconda/3-5.1.0/etc/profile.d/conda.sh",\
                "conda activate pyHiM",""])

    for folder in folders:

        outputFile = runParameters["HOME"] + os.sep + "logs" + os.sep + runParameters["dataset"] + "_" + os.path.basename(folder) + "_" + cmdName +".log"
        jobName = os.path.basename(folder)+jobNameExt

        print("Folder to run: {}".format(folder))
        print("Output logfile: {}".format(outputFile))

        pyHiM = (
            "pyHiM.py -F "
            + folder
            + CMD
            + threads
            + " > "
            + outputFile
            )

        if not runParameters["sbatch"]:
            pyHiM = pyHiM + " &"

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

        if runParameters["sbatch"]:
            SBATCH_list = list()
            SBATCH_list = SBATCH_list + SBATCH_header[0]
            SBATCH_list.append("#SBATCH --job-name={}".format(jobName))
            SBATCH_list = SBATCH_list + SBATCH_header[1]
            SBATCH_list.append("\n# dataset: {}".format(jobName))
            SBATCH_list.append(
                "srun "
                + pyHiM
            )

        if runParameters["xrun"]:
            os.system(pyHiM)
        elif runParameters["srun"]:
            os.system(SRUN)


        if not runParameters["sbatch"]:
            print("Command to run: {}".format(SRUN))
            print("-"*50)
        elif runParameters["sbatch"]:
            print("SBATCH script:\n{}".format("\n".join(SBATCH_list)))
            print("-"*80)

            fileName="sbatch_script_{}.bash".format(jobName)
            with open(fileName, 'w') as f:
                for item in SBATCH_list:
                    f.write("{}\n".format(item))

            BATCH_file.append("sbatch {}".format(fileName))

    if runParameters["sbatch"]:

        print("*"*80)
        BATCH_file.append("\n")
        BASHscriptName="batch_script_{}.bash".format(runParameters["dataset"])
        with open(BASHscriptName, 'w') as f:
            for item in BATCH_file:
                f.write("{}\n".format(item))

        print("\nTo run master bash script:\n$ bash {}".format(BASHscriptName))
