#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 18:05:08 2023

@author: marcnol

Usage:
    
    $ ls scan_*ROI.tif | run_cellpose.py --pipe 

or just for a single file

    run_cellpose.py --input scan_001_DAPI_001_ROI.tif 

"""
import os
import sys
import subprocess
import select
import argparse


def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Name of input trace file.")
    parser.add_argument("--cellprob", help="cellprob threshold. Default = -8.")
    parser.add_argument("--flow", help="flow threshold. Default = 10.")
    parser.add_argument("--stitch", help="stitch threshold. Default = 0.1.")
    parser.add_argument("--diam", help="diameter. Default = 50.")
    parser.add_argument("--pipe", help="inputs Trace file list from stdin (pipe)", action="store_true")

    p = {}

    args = parser.parse_args()

    if args.input:
        p["input"] = args.input
    else:
        p["input"] = None

    if args.cellprob:
        p["cellprob"] = float(args.cellprob)
    else:
        p["cellprob"] = -8

    if args.flow:
        p["flow"] = float(args.flow)
    else:
        p["flow"] = 10

    if args.stitch:
        p["stitch"] = float(args.stitch)
    else:
        p["stitch"] = 0.1

    if args.diam:
        p["diam"] = float(args.diam)
    else:
        p["diam"] = 50

    p["files"] = []
    if args.pipe:
        p["pipe"] = True
        if select.select(
            [
                sys.stdin,
            ],
            [],
            [],
            0.0,
        )[0]:
            p["files"] = [line.rstrip("\n") for line in sys.stdin]
        else:
            print("Nothing in stdin")
    else:
        p["pipe"] = False
        p["files"] = [p["input"]]

    return p


def run_cellpose(image_path, diam, cellprob, flow, stitch):
    save_folder = os.path.dirname(image_path)
    
    # if save_folder == "":
        # save_folder = "./"

    command = (
        f"cellpose --verbose "
        + f"--image_path {image_path} "
        + "--use_gpu "
        + f"--chan 0 --diameter {diam} "
        + f"--stitch_threshold {stitch} "
        + f"--flow_threshold {flow} "
        + f"--cellprob_threshold {cellprob}"
    )
    print(f"$ will run: {command}")
    subprocess.run(command, shell=True)

    # moves image to new location
    mask_name = image_path.split(".")[0] + "_seg.npy"
    print(f"$ Mask image saved at: {mask_name}")
    
    new_mask_name = (
        save_folder + "segmentedObjects" + os.sep + os.path.basename(image_path).split(".")[0] + "_Masks.npy"
    )
    
    if os.path.exists(new_mask_name):
        subprocess.run(f"mv {new_mask_name} {new_mask_name}_old", shell=True)
        print(f'Warning: File already exists, we moved it to: {new_mask_name}_old')

    move_img = f"mv -f {mask_name} {new_mask_name}"
    subprocess.run(move_img, shell=True)
    print(f"$ Moved segmentation file to : {new_mask_name}")

def process_images(cellprob=-8, flow=10, stitch=0.1, diam=50, files=list()):
    print(f"Parameters: diam={diam} | cellprob={cellprob} | flow={flow} | stitch={stitch}\n")

    if len(files) > 0:
        print("\n{} trace files to process= {}".format(len(files), "\n".join(map(str, files))))

        # iterates over traces in folder
        for file in files:
            print(f"> Analyzing image {file}")

            run_cellpose(file, diam, cellprob, flow, stitch)


# =============================================================================
# MAIN
# =============================================================================


def main():
    # [parsing arguments]
    p = parseArguments()

    print("Remember to activate environment: conda activate cellpose!\n")

    # [loops over lists of datafolders]
    process_images(
        cellprob=p["cellprob"],
        flow=p["flow"],
        stitch=p["stitch"],
        diam=p["diam"],
        files=p["files"],
    )

    print("Finished execution")


if __name__ == "__main__":
    main()
