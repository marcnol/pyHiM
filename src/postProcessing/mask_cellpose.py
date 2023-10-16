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
import numpy as np
from cellpose import models
from cellpose.io import imread


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

def run_cellpose_api(image_path, diam, cellprob, flow, stitch,gpu = True,):
    
    # model_type='cyto' or 'nuclei' or 'cyto2'
    model = models.Cellpose(gpu = gpu, model_type='cyto')

    # list of files
    files = [image_path]

    imgs = [imread(f) for f in files]

    # define CHANNELS to run segementation on
    channels = [[0,0]]

    # runs model    
    masks, flows, styles, diams = model.eval(imgs,
                                             channels=channels,
                                             diameter=None,
                                             cellprob_threshold = cellprob,
                                             flow_threshold=flow,
                                             stitch_threshold= stitch,
                                             )

    return masks

def run_cellpose(image_path, diam, cellprob, flow, stitch,folder_destination='segmentedObjects'):
    save_folder = os.path.dirname(image_path)

    '''
    command = (
        f"cellpose --verbose "
        + f"--image_path {image_path} "
        + "--use_gpu "
        + f"--chan 0 --diameter {diam} "
        + f"--stitch_threshold {stitch} "
        + f"--flow_threshold {flow} "
        + f"--cellprob_threshold {cellprob}"
    )
    '''
    
    command = (
        f"cellpose --verbose "
        + f"--image_path {image_path} --no_npy --save_tif "
        + "--use_gpu "
        + f"--chan 0 --diameter {diam} "
        + f"--stitch_threshold {stitch} "
        + f"--flow_threshold {flow} "
        + f"--cellprob_threshold {cellprob}"
    )
    
    print(f"$ will run: {command}")
    subprocess.run(command, shell=True)

    # moves image to new location
    mask_name = image_path.split(".")[0] + "_cp_masks.tif" #"_seg.npy"
    print(f"$ Mask image saved at: {mask_name}")
    
    new_mask_name = (
        save_folder + folder_destination  + os.sep + os.path.basename(image_path).split(".")[0] + "_Masks.tif"
    )
    
    if os.path.exists(new_mask_name):
        subprocess.run(f"mv {new_mask_name} {new_mask_name}_old", shell=True)
        print(f'Warning: File already exists, we moved it to: {new_mask_name}_old')

    move_img = f"mv -f {mask_name} {new_mask_name}"
    subprocess.run(move_img, shell=True)
    print(f"$ Moved segmentation file to : {new_mask_name}")

def process_images(cellprob=-8, flow=10, stitch=0.1, diam=50, files=list()):
    print(f"Parameters: diam={diam} | cellprob={cellprob} | flow={flow} | stitch={stitch}\n")
    
    folder_destination = "segmentedObjects"
    if ~os.path.exists(folder_destination):
        os.mkdir(folder_destination)
        
    if len(files) > 0:
        print("\n{} trace files to process= {}".format(len(files), "\n".join(map(str, files))))

        # iterates over traces in folder
        for file in files:
            print(f"> Analyzing image {file}")
            save_folder = os.path.dirname(file)

            # run_cellpose(file, diam, cellprob, flow, stitch,folder_destination=folder_destination)
            mask = run_cellpose_api(file, diam, cellprob, flow, stitch, gpu=True)

            new_mask_name = (
                save_folder + folder_destination  + os.sep + os.path.basename(file).split(".")[0] + "_Masks.npy"
            )
            
            np.save(new_mask_name,mask)

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
