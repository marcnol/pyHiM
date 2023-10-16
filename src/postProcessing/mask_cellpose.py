#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 18:05:08 2023

@author: marcnol

Usage:
    
    $ ls scan_*ROI.tif | mask_cellpose.py --gpu

or just for a single file

    mask_cellpose.py --input scan_001_DAPI_001_ROI.tif --gpu
    

"""
import os
import sys
import subprocess
import select
import argparse
import numpy as np
from cellpose import models, core
from cellpose.io import imread


def parseArguments():

    parser=argparse.ArgumentParser(
        description='''mask_cellpose.py will segment a TIF file using predermined parameters using cellpose. \n These parameters can be
        changed using arguments. See below ''',
        epilog="""All is well that ends well.""")
    
    parser.add_argument("--input", help="Name of input trace file.")
    parser.add_argument(
        "--gpu", help="If used it will use gpu mode", action="store_true"
    )
    parser.add_argument(
        "--cli", help="It will call the CLI command instead of the API (which sometimes crashes). Default = False", action="store_true"
    )
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

    if args.gpu:
        p["gpu"] = [True, None]
    else:
        p["gpu"] = [False, None]

    if args.cli:
        p["cli"] = True
    else:
        p["cli"] = False
        
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

def run_cellpose_api(image_path, diam, cellprob, flow, stitch,gpu = [False,None]):
    
    # model_type='cyto' or 'nuclei' or 'cyto2'
    if gpu[0]:
        model = models.Cellpose(gpu = gpu[1], model_type='cyto')
    else:
        model = models.Cellpose(model_type='cyto')

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
                                             batch_size=8,
                                             )

    return masks[0]

def run_cellpose(image_path, diam, cellprob, flow, stitch,folder_destination='segmentedObjects',gpu = [False,None]):
    save_folder = os.path.dirname(image_path)

    command = (
        f"cellpose --verbose "
        + f"--image_path {image_path} --no_npy --save_tif "
        + f"--chan 0 --diameter {diam} "
        + f"--stitch_threshold {stitch} "
        + f"--flow_threshold {flow} "
        + f"--cellprob_threshold {cellprob}"
    )

    if gpu[0]:
        command = command + " --use_gpu"
        
    print(f"$ will run: {command}")
    subprocess.run(command, shell=True)

    mask_name = image_path.split(".")[0] + "_cp_masks.tif" #"_seg.npy"
    print(f"$ Reads TIF mask image: {mask_name}")
    
    if os.path.exists(mask_name):
        img = imread(mask_name)
        return img        

def process_images(cellprob=-8, flow=10, stitch=0.1, diam=50, files=list(), gpu = [False,None], cli=False):
    print(f"Parameters: diam={diam} | cellprob={cellprob} | flow={flow} | stitch={stitch}\n")
    
    folder_destination = "segmentedObjects"
        
    try:
        os.mkdir(folder_destination)
    except FileExistsError:
        print(">>> Output folder exists")
        
    if len(files) > 0:
        print("\n{} trace files to process= {}".format(len(files), "\n".join(map(str, files))))

        # iterates over traces in folder
        for file in files:
            print(f"> Analyzing image {file}")
            save_folder = os.path.dirname(file)

            if cli:
                mask = run_cellpose(file, diam, cellprob, flow, stitch,folder_destination=folder_destination,gpu=gpu)
            else:
                mask = run_cellpose_api(file, diam, cellprob, flow, stitch, gpu=gpu)
    
            new_mask_name = (
                save_folder + folder_destination  + os.sep + os.path.basename(file).split(".")[0] + "_Masks.npy"
            )

            if os.path.exists(new_mask_name):
                subprocess.run(f"mv {new_mask_name} {new_mask_name}_old", shell=True)
                print(f'Warning: File already exists, we moved it to: {new_mask_name}_old')
        
            print(f"> Saving output mask image to {new_mask_name}")            
            np.save(new_mask_name,mask)

# =============================================================================
# MAIN
# =============================================================================


def main():
    # [parsing arguments]
    p = parseArguments()

    print("Remember to activate environment: conda activate cellpose!\n")

    if p["gpu"][0]:
        use_GPU = core.use_gpu()
        p["gpu"][1] = use_GPU
        print('>>> GPU activated? {}'.format(use_GPU))
    else:
        print('>>> will use CPU')

    # [loops over lists of datafolders]
    process_images(
        cellprob=p["cellprob"],
        flow=p["flow"],
        stitch=p["stitch"],
        diam=p["diam"],
        files=p["files"],
        gpu = p["gpu"],
        cli = p["cli"],
    )

    print("Finished execution")


if __name__ == "__main__":
    main()
