#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 14:49:22 2022

@author: marcnol

converts NPY to TIFF

"""

# from apifish.stack.io import read_image as read_image
import os
import glob
import sys
import select

import numpy as np
from tifffile import TiffWriter

files = sys.argv[1:]

if select.select([sys.stdin, ], [], [], 0.0)[0]:
    print("Have data!")

    piped_files = list()    
    for line in sys.stdin:
        piped_files.append(line.rstrip("\n"))    

    if len(piped_files):
        files = piped_files
        print(f"Read {piped_files} files list from stdin!\n")
else:
    print("No stdin")
  
  
    
    
if len(files)>0:
    
    print(f"{len(files)} files to process: ")
    for file in files:
        print(f"\n{os.path.basename(file)}")
        
        data = np.load(file)
        # data.astype('float')
        
        file_out = file.rsplit('.')[0]+'.tif'
        print(f"output file: {file_out}")
        
        with TiffWriter(file_out) as tif:
            tif.save(data.astype(np.uint16))        
        
        print(f"Finished processing {len(files)} files")
        
else:
    print("Please provide list of NPY files to process! none provided")
    
    print("\n Usage:\n")
    print("with find:")    
    print("npy_to_tiff $(find -name \"*ch0*_2d_registered.npy\")\n\n")
    print("with pipes:")
    print("ls *ch0*_2d_registered.npy | npy_to_tiff")        