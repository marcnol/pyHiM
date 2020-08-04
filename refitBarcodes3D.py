#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:50:49 2020

@author: marcnol

Purpose: Performs 3D Gaussian fits for barcodes

steps: 
    - iterate over ROIs
    - load mask for ROI <i>
    - check if mask available and load mask file 
    - iterate over cycles <i>
    - load 3D file for barcode <i>
    - iterate over masks <j>
    - 3D gaussian fit of barcodes <k> in mask <j>
    - determine if we keep or not
    - store in database.
    


"""


