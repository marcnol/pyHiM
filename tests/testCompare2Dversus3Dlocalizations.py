#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:42:24 2021

@author: marcnol
"""

from astropy.table import Table
import matplotlib.pyplot as plt

rootFolder = "/home/marcnol/data/Embryo_debug_dataset/test_dataset/segmentedObjects/"
file1 = rootFolder + "segmentedObjects_barcode.dat"
file2 = rootFolder + "segmentedObjects_3D_barcode.dat"


barcodeMap2D = Table.read(file1, format="ascii.ecsv")
barcodeMap3D = Table.read(file2, format="ascii.ecsv")

names=['zcentroid','ycentroid','xcentroid']
zyx2D = [barcodeMap2D[x] for x in names]
zyx3D = [barcodeMap3D[x] for x in names]


plt.plot(zyx2D[1],zyx2D[2],'or')
plt.plot(zyx3D[1],zyx3D[2],'+b')