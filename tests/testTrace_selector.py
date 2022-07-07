#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 12:48:02 2022

@author: marcnol
"""

from astropy.table import Table
import collections
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 30}

matplotlib.rc('font', **font)

file = '/home/marcnol/grey/users/marcnol/test_HiM/testDataset_17RTs/001_ROI' + '/buildsPWDmatrix/Trace_3D_barcode_mask:mask1_ROI:1.ecsv'

trace_table = Table.read(file, format="ascii.ecsv")
trace_table

print(f"$ Number of lines in trace: {len(trace_table)}")

trace_by_ID = trace_table.group_by('Trace_ID')
print(trace_by_ID)


trace_lengths = list()
trace_unique_barcodes = list()
trace_repeated_barcodes = list()
number_unique_barcodes = list()
number_repeated_barcodes = list()

for sub_trace_table in trace_by_ID.groups:
    trace_lengths.append(len(sub_trace_table))
    
    unique_barcodes = list(set(sub_trace_table['Barcode #']))
    trace_unique_barcodes.append(unique_barcodes)
    number_unique_barcodes.append(len(unique_barcodes))
    
    repeated_barcodes = [item for item, count in collections.Counter(sub_trace_table['Barcode #']).items() if count > 1]
    trace_repeated_barcodes.append(repeated_barcodes)
    number_repeated_barcodes.append(len(repeated_barcodes))
    
    print(f"\ntrace length {len(sub_trace_table)}")
    print(f"number unique barcodes {len(unique_barcodes)}")
    print(f"number unique barcodes {len(repeated_barcodes)}: {repeated_barcodes}")
    

distributions = [trace_lengths,number_unique_barcodes,number_repeated_barcodes]
axis_x_labels = ['number of barcodes','number of unique barcodes','number of repeated barcodes']
number_plots = len(distributions)

fig = plt.figure(constrained_layout=True)
im_size = 10
fig.set_size_inches(( im_size*number_plots,im_size))
gs = fig.add_gridspec(1,number_plots)
axes = [fig.add_subplot(gs[0,i]) for i in range(number_plots)]


for axis, distribution,xlabel in zip(axes,distributions,axis_x_labels):
    axis.hist(distribution, alpha=.3)
    axis.set_xlabel(xlabel)
    axis.set_ylabel('counts')
    axis.set_title('median = '+str(np.median(distribution)))
    
plt.savefig('test.png')
