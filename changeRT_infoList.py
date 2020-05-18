#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:49:07 2020

@author: marcnol
"""

import os
import sys
import argparse

nArgs=len(sys.argv)
print("Total arguments passed: {}".format(nArgs))

oldRT='RT33'
newRT='RT95'


labels2Process = [{'label':'fiducial', 'parameterFile': 'infoList_fiducial.json'},
                  {'label': 'barcode', 'parameterFile': 'infoList_barcode.json'},
                  {'label':'DAPI', 'parameterFile': 'infoList_DAPI.json'},
                  {'label':'RNA', 'parameterFile': 'infoList_RNA.json'}]
                    
if nArgs == 3: 
    oldRT=sys.argv[1]
    newRT=sys.argv[2]
    print("Changing {} to {}".format(oldRT, newRT))
    
    for ilabel in range(len(labels2Process)):
        label=labels2Process[ilabel]['label']
        labelParameterFile=labels2Process[ilabel]['parameterFile']
        print("**Modifying label {}: {}**".format(label,labelParameterFile))
        command2Run1='sed -i \'' + 's+' + oldRT + '+' + newRT + '+g\' ' + labelParameterFile
        print('Command: {}'.format(command2Run1))

        os.system(command2Run1)
        
        
        
else:
    print('not enough arguments.')