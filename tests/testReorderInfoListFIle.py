#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 10:47:24 2021

@author: marcnol
"""
import json

file_name = '/home/marcnol/data/Embryo_debug_dataset/test_dataset/infoList_model.json'
with open(file_name, encoding="utf-8") as json_file:
    param0 = json.load(json_file)
current_param = param0["common"]
labels = param0["labels"]

ordered_list = [" "]*len(labels.keys())
for i,label in enumerate(labels.keys()):
    order = labels[label]["order"]
    print('order={}'.format(order))
    ordered_list[order-1]=label


current_param["labels"] = ordered_list


''' need to:
    - incorporate changes in labels list to the current_param list that will be returned.
    - call on Parameters to read the new parameters file before starting the loop in pyHiM
    - iterate over 'current_param['labels']' within pyHiM
    - refactor Parameters() so that we can 'reset' the parameters structure depending on the label used.
    - this may break dependencies as before Parameters() was receiving a filename for a label, now it will receive a filename + label...
    - this will  require possibly  minor changes but I am not positive...:
        - cleanHiM
        - zipHiMrun
        - alignBarcodesMasks.py
        - processSNDchannel.py


'''