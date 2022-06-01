#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 10:55:13 2020

@author: marcnol
"""

import os, glob


def makeLink(new_folder, x, list1, ext):
    y = new_folder + os.sep + "_".join(list1) + "." + ext
    os.system("ln -s " + x + " " + y)
    return 1


rootFolderDAPI = "/home/marcnol/grey/rawData_2020/Exp_Combinatory_3_Tof/DAPI/Dapi_deconvolved/"
rootFolderRT = "/home/marcnol/grey/rawData_2020/Exp_Combinatory_3_Tof/RT/RT_deconvolved"
dest_folder = "/home/marcnol/grey/rawData_2020/Exp_Combinatory_3_Tof/rois"

filesDAPI = glob.glob(rootFolderDAPI + os.sep + "*.tif")
filesRT = glob.glob(rootFolderRT + os.sep + "*.tif")


listROIs = []
for x in filesDAPI:
    ROI = os.path.basename(x).split(".")[0].split("_")[3]
    if ROI not in listROIs:
        listROIs.append(ROI)

for ROI in listROIs:

    new_folder = dest_folder + os.sep + "ROI" + ROI
    if not os.path.isdir(new_folder):
        os.mkdir(new_folder)

    filesDAPIFiltered = [x for x in filesDAPI if os.path.basename(x).split(".")[0].split("_")[3] == ROI]

    filesRTFiltered = [x for x in filesRT if os.path.basename(x).split(".")[0].split("_")[3] == ROI]

    filesRTFiltered = [x for x in filesRT if os.path.basename(x).split(".")[0].split("_")[5] == ROI]

    # print([os.path.basename(x) for x in filesDAPIFiltered])
    # print([os.path.basename(x) for x in filesRTFiltered])

    for x in filesDAPIFiltered:
        y = new_folder + os.sep + os.path.basename(x)
        os.system("ln -s " + x + " " + y)

    counter = 0
    for x in filesRTFiltered:
        filename = os.path.basename(x)
        list0 = filename.split(".")[0].split("_")
        ext = filename.split(".")[1]

        if list0[-1] == "ch01":
            list1 = list0[0:2] + ["".join(list0[2] + list0[3])] + list0[5:]
            counter += makeLink(new_folder, x, list1, ext)

        elif list0[-1] == "ch02":
            list1 = list0[0:2] + ["".join(list0[2] + list0[4])] + list0[5:]
            list1 = [x if x != "ch02" else "ch01" for x in list1]
            counter += makeLink(new_folder, x, list1, ext)

        elif list0[-1] == "ch00":
            list1 = list0[0:2] + ["".join(list0[2] + list0[3])] + list0[5:]
            counter += makeLink(new_folder, x, list1, ext)

            list2 = list0[0:2] + ["".join(list0[2] + list0[4])] + list0[5:]
            counter += makeLink(new_folder, x, list2, ext)

    print("{} RT links made for ROI {}".format(counter, ROI))
