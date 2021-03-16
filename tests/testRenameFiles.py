#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 10:55:13 2020

@author: marcnol
"""

import os, glob


def makeLink(newFolder, x, list1, ext):
    y = newFolder + os.sep + "_".join(list1) + "." + ext
    os.system("ln -s " + x + " " + y)
    return 1


rootFolderDAPI = "/home/marcnol/grey/rawData_2020/Exp_Combinatory_3_Tof/DAPI/Dapi_deconvolved/"
rootFolderRT = "/home/marcnol/grey/rawData_2020/Exp_Combinatory_3_Tof/RT/RT_deconvolved"
destFolder = "/home/marcnol/grey/rawData_2020/Exp_Combinatory_3_Tof/ROIs"

filesDAPI = glob.glob(rootFolderDAPI + os.sep + "*.tif")
filesRT = glob.glob(rootFolderRT + os.sep + "*.tif")


listROIs = []
for x in filesDAPI:
    ROI = os.path.basename(x).split(".")[0].split("_")[3]
    if ROI not in listROIs:
        listROIs.append(ROI)

for ROI in listROIs:

    newFolder = destFolder + os.sep + "ROI" + ROI
    if not os.path.isdir(newFolder):
        os.mkdir(newFolder)

    filesDAPIFiltered = [x for x in filesDAPI if os.path.basename(x).split(".")[0].split("_")[3] == ROI]

    filesRTFiltered = [x for x in filesRT if os.path.basename(x).split(".")[0].split("_")[3] == ROI]

    filesRTFiltered = [x for x in filesRT if os.path.basename(x).split(".")[0].split("_")[5] == ROI]

    # print([os.path.basename(x) for x in filesDAPIFiltered])
    # print([os.path.basename(x) for x in filesRTFiltered])

    for x in filesDAPIFiltered:
        y = newFolder + os.sep + os.path.basename(x)
        os.system("ln -s " + x + " " + y)

    counter = 0
    for x in filesRTFiltered:
        filename = os.path.basename(x)
        list0 = filename.split(".")[0].split("_")
        ext = filename.split(".")[1]

        if list0[-1] == "ch01":
            list1 = list0[0:2] + ["".join(list0[2] + list0[3])] + list0[5:]
            counter += makeLink(newFolder, x, list1, ext)

        elif list0[-1] == "ch02":
            list1 = list0[0:2] + ["".join(list0[2] + list0[4])] + list0[5:]
            list1 = [x if x != "ch02" else "ch01" for x in list1]
            counter += makeLink(newFolder, x, list1, ext)

        elif list0[-1] == "ch00":
            list1 = list0[0:2] + ["".join(list0[2] + list0[3])] + list0[5:]
            counter += makeLink(newFolder, x, list1, ext)

            list2 = list0[0:2] + ["".join(list0[2] + list0[4])] + list0[5:]
            counter += makeLink(newFolder, x, list2, ext)

    print("{} RT links made for ROI {}".format(counter, ROI))
