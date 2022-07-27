#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import ttk
from tkinter import END
from functools import partial
from info_parameters import *
from function_parameters import *
from tkinter import messagebox
import os
import sys

# initialised when entries are created, dictionary that store the class entry for the parameter,
# the value and the type of the value
entries_dic = {}

# dictionary that stores the user parameters entered in the entries
user_values_dic = {}


main_window = tk.Tk()
main_window.title("pyHiM parameters")
main_window.minsize(width=1000, height=480)
# ------------------------------------Looking for current and script directories--------------------------------------------

current_dir = os.getcwd()
script_dir = os.path.dirname(os.path.realpath(__file__))


# ---------------------------------------------Functions------------------------------------------------


def display_help(param):
    message = help_dic[param]
    help_label.config(text=message)


def restore_setting():
    default_Entry_dic, _ =import_parameters(script_dir)
    dapi_ch.set(default_Entry_dic["dapi_ch"])
    dapiFid_ch.set(default_Entry_dic["dapiFid_ch"])
    barcode_ch.set(default_Entry_dic["barcode_ch"])
    barcodeFid_ch.set(default_Entry_dic["barcodeFid_ch"])
    mask_ch.set(default_Entry_dic["mask_ch"])
    maskFid_ch.set(default_Entry_dic["maskFid_ch"])
    rna_ch.set(default_Entry_dic["rna_ch"])
    # rnaFid_ch.set(default_Entry_dic["rnaFid_ch"])
    for key, list_values in entries_dic.items():
        if list_values[0].get() != "ch00" and list_values[0].get() != "ch01" and list_values[0].get() != "ch02":
            list_values[0].delete(first=0, last=END)
            list_values[0].insert(0, string=default_Entry_dic[key])
            # list_values[0].insert(0, string=list_values[1])


def save_setting(entries_dic, user_values_dic):
    if check_settings(entries_dic):
        # update user values in entries_dic:
        for key, list_values in entries_dic.items():
            entered_value = list_values[0].get()
            user_values_dic[key] = entered_value
        update_infoList(user_values_dic, infoList_full, current_dir)
        messagebox.showinfo("Info Message", "The settings have been saved.")


# --------------------------------import values parameters from a infoList.json-----------------------------------------
infolist_partial, infoList_full = import_parameters(script_dir, current_dir)

# -----------------------------------------Save and restore button------------------------------------------------------

# Restore button
restore_button = tk.Button(main_window, text='Restore default settings', command=restore_setting)
restore_button.grid(row=27, column=0, columnspan=3, pady=5, padx=5)

# Save button
save_button = tk.Button(main_window, text='Save settings', command=partial(save_setting, entries_dic, user_values_dic))
save_button.grid(row=27, column=2, columnspan=3, pady=5, padx=5)


# --------------------------------------------Tab creation----------------------------------------------------
notebook = ttk.Notebook(main_window)
notebook.grid(row=0, column=0, rowspan=11, columnspan=5, pady=10)

tab1 = ttk.Frame(notebook, width=700, height=700)
tab1.grid(row=0, column=0, rowspan=10, columnspan=5)

tab2 = ttk.Frame(notebook, width=700, height=700)
tab2.grid(row=0, column=0, rowspan=10, columnspan=5)

notebook.add(tab1, text='Standard Settings')
notebook.add(tab2, text='Expert Settings')


# -------------------------Label Frame--------------------------------
# label frame (box) for quick help box
help_LabelFrame = tk.LabelFrame(main_window, text="Quick Help")
help_LabelFrame.grid(row=0, column=6, pady=5, padx=5)

# label frame (box) for Acquisition parameters (tab1)
acquisition_labelFrame = tk.LabelFrame(tab1, text="1. Acquisition parameters")
acquisition_labelFrame.grid(row=0, column=0, columnspan=6, pady=5, padx=5)

# label frame (box) for AlignImages parameters (tab1)
alignImages_LabelFrame = tk.LabelFrame(tab1, text="2. AlignImages parameters")
alignImages_LabelFrame.grid(row=11, column=0, columnspan=6, pady=5, padx=5, sticky='WE')

# label frame (box) for AlignImages parameters (tab2)
alignImages_LabelFrame2 = tk.LabelFrame(tab2, text="1. AlignImages parameters")
alignImages_LabelFrame2.grid(row=0, column=0, columnspan=6, pady=5, padx=5, sticky='WE')

# label frame (box) for BuildsPWDmatrix parameters (tab2)
buildsPWDmatrix_LabelFrame = tk.LabelFrame(tab2, text="2. BuildsPWDmatrix parameters")
buildsPWDmatrix_LabelFrame.grid(row=1, column=0, columnspan=6, pady=5, padx=5, sticky='WE')
#
# label frame (box) for SegmentedObjects parameters (tab2)
segmentedObjects_LabelFrame = tk.LabelFrame(tab2, text="3. SegmentedObjects parameters")
segmentedObjects_LabelFrame.grid(row=7, column=0, columnspan=6, pady=5, padx=5, sticky='WE')

# label frame (box) for Labels parameters (tab2)
labels_LabelFrame = tk.LabelFrame(tab2, text="4. Labels parameters")
labels_LabelFrame.grid(row=9, column=0, columnspan=6, pady=5, padx=5, sticky='WE')

# ---------------------------Help Box-----------------------------------
help_label = tk.Label(help_LabelFrame,
                      text="Help will come to those who ask for it",
                      height=20,
                      width=50,
                      justify='center',
                      wraplength=350)
help_label.grid(row=0, column=6, rowspan=11)

# ---------------------------Image pyHiM----------------------------------------------------------
img_path = script_dir + os.sep + "pyHiM_image.png"
img = tk.PhotoImage(file = img_path)
img = img.zoom(10)
img = img.subsample(15)
label_img = tk.Label(main_window, image=img)
label_img.grid(row=10, column=6, pady=10)


# -----------------------------------------Help Button--------------------------------------------
# Dapi Channel Help Button (tab1):
dapi_ch_HelpButton = tk.Button(acquisition_labelFrame, text="?", command=partial(display_help, "DAPI_channel"))
dapi_ch_HelpButton.grid(row=0, column=2)

# Dapi Fiducial Channel Help Button (tab1):
dapiFid_ch_HelpButton = tk.Button(acquisition_labelFrame, text="?",
                                  command=partial(display_help, "fiducialDAPI_channel"))
dapiFid_ch_HelpButton.grid(row=0, column=5)

# Barcode Channel Help Button (tab1):
barcode_ch_HelpButton = tk.Button(acquisition_labelFrame, text="?", command=partial(display_help, "barcode_channel"))
barcode_ch_HelpButton.grid(row=3, column=2)

# Barcode Fiducial Channel Help Button (tab1):
barcodeFid_ch_HelpButton = tk.Button(acquisition_labelFrame, text="?",
                                     command=partial(display_help, "fiducialBarcode_channel"))
barcodeFid_ch_HelpButton.grid(row=3, column=5)

# Mask Channel Help Button (tab1):
mask_ch_HelpButton = tk.Button(acquisition_labelFrame, text="?",
                               command=partial(display_help, "mask_channel"))
mask_ch_HelpButton.grid(row=5, column=2)

# Barcode Fiducial Channel Help Button (tab1):
maskFid_ch_HelpButton = tk.Button(acquisition_labelFrame, text="?",
                                  command=partial(display_help, "fiducialMask_channel"))
maskFid_ch_HelpButton.grid(row=5, column=5)

# RNA Channel Help Button (tab1):
rna_ch_HelpButton = tk.Button(acquisition_labelFrame, text="?", command=partial(display_help, "RNA_channel"))
rna_ch_HelpButton.grid(row=7, column=2)


# pixelSizeXY Help Button (tab1):
pixelSizeXY_HelpButton = tk.Button(acquisition_labelFrame, text="?", command=partial(display_help, "pixelSizeXY"))
pixelSizeXY_HelpButton.grid(row=10, column=2)

# pixelSizeZ Help Button (tab1):
pixelSizeZ_HelpButton = tk.Button(acquisition_labelFrame, text="?", command=partial(display_help, "pixelSizeZ"))
pixelSizeZ_HelpButton.grid(row=10, column=5)

# ReferenceFiducial Help Button (tab1):
referenceFiducial_HelpButton = tk.Button(alignImages_LabelFrame, text="?",
                                         command=partial(display_help, "referenceFiducial"))
referenceFiducial_HelpButton.grid(row=12, column=2)

# BlockSize Help Button (tab2):
blockSize_HelpButton = tk.Button(alignImages_LabelFrame2, text="?", command=partial(display_help, "blockSize"))
blockSize_HelpButton.grid(row=0, column=2)

# flux_min Help Button (tab2):
flux_min_HelpButton = tk.Button(buildsPWDmatrix_LabelFrame, text="?", command=partial(display_help, "flux_min"))
flux_min_HelpButton.grid(row=1, column=2)

# flux_min_3D Help Button (tab2):
flux_min_3D_HelpButton = tk.Button(buildsPWDmatrix_LabelFrame, text="?", command=partial(display_help, "flux_min_3D"))
flux_min_3D_HelpButton.grid(row=1, column=5)

# toleranceDrift Help Button (tab2):
toleranceDrift_HelpButton = tk.Button(buildsPWDmatrix_LabelFrame, text="?",
                                      command=partial(display_help, "toleranceDrift"))
toleranceDrift_HelpButton.grid(row=2, column=2)

# mask_expansion Help Button (tab2):
mask_expansion_HelpButton = tk.Button(buildsPWDmatrix_LabelFrame, text="?",
                                      command=partial(display_help, "mask_expansion"))
mask_expansion_HelpButton.grid(row=2, column=5)

# folder Help Button (tab2):
mask_expansion_HelpButton = tk.Button(buildsPWDmatrix_LabelFrame, text="?", command=partial(display_help, "folder"))
mask_expansion_HelpButton.grid(row=3, column=2)

# masks2process Help Button (tab2):
masks2process_HelpButton = tk.Button(buildsPWDmatrix_LabelFrame, text="?",
                                     command=partial(display_help, "masks2process"))
masks2process_HelpButton.grid(row=4, column=2)

# tracing_method Help Button (tab2):
tracing_method_HelpButton = tk.Button(buildsPWDmatrix_LabelFrame, text="?",
                                      command=partial(display_help, "tracing_method"))
tracing_method_HelpButton.grid(row=5, column=2)

# KDtree_distance_threshold_mum Help Button (tab2):
tracing_method_HelpButton = tk.Button(buildsPWDmatrix_LabelFrame, text="?",
                                      command=partial(display_help, "KDtree_distance_threshold_mum"))
tracing_method_HelpButton.grid(row=6, column=2)

# stardist_basename Help Button (tab2):
stardist_HelpButton = tk.Button(segmentedObjects_LabelFrame, text="?",
                                command=partial(display_help, "Stardist_basename"))
stardist_HelpButton.grid(row=7, column=2)

# brightest Help Button (tab2):
brightest_HelpButton = tk.Button(segmentedObjects_LabelFrame, text="?", command=partial(display_help, "brightest"))
brightest_HelpButton.grid(row=8, column=2)

# segmentObject_Labels_aeramax Help Button (tab2):
segmentObject_Labels_aeraMax_HelpButton = tk.Button(labels_LabelFrame, text="?",
                                                    command=partial(display_help, "segmentObject_Labels_dapi_aeraMax"))
segmentObject_Labels_aeraMax_HelpButton.grid(row=10, column=2)
# segmentObject_Labels_aeramin Help Button (tab2):
segmentObject_Labels_aeraMin_HelpButton = tk.Button(labels_LabelFrame, text="?",
                                                    command=partial(display_help, "segmentObject_Labels_dapi_aeraMin"))
segmentObject_Labels_aeraMin_HelpButton.grid(row=11, column=2)

# ZProject_dapi_zmax Help Button (tab2):
zProject_zmax_HelpButton = tk.Button(labels_LabelFrame, text="?", command=partial(display_help, "zProject_dapi_zmax"))
zProject_zmax_HelpButton.grid(row=10, column=5)
# ZProject_dapi_zmin Help Button (tab2):
zProject_zmin_HelpButton = tk.Button(labels_LabelFrame, text="?", command=partial(display_help, "zProject_dapi_zmin"))
zProject_zmin_HelpButton.grid(row=11, column=5)

# ZProject_Bcd_zmax Help Button (tab2):
zProject_Bcd_zmax_HelpButton = tk.Button(labels_LabelFrame, text="?",
                                         command=partial(display_help, "zProject_Bcd_zmax"))
zProject_Bcd_zmax_HelpButton.grid(row=13, column=2)
# ZProject_Bcd_zmin Help Button (tab2):
zProject_Bcd_zmin_HelpButton = tk.Button(labels_LabelFrame, text="?",
                                         command=partial(display_help, "zProject_Bcd_zmin"))
zProject_Bcd_zmin_HelpButton.grid(row=14, column=2)

# ZProject_Mask_zmax Help Button (tab2):
zProject_Mask_zmax_HelpButton = tk.Button(labels_LabelFrame, text="?",
                                          command=partial(display_help, "zProject_Mask_zmax"))
zProject_Mask_zmax_HelpButton.grid(row=13, column=5)
# ZProject_Mask_zmin Help Button (tab2):
zProject_Mask_zmin_HelpButton = tk.Button(labels_LabelFrame, text="?",
                                          command=partial(display_help, "zProject_Mask_zmin"))
zProject_Mask_zmin_HelpButton.grid(row=14, column=5)


# --------------------------Label of different channels for acquisition_labelFrame---------------------
# Dapi Channel Label:
dapi_ch_label = tk.Label(acquisition_labelFrame, text="DAPI Channel :")
dapi_ch_label.grid(row=0, column=0)
# Dapi Fiducial Channel Label:
dapiFid_ch_label = tk.Label(acquisition_labelFrame, text="DAPI Fiducial Channel :")
dapiFid_ch_label.grid(row=0, column=3)
# Barcode Channel Label:
barcode_ch_label = tk.Label(acquisition_labelFrame, text="Barcode/RT Channel :")
barcode_ch_label.grid(row=3, column=0)
# Barcode Fiducial Channel Label:
barcodeFid_ch_label = tk.Label(acquisition_labelFrame, text="Barcode/RT Fiducial Channel :")
barcodeFid_ch_label.grid(row=3, column=3)
# Mask Channel Label:
mask_ch_label = tk.Label(acquisition_labelFrame, text="Mask Channel :")
mask_ch_label.grid(row=5, column=0)
# Barcode Fiducial Channel Label:
maskFid_ch_label = tk.Label(acquisition_labelFrame, text="Mask Fiducial Channel :")
maskFid_ch_label.grid(row=5, column=3)
# RNA Channel Label:
rna_ch_label = tk.Label(acquisition_labelFrame, text="RNA Channel :")
rna_ch_label.grid(row=7, column=0)
# # RNA Fiducial Channel Label:
# rnaFid_ch_label = tk.Label(acquisition_labelFrame, text="RNA Fiducial Channel :")
# rnaFid_ch_label.grid(row=6, column=4)

# pixelSizeXY Label:
pixelSizeXY_label = tk.Label(acquisition_labelFrame, text="Pixel size (XY) :")
pixelSizeXY_label.grid(row=10, column=0)

# pixelSizeZ Label:
pixelSizeZ_label = tk.Label(acquisition_labelFrame, text="Pixel size (Z) :")
pixelSizeZ_label.grid(row=10, column=3)

# --------------------------Entry for acquisition_labelFrame---------------------
# pixelSizeXY Entry
pixelSizeXY_Entry = tk.Entry(acquisition_labelFrame)
value = infolist_partial["pixelSizeXY_Entry"]
pixelSizeXY_Entry.insert(0, string=value)
entries_dic["pixelSizeXY_Entry"] = [pixelSizeXY_Entry, value, type(value)]
pixelSizeXY_Entry.grid(row=10, column=1)

# pixelSizeZ Entry
pixelSizeZ_Entry = tk.Entry(acquisition_labelFrame)
value = infolist_partial["pixelSizeZ_Entry"]
pixelSizeZ_Entry.insert(0, string=value)
entries_dic["pixelSizeZ_Entry"] = [pixelSizeZ_Entry, value, type(value)]
pixelSizeZ_Entry.grid(row=10, column=4)

# --------------------------Radiobutton of different channels for acquisition_labelFrame---------------------
# Dapi Channel Radiobutton:
dapi_ch = tk.StringVar()
value = infolist_partial["dapi_ch"]
dapi_ch.set(value)
entries_dic["dapi_ch"] = [dapi_ch, value, type(value)]
dapi_ch_Radiobutton1 = tk.Radiobutton(acquisition_labelFrame, text="Ch00", value="ch00", variable=dapi_ch)
dapi_ch_Radiobutton1.grid(row=0, column=1)
dapi_ch_Radiobutton2 = tk.Radiobutton(acquisition_labelFrame, text="Ch01", value="ch01", variable=dapi_ch)
dapi_ch_Radiobutton2.grid(row=1, column=1)
dapi_ch_Radiobutton3 = tk.Radiobutton(acquisition_labelFrame, text="Ch02", value="ch02", variable=dapi_ch)
dapi_ch_Radiobutton3.grid(row=2, column=1)

# Dapi Fiducial Channel Radiobutton:
dapiFid_ch = tk.StringVar()
value = infolist_partial["dapiFid_ch"]
dapiFid_ch.set(value)
entries_dic["dapiFid_ch"] = [dapiFid_ch, value, type(value)]
dapiFid_ch_Radiobutton1 = tk.Radiobutton(acquisition_labelFrame, text="Ch00", value="ch00", variable=dapiFid_ch)
dapiFid_ch_Radiobutton1.grid(row=0, column=4)
dapiFid_ch_Radiobutton2 = tk.Radiobutton(acquisition_labelFrame, text="Ch01", value="ch01", variable=dapiFid_ch)
dapiFid_ch_Radiobutton2.grid(row=1, column=4)
dapiFid_ch_Radiobutton3 = tk.Radiobutton(acquisition_labelFrame, text="Ch02", value="ch02", variable=dapiFid_ch)
dapiFid_ch_Radiobutton3.grid(row=2, column=4)

# Barcode Channel Radiobutton:
barcode_ch = tk.StringVar()
value = infolist_partial["barcode_ch"]
barcode_ch.set(value)
entries_dic["barcode_ch"] = [barcode_ch, value, type(value)]
barcode_ch_Radiobutton1 = tk.Radiobutton(acquisition_labelFrame, text="Ch00", value="ch00", variable=barcode_ch)
barcode_ch_Radiobutton1.grid(row=3, column=1)
barcode_ch_Radiobutton2 = tk.Radiobutton(acquisition_labelFrame, text="Ch01", value="ch01", variable=barcode_ch)
barcode_ch_Radiobutton2.grid(row=4, column=1)

# Barcode Fiducial Channel Radiobutton:
barcodeFid_ch = tk.StringVar()
value = infolist_partial["barcodeFid_ch"]
barcodeFid_ch.set(value)
entries_dic["barcodeFid_ch"] = [barcodeFid_ch, value, type(value)]
barcodeFid_ch_Radiobutton1 = tk.Radiobutton(acquisition_labelFrame, text="Ch00", value="ch00", variable=barcodeFid_ch)
barcodeFid_ch_Radiobutton1.grid(row=3, column=4)
barcodeFid_ch_Radiobutton2 = tk.Radiobutton(acquisition_labelFrame, text="Ch01", value="ch01", variable=barcodeFid_ch)
barcodeFid_ch_Radiobutton2.grid(row=4, column=4)

# Mask Channel Radiobutton:
mask_ch = tk.StringVar()
value = infolist_partial["mask_ch"]
mask_ch.set(value)
entries_dic["mask_ch"] = [mask_ch, value, type(value)]
mask_ch_Radiobutton1 = tk.Radiobutton(acquisition_labelFrame, text="Ch00", value="ch00", variable=mask_ch)
mask_ch_Radiobutton1.grid(row=5, column=1)
mask_ch_Radiobutton2 = tk.Radiobutton(acquisition_labelFrame, text="Ch01", value="ch01", variable=mask_ch)
mask_ch_Radiobutton2.grid(row=6, column=1)

# # Barcode Fiducial Channel Radiobutton:
maskFid_ch = tk.StringVar()
value = infolist_partial["maskFid_ch"]
maskFid_ch.set(value)
entries_dic["maskFid_ch"] = [maskFid_ch, value, type(value)]
maskFid_ch_Radiobutton1 = tk.Radiobutton(acquisition_labelFrame, text="Ch00", value="ch00", variable=maskFid_ch)
maskFid_ch_Radiobutton1.grid(row=5, column=4)
maskFid_ch_Radiobutton2 = tk.Radiobutton(acquisition_labelFrame, text="Ch01", value="ch01", variable=maskFid_ch)
maskFid_ch_Radiobutton2.grid(row=6, column=4)

# # RNA Channel Radiobutton:
# rna_ch_Radiobutton
rna_ch = tk.StringVar()
value = infolist_partial["rna_ch"]
rna_ch.set(value)
entries_dic["rna_ch"] = [rna_ch, value, type(value)]
rna_ch_Radiobutton1 = tk.Radiobutton(acquisition_labelFrame, text="Ch00", value="ch00", variable=rna_ch)
rna_ch_Radiobutton1.grid(row=7, column=1)
rna_ch_Radiobutton2 = tk.Radiobutton(acquisition_labelFrame, text="Ch01", value="ch01", variable=rna_ch)
rna_ch_Radiobutton2.grid(row=8, column=1)
rna_ch_Radiobutton3 = tk.Radiobutton(acquisition_labelFrame, text="Ch02", value="ch02", variable=rna_ch)
rna_ch_Radiobutton3.grid(row=9, column=1)



# --------------------------------alignImages parameters Reference fiducial (tab1)-------------------------------

# ReferenceFiducial Label
referenceFiducial_label = tk.Label(alignImages_LabelFrame, text="Reference Fiducial :")
referenceFiducial_label.grid(row=12, column=0)

# ReferenceFiducial Entry
referenceFiducial_Entry = tk.Entry(alignImages_LabelFrame)
value = infolist_partial["referenceFiducial_Entry"]
referenceFiducial_Entry.insert(0, string=value)
entries_dic["referenceFiducial_Entry"] = [referenceFiducial_Entry, value, type(value)]
referenceFiducial_Entry.grid(row=12, column=1, pady=15)

# --------------------------------alignImages parameters Block Size (tab2)-------------------------------

# BlockSize Label
blockSize_label = tk.Label(alignImages_LabelFrame2, text="Block Size :")
blockSize_label.grid(row=0, column=0)

# BlockSize Entry
blockSize_Entry = tk.Entry(alignImages_LabelFrame2)
value = infolist_partial["blockSize_Entry"]
blockSize_Entry.insert(0, string=value)
entries_dic["blockSize_Entry"] = [blockSize_Entry, value, type(value)]
blockSize_Entry.grid(row=0, column=1, pady=15)

# -----------------------------buildsPWDmatrix parameters tab2)-------------------------------

# flux_min Label
flux_min_label = tk.Label(buildsPWDmatrix_LabelFrame, text="flux_min :")
flux_min_label.grid(row=1, column=0)

# flux_min_3D label
flux_min_3D_label = tk.Label(buildsPWDmatrix_LabelFrame, text="flux_min_3D :")
flux_min_3D_label.grid(row=1, column=3)

# toleranceDrift label
toleranceDrift_label = tk.Label(buildsPWDmatrix_LabelFrame, text="toleranceDrift :")
toleranceDrift_label.grid(row=2, column=0)

# mask_expansion label
mask_expansion_label = tk.Label(buildsPWDmatrix_LabelFrame, text="mask_expansion :")
mask_expansion_label.grid(row=2, column=3)

# folder label
folder_label = tk.Label(buildsPWDmatrix_LabelFrame, text="folder :")
folder_label.grid(row=3, column=0)

# masks2process label
masks2process_label = tk.Label(buildsPWDmatrix_LabelFrame, text="masks2process :")
masks2process_label.grid(row=4, column=0)

# tracing_method label
tracing_method_label = tk.Label(buildsPWDmatrix_LabelFrame, text="tracing_method :")
tracing_method_label.grid(row=5, column=0)

# KDtree_distance_threshold_mum label
KDtree_distance_threshold_mum_label = tk.Label(buildsPWDmatrix_LabelFrame, text="KDtree_distance_threshold_mum :")
KDtree_distance_threshold_mum_label.grid(row=6, column=0)

# flux_min Entry
flux_min_Entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
value = infolist_partial["flux_min_Entry"]
flux_min_Entry.insert(0, string=value)
entries_dic["flux_min_Entry"] = [flux_min_Entry, value, type(value)]
flux_min_Entry.grid(row=1, column=1)

# flux_min_3D Entry
flux_min_3D_Entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=10)
value = infolist_partial["flux_min_3D_Entry"]
flux_min_3D_Entry.insert(0, string=value)
entries_dic["flux_min_3D_Entry"] = [flux_min_3D_Entry, value, type(value)]
flux_min_3D_Entry.grid(row=1, column=4)

# toleranceDrift Entry
toleranceDrift_Entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
value = infolist_partial["toleranceDrift_Entry"]
toleranceDrift_Entry.insert(0, string=value)
entries_dic["toleranceDrift_Entry"] = [toleranceDrift_Entry, value, type(value)]
toleranceDrift_Entry.grid(row=2, column=1)


# mask_expansion Entry
mask_expansion_Entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=10)
value = infolist_partial["mask_expansion_Entry"]
mask_expansion_Entry.insert(0, string=value)
entries_dic["mask_expansion_Entry"] = [mask_expansion_Entry, value, type(value)]
mask_expansion_Entry.grid(row=2, column=4)

# folder Entry
folder_Entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
value = infolist_partial["folder_Entry"]
folder_Entry.insert(0, string=value)
entries_dic["folder_Entry"] = [folder_Entry, value, type(value)]
folder_Entry.grid(row=3, column=1)

# masks2process Entry
masks2process_Entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
value = infolist_partial["masks2process_Entry"]
masks2process_Entry.insert(0, string=value)
entries_dic["masks2process_Entry"] = [masks2process_Entry, value, type(value)]
masks2process_Entry.grid(row=4, column=1)

# tracing_method Entry
tracing_method_Entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
value = infolist_partial["tracing_method_Entry"]
tracing_method_Entry.insert(0, string=value)
entries_dic["tracing_method_Entry"] = [tracing_method_Entry, value, type(value)]
tracing_method_Entry.grid(row=5, column=1)

# KDtree_distance_threshold_mum Entry
KDtree_distance_threshold_mum_Entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
value = infolist_partial["KDtree_distance_threshold_mum_Entry"]
KDtree_distance_threshold_mum_Entry.insert(0, string=value)
entries_dic["KDtree_distance_threshold_mum_Entry"] = [KDtree_distance_threshold_mum_Entry, value, type(value)]
KDtree_distance_threshold_mum_Entry.grid(row=6, column=1)

# ----------------------SegmentedObjects parameters (tab2)-----------------------------------
# stardist_basename Label
stardist_label = tk.Label(segmentedObjects_LabelFrame, text="Stardist_basename :")
stardist_label.grid(row=7, column=0)

# brightest Label
brightest_label = tk.Label(segmentedObjects_LabelFrame, text="Brightest :")
brightest_label.grid(row=8, column=0)

# stardist_basename Entry
stardist_Entry = tk.Entry(segmentedObjects_LabelFrame, width=65)
value = infolist_partial["stardist_Entry"]
stardist_Entry.insert(0, string=value)
entries_dic["stardist_Entry"] = [stardist_Entry, value, type(value)]
stardist_Entry.grid(row=7, column=1)

# brightest Entry
brightest_Entry = tk.Entry(segmentedObjects_LabelFrame, width=25)
value = infolist_partial["brightest_Entry"]
brightest_Entry.insert(0, string=value)
entries_dic["brightest_Entry"] = [brightest_Entry, value, type(value)]
brightest_Entry.grid(row=8, column=1)


# -----------------------------------Labels parameters (tab2)-------------------------------------------

# segmentedObjects Aera Label
segmentedObjects_label = tk.Label(labels_LabelFrame, text="Dapi_SegmentedObjects :")
segmentedObjects_label.grid(row=9, column=0)
#       Aera_max Label
aeraMaxSegObjt_label = tk.Label(labels_LabelFrame, text="Aera_max :")
aeraMaxSegObjt_label.grid(row=10, column=0)
#       Aera_min Label
aeraMinSegObjt_label = tk.Label(labels_LabelFrame, text="Aera_min :")
aeraMinSegObjt_label.grid(row=11, column=0)

# zProject for Dapi Label
zProjectDapi_label = tk.Label(labels_LabelFrame, text="zProject for Dapi :")
zProjectDapi_label.grid(row=9, column=4)
#       zmax Dapi Label
zmaxZProjct_dapi_label = tk.Label(labels_LabelFrame, text="zmax :")
zmaxZProjct_dapi_label.grid(row=10, column=4)
#       zmin Dapi Label
zminZProjct_dapi_label = tk.Label(labels_LabelFrame, text="zmin :")
zminZProjct_dapi_label.grid(row=11, column=4)

# zProject barcode Label
zProjectBcd_label = tk.Label(labels_LabelFrame, text="zProject for Barcode :")
zProjectBcd_label.grid(row=12, column=0)
#       zmax Label
zmaxZProjct_Bcd_label = tk.Label(labels_LabelFrame, text="zmax :")
zmaxZProjct_Bcd_label.grid(row=13, column=0)
#       zmin Label
zminZProjct_Bcd_label = tk.Label(labels_LabelFrame, text="zmin :")
zminZProjct_Bcd_label.grid(row=14, column=0)

# zProject mask Label
zProjectMask_label = tk.Label(labels_LabelFrame, text="zProject for Mask :")
zProjectMask_label.grid(row=12, column=4)
#       zmax Label
zmaxZProjct_Mask_label = tk.Label(labels_LabelFrame, text="zmax :")
zmaxZProjct_Mask_label.grid(row=13, column=4)
#       zmin Label
zminZProjct_Mask_label = tk.Label(labels_LabelFrame, text="zmin :")
zminZProjct_Mask_label.grid(row=14, column=4)

# segmentedObjects Aera_max Entry
aeraMax_dapi_SegObjt_Entry = tk.Entry(labels_LabelFrame, width=10)
value = infolist_partial["aeraMax_dapi_SegObjt_Entry"]
aeraMax_dapi_SegObjt_Entry.insert(0, string=value)
entries_dic["aeraMax_dapi_SegObjt_Entry"] = [aeraMax_dapi_SegObjt_Entry, value, type(value)]
aeraMax_dapi_SegObjt_Entry.grid(row=10, column=1)
# segmentedObjects Aera_min Entry
aeraMin_dapi_SegObjt_Entry = tk.Entry(labels_LabelFrame, width=10)
value = infolist_partial["aeraMin_dapi_SegObjt_Entry"]
aeraMin_dapi_SegObjt_Entry.insert(0, string=value)
entries_dic["aeraMin_dapi_SegObjt_Entry"] = [aeraMin_dapi_SegObjt_Entry, value, type(value)]
aeraMin_dapi_SegObjt_Entry.grid(row=11, column=1)

# zProject zmax Entry
zProject_Dapi_zmax_Entry = tk.Entry(labels_LabelFrame, width=10)
value = infolist_partial["zProject_Dapi_zmax_Entry"]
zProject_Dapi_zmax_Entry.insert(0, string=value)
entries_dic["zProject_Dapi_zmax_Entry"] = [zProject_Dapi_zmax_Entry, value, type(value)]
zProject_Dapi_zmax_Entry.grid(row=10, column=4)

# zProject zmin Entry
zProject_Dapi_zmin_Entry = tk.Entry(labels_LabelFrame, width=10)
value = infolist_partial["zProject_Dapi_zmin_Entry"]
zProject_Dapi_zmin_Entry.insert(0, string=value)
entries_dic["zProject_Dapi_zmin_Entry"] = [zProject_Dapi_zmin_Entry, value, type(value)]
zProject_Dapi_zmin_Entry.grid(row=11, column=4)


# zProject barcode zmax Entry
zProject_Bcd_zmax_Entry = tk.Entry(labels_LabelFrame, width=10)
value = infolist_partial["zProject_Bcd_zmax_Entry"]
zProject_Bcd_zmax_Entry.insert(0, string=value)
entries_dic["zProject_Bcd_zmax_Entry"] = [zProject_Bcd_zmax_Entry, value, type(value)]
zProject_Bcd_zmax_Entry.grid(row=13, column=1)
# zProject barcode zmin Entry
zProject_Bcd_zmin_Entry = tk.Entry(labels_LabelFrame, width=10)
value = infolist_partial["zProject_Bcd_zmin_Entry"]
zProject_Bcd_zmin_Entry.insert(0, string=value)
entries_dic["zProject_Bcd_zmin_Entry"] = [zProject_Bcd_zmin_Entry, value, type(value)]
zProject_Bcd_zmin_Entry.grid(row=14, column=1)

# zProject mask zmax Entry
zProject_Mask_zmax_Entry = tk.Entry(labels_LabelFrame, width=10)
value = infolist_partial["zProject_Mask_zmax_Entry"]
zProject_Mask_zmax_Entry.insert(0, string=value)
entries_dic["zProject_Mask_zmax_Entry"] = [zProject_Mask_zmax_Entry, value, type(value)]
zProject_Mask_zmax_Entry.grid(row=13, column=4)
# zProject mask zmin Entry
zProject_Mask_zmin_Entry = tk.Entry(labels_LabelFrame, width=10)
value = infolist_partial["zProject_Mask_zmin_Entry"]
zProject_Mask_zmin_Entry.insert(0, string=value)
entries_dic["zProject_Mask_zmin_Entry"] = [zProject_Mask_zmin_Entry, value, type(value)]
zProject_Mask_zmin_Entry.grid(row=14, column=4)

main_window.mainloop()
