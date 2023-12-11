#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import tkinter as tk
from functools import partial
from tkinter import END, messagebox, ttk

from toolbox.parameter_file.function_parameters import (
    check_settings,
    import_parameters,
    update_parameters,
)
from toolbox.parameter_file.info_parameters import help_dic


def main():
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
        default_entry_dic, _ = import_parameters(script_dir)
        dapi_ch.set(default_entry_dic["dapi_ch"])
        dapiFid_ch.set(default_entry_dic["dapiFid_ch"])
        barcode_ch.set(default_entry_dic["barcode_ch"])
        barcodeFid_ch.set(default_entry_dic["barcodeFid_ch"])
        mask_ch.set(default_entry_dic["mask_ch"])
        maskFid_ch.set(default_entry_dic["maskFid_ch"])
        rna_ch.set(default_entry_dic["rna_ch"])
        # rnaFid_ch.set(default_entry_dic["rnaFid_ch"])
        for key, list_values in entries_dic.items():
            if (
                list_values[0].get() != "ch00"
                and list_values[0].get() != "ch01"
                and list_values[0].get() != "ch02"
            ):
                list_values[0].delete(first=0, last=END)
                list_values[0].insert(0, string=default_entry_dic[key])
                # list_values[0].insert(0, string=list_values[1])

    def save_setting(entries_dic, user_values_dic):
        if check_settings(entries_dic):
            # update user values in entries_dic:
            for key, list_values in entries_dic.items():
                entered_value = list_values[0].get()
                user_values_dic[key] = entered_value
            update_parameters(user_values_dic, parameters_full, current_dir)
            messagebox.showinfo("Info Message", "The settings have been saved.")

    # --------------------------------import values parameters from a parameters.json-----------------------------------------
    parameters_partial, parameters_full = import_parameters(script_dir, current_dir)

    # -----------------------------------------Save and restore button------------------------------------------------------

    # Restore button
    restore_button = tk.Button(
        main_window, text="Restore default settings", command=restore_setting
    )
    restore_button.grid(row=27, column=0, columnspan=3, pady=5, padx=5)

    # Save button
    save_button = tk.Button(
        main_window,
        text="Save settings",
        command=partial(save_setting, entries_dic, user_values_dic),
    )
    save_button.grid(row=27, column=2, columnspan=3, pady=5, padx=5)

    # --------------------------------------------Tab creation----------------------------------------------------
    notebook = ttk.Notebook(main_window)
    notebook.grid(row=0, column=0, rowspan=11, columnspan=5, pady=10)

    tab1 = ttk.Frame(notebook, width=700, height=700)
    tab1.grid(row=0, column=0, rowspan=10, columnspan=5)

    tab2 = ttk.Frame(notebook, width=700, height=700)
    tab2.grid(row=0, column=0, rowspan=10, columnspan=5)

    notebook.add(tab1, text="Standard Settings")
    notebook.add(tab2, text="Expert Settings")

    # -------------------------Label Frame--------------------------------
    # label frame (box) for quick help box
    help_LabelFrame = tk.LabelFrame(main_window, text="Quick Help")
    help_LabelFrame.grid(row=0, column=6, pady=5, padx=5)

    # label frame (box) for Acquisition parameters (tab1)
    acquisition_labelFrame = tk.LabelFrame(tab1, text="1. Acquisition parameters")
    acquisition_labelFrame.grid(row=0, column=0, columnspan=6, pady=5, padx=5)

    # label frame (box) for AlignImages parameters (tab1)
    alignImages_LabelFrame = tk.LabelFrame(tab1, text="2. AlignImages parameters")
    alignImages_LabelFrame.grid(
        row=11, column=0, columnspan=6, pady=5, padx=5, sticky="WE"
    )

    # label frame (box) for AlignImages parameters (tab2)
    alignImages_LabelFrame2 = tk.LabelFrame(tab2, text="1. AlignImages parameters")
    alignImages_LabelFrame2.grid(
        row=0, column=0, columnspan=6, pady=5, padx=5, sticky="WE"
    )

    # label frame (box) for BuildsPWDmatrix parameters (tab2)
    buildsPWDmatrix_LabelFrame = tk.LabelFrame(
        tab2, text="2. BuildsPWDmatrix parameters"
    )
    buildsPWDmatrix_LabelFrame.grid(
        row=1, column=0, columnspan=6, pady=5, padx=5, sticky="WE"
    )
    #
    # label frame (box) for SegmentedObjects parameters (tab2)
    segmentedObjects_LabelFrame = tk.LabelFrame(
        tab2, text="3. SegmentedObjects parameters"
    )
    segmentedObjects_LabelFrame.grid(
        row=7, column=0, columnspan=6, pady=5, padx=5, sticky="WE"
    )

    # label frame (box) for Labels parameters (tab2)
    labels_LabelFrame = tk.LabelFrame(tab2, text="4. Labels parameters")
    labels_LabelFrame.grid(row=9, column=0, columnspan=6, pady=5, padx=5, sticky="WE")

    # ---------------------------Help Box-----------------------------------
    help_label = tk.Label(
        help_LabelFrame,
        text="Help will come to those who ask for it",
        height=20,
        width=50,
        justify="center",
        wraplength=350,
    )
    help_label.grid(row=0, column=6, rowspan=11)

    # ---------------------------Image pyHiM----------------------------------------------------------
    img_path = script_dir + os.sep + "pyhim_image.png"
    img = tk.PhotoImage(file=img_path)
    img = img.zoom(10)
    img = img.subsample(15)
    label_img = tk.Label(main_window, image=img)
    label_img.grid(row=10, column=6, pady=10)

    # -----------------------------------------Help Button--------------------------------------------
    # Dapi Channel Help Button (tab1):
    dapi_ch_help_button = tk.Button(
        acquisition_labelFrame, text="?", command=partial(display_help, "DAPI_channel")
    )
    dapi_ch_help_button.grid(row=0, column=2)

    # Dapi Fiducial Channel Help Button (tab1):
    dapiFid_ch_help_button = tk.Button(
        acquisition_labelFrame,
        text="?",
        command=partial(display_help, "fiducialDAPI_channel"),
    )
    dapiFid_ch_help_button.grid(row=0, column=5)

    # Barcode Channel Help Button (tab1):
    barcode_ch_help_button = tk.Button(
        acquisition_labelFrame,
        text="?",
        command=partial(display_help, "barcode_channel"),
    )
    barcode_ch_help_button.grid(row=3, column=2)

    # Barcode Fiducial Channel Help Button (tab1):
    barcodeFid_ch_help_button = tk.Button(
        acquisition_labelFrame,
        text="?",
        command=partial(display_help, "fiducialBarcode_channel"),
    )
    barcodeFid_ch_help_button.grid(row=3, column=5)

    # Mask Channel Help Button (tab1):
    mask_ch_help_button = tk.Button(
        acquisition_labelFrame, text="?", command=partial(display_help, "mask_channel")
    )
    mask_ch_help_button.grid(row=5, column=2)

    # Barcode Fiducial Channel Help Button (tab1):
    maskFid_ch_help_button = tk.Button(
        acquisition_labelFrame,
        text="?",
        command=partial(display_help, "fiducialMask_channel"),
    )
    maskFid_ch_help_button.grid(row=5, column=5)

    # RNA Channel Help Button (tab1):
    rna_ch_help_button = tk.Button(
        acquisition_labelFrame, text="?", command=partial(display_help, "RNA_channel")
    )
    rna_ch_help_button.grid(row=7, column=2)

    # pixelSizeXY Help Button (tab1):
    pixelSizeXY_help_button = tk.Button(
        acquisition_labelFrame, text="?", command=partial(display_help, "pixelSizeXY")
    )
    pixelSizeXY_help_button.grid(row=10, column=2)

    # pixelSizeZ Help Button (tab1):
    pixelSizeZ_help_button = tk.Button(
        acquisition_labelFrame, text="?", command=partial(display_help, "pixelSizeZ")
    )
    pixelSizeZ_help_button.grid(row=10, column=5)

    # ReferenceFiducial Help Button (tab1):
    referenceFiducial_help_button = tk.Button(
        alignImages_LabelFrame,
        text="?",
        command=partial(display_help, "referenceFiducial"),
    )
    referenceFiducial_help_button.grid(row=12, column=2)

    # BlockSize Help Button (tab2):
    blockSize_help_button = tk.Button(
        alignImages_LabelFrame2, text="?", command=partial(display_help, "blockSize")
    )
    blockSize_help_button.grid(row=0, column=2)

    # flux_min Help Button (tab2):
    flux_min_help_button = tk.Button(
        buildsPWDmatrix_LabelFrame, text="?", command=partial(display_help, "flux_min")
    )
    flux_min_help_button.grid(row=1, column=2)

    # flux_min_3D Help Button (tab2):
    flux_min_3D_help_button = tk.Button(
        buildsPWDmatrix_LabelFrame,
        text="?",
        command=partial(display_help, "flux_min_3D"),
    )
    flux_min_3D_help_button.grid(row=1, column=5)

    # toleranceDrift Help Button (tab2):
    toleranceDrift_help_button = tk.Button(
        buildsPWDmatrix_LabelFrame,
        text="?",
        command=partial(display_help, "toleranceDrift"),
    )
    toleranceDrift_help_button.grid(row=2, column=2)

    # mask_expansion Help Button (tab2):
    mask_expansion_help_button = tk.Button(
        buildsPWDmatrix_LabelFrame,
        text="?",
        command=partial(display_help, "mask_expansion"),
    )
    mask_expansion_help_button.grid(row=2, column=5)

    # folder Help Button (tab2):
    mask_expansion_help_button = tk.Button(
        buildsPWDmatrix_LabelFrame, text="?", command=partial(display_help, "folder")
    )
    mask_expansion_help_button.grid(row=3, column=2)

    # masks2process Help Button (tab2):
    masks2process_help_button = tk.Button(
        buildsPWDmatrix_LabelFrame,
        text="?",
        command=partial(display_help, "masks2process"),
    )
    masks2process_help_button.grid(row=4, column=2)

    # tracing_method Help Button (tab2):
    tracing_method_help_button = tk.Button(
        buildsPWDmatrix_LabelFrame,
        text="?",
        command=partial(display_help, "tracing_method"),
    )
    tracing_method_help_button.grid(row=5, column=2)

    # KDtree_distance_threshold_mum Help Button (tab2):
    tracing_method_help_button = tk.Button(
        buildsPWDmatrix_LabelFrame,
        text="?",
        command=partial(display_help, "KDtree_distance_threshold_mum"),
    )
    tracing_method_help_button.grid(row=6, column=2)

    # stardist_basename Help Button (tab2):
    stardist_help_button = tk.Button(
        segmentedObjects_LabelFrame,
        text="?",
        command=partial(display_help, "Stardist_basename"),
    )
    stardist_help_button.grid(row=7, column=2)

    # brightest Help Button (tab2):
    brightest_help_button = tk.Button(
        segmentedObjects_LabelFrame,
        text="?",
        command=partial(display_help, "brightest"),
    )
    brightest_help_button.grid(row=8, column=2)

    # segmentObject_Labels_aeramax Help Button (tab2):
    segmentObject_Labels_aeraMax_help_button = tk.Button(
        labels_LabelFrame,
        text="?",
        command=partial(display_help, "segmentObject_Labels_dapi_aeraMax"),
    )
    segmentObject_Labels_aeraMax_help_button.grid(row=10, column=2)
    # segmentObject_Labels_aeramin Help Button (tab2):
    segmentObject_Labels_aeraMin_help_button = tk.Button(
        labels_LabelFrame,
        text="?",
        command=partial(display_help, "segmentObject_Labels_dapi_aeraMin"),
    )
    segmentObject_Labels_aeraMin_help_button.grid(row=11, column=2)

    # ZProject_dapi_zmax Help Button (tab2):
    zProject_zmax_help_button = tk.Button(
        labels_LabelFrame, text="?", command=partial(display_help, "zProject_dapi_zmax")
    )
    zProject_zmax_help_button.grid(row=10, column=5)
    # ZProject_dapi_zmin Help Button (tab2):
    zProject_zmin_help_button = tk.Button(
        labels_LabelFrame, text="?", command=partial(display_help, "zProject_dapi_zmin")
    )
    zProject_zmin_help_button.grid(row=11, column=5)

    # ZProject_Bcd_zmax Help Button (tab2):
    zProject_Bcd_zmax_help_button = tk.Button(
        labels_LabelFrame, text="?", command=partial(display_help, "zProject_Bcd_zmax")
    )
    zProject_Bcd_zmax_help_button.grid(row=13, column=2)
    # ZProject_Bcd_zmin Help Button (tab2):
    zProject_Bcd_zmin_help_button = tk.Button(
        labels_LabelFrame, text="?", command=partial(display_help, "zProject_Bcd_zmin")
    )
    zProject_Bcd_zmin_help_button.grid(row=14, column=2)

    # ZProject_Mask_zmax Help Button (tab2):
    zProject_Mask_zmax_help_button = tk.Button(
        labels_LabelFrame, text="?", command=partial(display_help, "zProject_Mask_zmax")
    )
    zProject_Mask_zmax_help_button.grid(row=13, column=5)
    # ZProject_Mask_zmin Help Button (tab2):
    zProject_Mask_zmin_help_button = tk.Button(
        labels_LabelFrame, text="?", command=partial(display_help, "zProject_Mask_zmin")
    )
    zProject_Mask_zmin_help_button.grid(row=14, column=5)

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
    barcodeFid_ch_label = tk.Label(
        acquisition_labelFrame, text="Barcode/RT Fiducial Channel :"
    )
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
    pixelSizeXY_entry = tk.Entry(acquisition_labelFrame)
    value = parameters_partial["pixelSizeXY_entry"]
    pixelSizeXY_entry.insert(0, string=value)
    entries_dic["pixelSizeXY_entry"] = [pixelSizeXY_entry, value, type(value)]
    pixelSizeXY_entry.grid(row=10, column=1)

    # pixelSizeZ Entry
    pixelSizeZ_entry = tk.Entry(acquisition_labelFrame)
    value = parameters_partial["pixelSizeZ_entry"]
    pixelSizeZ_entry.insert(0, string=value)
    entries_dic["pixelSizeZ_entry"] = [pixelSizeZ_entry, value, type(value)]
    pixelSizeZ_entry.grid(row=10, column=4)

    # --------------------------Radiobutton of different channels for acquisition_labelFrame---------------------
    # Dapi Channel Radiobutton:
    dapi_ch = tk.StringVar()
    value = parameters_partial["dapi_ch"]
    dapi_ch.set(value)
    entries_dic["dapi_ch"] = [dapi_ch, value, type(value)]
    dapi_ch_Radiobutton1 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch00", value="ch00", variable=dapi_ch
    )
    dapi_ch_Radiobutton1.grid(row=0, column=1)
    dapi_ch_Radiobutton2 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch01", value="ch01", variable=dapi_ch
    )
    dapi_ch_Radiobutton2.grid(row=1, column=1)
    dapi_ch_Radiobutton3 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch02", value="ch02", variable=dapi_ch
    )
    dapi_ch_Radiobutton3.grid(row=2, column=1)

    # Dapi Fiducial Channel Radiobutton:
    dapiFid_ch = tk.StringVar()
    value = parameters_partial["dapiFid_ch"]
    dapiFid_ch.set(value)
    entries_dic["dapiFid_ch"] = [dapiFid_ch, value, type(value)]
    dapiFid_ch_Radiobutton1 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch00", value="ch00", variable=dapiFid_ch
    )
    dapiFid_ch_Radiobutton1.grid(row=0, column=4)
    dapiFid_ch_Radiobutton2 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch01", value="ch01", variable=dapiFid_ch
    )
    dapiFid_ch_Radiobutton2.grid(row=1, column=4)
    dapiFid_ch_Radiobutton3 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch02", value="ch02", variable=dapiFid_ch
    )
    dapiFid_ch_Radiobutton3.grid(row=2, column=4)

    # Barcode Channel Radiobutton:
    barcode_ch = tk.StringVar()
    value = parameters_partial["barcode_ch"]
    barcode_ch.set(value)
    entries_dic["barcode_ch"] = [barcode_ch, value, type(value)]
    barcode_ch_Radiobutton1 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch00", value="ch00", variable=barcode_ch
    )
    barcode_ch_Radiobutton1.grid(row=3, column=1)
    barcode_ch_Radiobutton2 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch01", value="ch01", variable=barcode_ch
    )
    barcode_ch_Radiobutton2.grid(row=4, column=1)

    # Barcode Fiducial Channel Radiobutton:
    barcodeFid_ch = tk.StringVar()
    value = parameters_partial["barcodeFid_ch"]
    barcodeFid_ch.set(value)
    entries_dic["barcodeFid_ch"] = [barcodeFid_ch, value, type(value)]
    barcodeFid_ch_Radiobutton1 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch00", value="ch00", variable=barcodeFid_ch
    )
    barcodeFid_ch_Radiobutton1.grid(row=3, column=4)
    barcodeFid_ch_Radiobutton2 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch01", value="ch01", variable=barcodeFid_ch
    )
    barcodeFid_ch_Radiobutton2.grid(row=4, column=4)

    # Mask Channel Radiobutton:
    mask_ch = tk.StringVar()
    value = parameters_partial["mask_ch"]
    mask_ch.set(value)
    entries_dic["mask_ch"] = [mask_ch, value, type(value)]
    mask_ch_Radiobutton1 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch00", value="ch00", variable=mask_ch
    )
    mask_ch_Radiobutton1.grid(row=5, column=1)
    mask_ch_Radiobutton2 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch01", value="ch01", variable=mask_ch
    )
    mask_ch_Radiobutton2.grid(row=6, column=1)

    # # Barcode Fiducial Channel Radiobutton:
    maskFid_ch = tk.StringVar()
    value = parameters_partial["maskFid_ch"]
    maskFid_ch.set(value)
    entries_dic["maskFid_ch"] = [maskFid_ch, value, type(value)]
    maskFid_ch_Radiobutton1 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch00", value="ch00", variable=maskFid_ch
    )
    maskFid_ch_Radiobutton1.grid(row=5, column=4)
    maskFid_ch_Radiobutton2 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch01", value="ch01", variable=maskFid_ch
    )
    maskFid_ch_Radiobutton2.grid(row=6, column=4)

    # # RNA Channel Radiobutton:
    # rna_ch_Radiobutton
    rna_ch = tk.StringVar()
    value = parameters_partial["rna_ch"]
    rna_ch.set(value)
    entries_dic["rna_ch"] = [rna_ch, value, type(value)]
    rna_ch_Radiobutton1 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch00", value="ch00", variable=rna_ch
    )
    rna_ch_Radiobutton1.grid(row=7, column=1)
    rna_ch_Radiobutton2 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch01", value="ch01", variable=rna_ch
    )
    rna_ch_Radiobutton2.grid(row=8, column=1)
    rna_ch_Radiobutton3 = tk.Radiobutton(
        acquisition_labelFrame, text="Ch02", value="ch02", variable=rna_ch
    )
    rna_ch_Radiobutton3.grid(row=9, column=1)

    # --------------------------------alignImages parameters Reference fiducial (tab1)-------------------------------

    # ReferenceFiducial Label
    referenceFiducial_label = tk.Label(
        alignImages_LabelFrame, text="Reference Fiducial :"
    )
    referenceFiducial_label.grid(row=12, column=0)

    # ReferenceFiducial Entry
    referenceFiducial_entry = tk.Entry(alignImages_LabelFrame)
    value = parameters_partial["referenceFiducial_entry"]
    referenceFiducial_entry.insert(0, string=value)
    entries_dic["referenceFiducial_entry"] = [
        referenceFiducial_entry,
        value,
        type(value),
    ]
    referenceFiducial_entry.grid(row=12, column=1, pady=15)

    # --------------------------------alignImages parameters Block Size (tab2)-------------------------------

    # BlockSize Label
    blockSize_label = tk.Label(alignImages_LabelFrame2, text="Block Size :")
    blockSize_label.grid(row=0, column=0)

    # BlockSize Entry
    blockSize_entry = tk.Entry(alignImages_LabelFrame2)
    value = parameters_partial["blockSize_entry"]
    blockSize_entry.insert(0, string=value)
    entries_dic["blockSize_entry"] = [blockSize_entry, value, type(value)]
    blockSize_entry.grid(row=0, column=1, pady=15)

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
    KDtree_distance_threshold_mum_label = tk.Label(
        buildsPWDmatrix_LabelFrame, text="KDtree_distance_threshold_mum :"
    )
    KDtree_distance_threshold_mum_label.grid(row=6, column=0)

    # flux_min Entry
    flux_min_entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
    value = parameters_partial["flux_min_entry"]
    flux_min_entry.insert(0, string=value)
    entries_dic["flux_min_entry"] = [flux_min_entry, value, type(value)]
    flux_min_entry.grid(row=1, column=1)

    # flux_min_3D Entry
    flux_min_3D_entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=10)
    value = parameters_partial["flux_min_3D_entry"]
    flux_min_3D_entry.insert(0, string=value)
    entries_dic["flux_min_3D_entry"] = [flux_min_3D_entry, value, type(value)]
    flux_min_3D_entry.grid(row=1, column=4)

    # toleranceDrift Entry
    toleranceDrift_entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
    value = parameters_partial["toleranceDrift_entry"]
    toleranceDrift_entry.insert(0, string=value)
    entries_dic["toleranceDrift_entry"] = [toleranceDrift_entry, value, type(value)]
    toleranceDrift_entry.grid(row=2, column=1)

    # mask_expansion Entry
    mask_expansion_entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=10)
    value = parameters_partial["mask_expansion_entry"]
    mask_expansion_entry.insert(0, string=value)
    entries_dic["mask_expansion_entry"] = [mask_expansion_entry, value, type(value)]
    mask_expansion_entry.grid(row=2, column=4)

    # folder Entry
    folder_entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
    value = parameters_partial["folder_entry"]
    folder_entry.insert(0, string=value)
    entries_dic["folder_entry"] = [folder_entry, value, type(value)]
    folder_entry.grid(row=3, column=1)

    # masks2process Entry
    masks2process_entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
    value = parameters_partial["masks2process_entry"]
    masks2process_entry.insert(0, string=value)
    entries_dic["masks2process_entry"] = [masks2process_entry, value, type(value)]
    masks2process_entry.grid(row=4, column=1)

    # tracing_method Entry
    tracing_method_entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
    value = parameters_partial["tracing_method_entry"]
    tracing_method_entry.insert(0, string=value)
    entries_dic["tracing_method_entry"] = [tracing_method_entry, value, type(value)]
    tracing_method_entry.grid(row=5, column=1)

    # KDtree_distance_threshold_mum Entry
    KDtree_distance_threshold_mum_entry = tk.Entry(buildsPWDmatrix_LabelFrame, width=25)
    value = parameters_partial["KDtree_distance_threshold_mum_entry"]
    KDtree_distance_threshold_mum_entry.insert(0, string=value)
    entries_dic["KDtree_distance_threshold_mum_entry"] = [
        KDtree_distance_threshold_mum_entry,
        value,
        type(value),
    ]
    KDtree_distance_threshold_mum_entry.grid(row=6, column=1)

    # ----------------------SegmentedObjects parameters (tab2)-----------------------------------
    # stardist_basename Label
    stardist_label = tk.Label(segmentedObjects_LabelFrame, text="Stardist_basename :")
    stardist_label.grid(row=7, column=0)

    # brightest Label
    brightest_label = tk.Label(segmentedObjects_LabelFrame, text="Brightest :")
    brightest_label.grid(row=8, column=0)

    # stardist_basename Entry
    stardist_entry = tk.Entry(segmentedObjects_LabelFrame, width=65)
    value = parameters_partial["stardist_entry"]
    stardist_entry.insert(0, string=value)
    entries_dic["stardist_entry"] = [stardist_entry, value, type(value)]
    stardist_entry.grid(row=7, column=1)

    # brightest Entry
    brightest_entry = tk.Entry(segmentedObjects_LabelFrame, width=25)
    value = parameters_partial["brightest_entry"]
    brightest_entry.insert(0, string=value)
    entries_dic["brightest_entry"] = [brightest_entry, value, type(value)]
    brightest_entry.grid(row=8, column=1)

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
    aeraMax_dapi_SegObjt_entry = tk.Entry(labels_LabelFrame, width=10)
    value = parameters_partial["aeraMax_dapi_SegObjt_entry"]
    aeraMax_dapi_SegObjt_entry.insert(0, string=value)
    entries_dic["aeraMax_dapi_SegObjt_entry"] = [
        aeraMax_dapi_SegObjt_entry,
        value,
        type(value),
    ]
    aeraMax_dapi_SegObjt_entry.grid(row=10, column=1)
    # segmentedObjects Aera_min Entry
    aeraMin_dapi_SegObjt_entry = tk.Entry(labels_LabelFrame, width=10)
    value = parameters_partial["aeraMin_dapi_SegObjt_entry"]
    aeraMin_dapi_SegObjt_entry.insert(0, string=value)
    entries_dic["aeraMin_dapi_SegObjt_entry"] = [
        aeraMin_dapi_SegObjt_entry,
        value,
        type(value),
    ]
    aeraMin_dapi_SegObjt_entry.grid(row=11, column=1)

    # zProject zmax Entry
    zProject_Dapi_zmax_entry = tk.Entry(labels_LabelFrame, width=10)
    value = parameters_partial["zProject_Dapi_zmax_entry"]
    zProject_Dapi_zmax_entry.insert(0, string=value)
    entries_dic["zProject_Dapi_zmax_entry"] = [
        zProject_Dapi_zmax_entry,
        value,
        type(value),
    ]
    zProject_Dapi_zmax_entry.grid(row=10, column=4)

    # zProject zmin Entry
    zProject_Dapi_zmin_entry = tk.Entry(labels_LabelFrame, width=10)
    value = parameters_partial["zProject_Dapi_zmin_entry"]
    zProject_Dapi_zmin_entry.insert(0, string=value)
    entries_dic["zProject_Dapi_zmin_entry"] = [
        zProject_Dapi_zmin_entry,
        value,
        type(value),
    ]
    zProject_Dapi_zmin_entry.grid(row=11, column=4)

    # zProject barcode zmax Entry
    zProject_Bcd_zmax_entry = tk.Entry(labels_LabelFrame, width=10)
    value = parameters_partial["zProject_Bcd_zmax_entry"]
    zProject_Bcd_zmax_entry.insert(0, string=value)
    entries_dic["zProject_Bcd_zmax_entry"] = [
        zProject_Bcd_zmax_entry,
        value,
        type(value),
    ]
    zProject_Bcd_zmax_entry.grid(row=13, column=1)
    # zProject barcode zmin Entry
    zProject_Bcd_zmin_entry = tk.Entry(labels_LabelFrame, width=10)
    value = parameters_partial["zProject_Bcd_zmin_entry"]
    zProject_Bcd_zmin_entry.insert(0, string=value)
    entries_dic["zProject_Bcd_zmin_entry"] = [
        zProject_Bcd_zmin_entry,
        value,
        type(value),
    ]
    zProject_Bcd_zmin_entry.grid(row=14, column=1)

    # zProject mask zmax Entry
    zProject_Mask_zmax_entry = tk.Entry(labels_LabelFrame, width=10)
    value = parameters_partial["zProject_Mask_zmax_entry"]
    zProject_Mask_zmax_entry.insert(0, string=value)
    entries_dic["zProject_Mask_zmax_entry"] = [
        zProject_Mask_zmax_entry,
        value,
        type(value),
    ]
    zProject_Mask_zmax_entry.grid(row=13, column=4)
    # zProject mask zmin Entry
    zProject_Mask_zmin_entry = tk.Entry(labels_LabelFrame, width=10)
    value = parameters_partial["zProject_Mask_zmin_entry"]
    zProject_Mask_zmin_entry.insert(0, string=value)
    entries_dic["zProject_Mask_zmin_entry"] = [
        zProject_Mask_zmin_entry,
        value,
        type(value),
    ]
    zProject_Mask_zmin_entry.grid(row=14, column=4)

    main_window.mainloop()


if __name__ == "__main__":
    main()
