#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:23:36 2020

@author: marcnol

test fitting barcode spots to masks

"""

# =============================================================================
# IMPORTS
# =============================================================================


import glob, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from imageProcessing.imageProcessing import Image
from fileManagement import Folders
from fileManagement import Session, write_string_to_file

from astropy.table import Table, vstack, Column
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import sigma_clip, sigma_clipped_stats

from photutils.segmentation import SegmentationImage
from sklearn.metrics import pairwise_distances

# =============================================================================
# CLASSES
# =============================================================================


class CellID:
    def __init__(self, barcode_map_roi, masks, ROI):
        self.barcode_map_roi = barcode_map_roi
        self.masks = masks
        self.n_cells_assigned = 0
        self.n_cells_unassigned = 0
        self.n_barcodes_in_mask = 0

        self.segmentation_mask = SegmentationImage(self.masks)
        self.number_masks = self.segmentation_mask.nlabels
        self.ROI = ROI

        self.barcodes_in_mask = {}
        for mask in range(self.number_masks + 1):
            self.barcodes_in_mask["maskID_" + str(mask)] = []

    def visualize(self):

        imageBarcodes = np.zeros([2048, 2048])
        MasksBarcodes = masks
        R = []

        for i in range(len(self.barcode_map_roi.groups[0])):
            y_int = int(self.barcode_map_roi.groups[0]["xcentroid"][i])
            x_int = int(self.barcode_map_roi.groups[0]["ycentroid"][i])
            barcode_id = self.barcode_map_roi.groups[0]["Barcode #"][i]
            imageBarcodes[x_int][y_int] = barcode_id
            MasksBarcodes[x_int][y_int] += 20 * barcode_id
            R.append([y_int, x_int, barcode_id])

        # Shows results
        Ra = np.array(R)
        plt.imshow(masks, origin="lower", cmap="jet")
        plt.scatter(Ra[:, 0], Ra[:, 1], s=5, c=Ra[:, 2], alpha=0.5)

    def align_by_masking(self):
        # [ Assigns barcodes to masks and creates <n_barcodes_in_mask> ]
        n_barcodes_in_mask = np.zeros(self.number_masks + 2)
        print("ROI:{}".format(self.ROI))
        for i in range(len(self.barcode_map_roi.groups[0])):
            y_int = int(self.barcode_map_roi.groups[0]["xcentroid"][i])
            x_int = int(self.barcode_map_roi.groups[0]["ycentroid"][i])
            # barcode_id = self.barcode_map_roi.groups[ROI]['Barcode #'][i]
            mask_id = self.masks[x_int][y_int]
            self.barcode_map_roi["CellID #"][i] = mask_id
            if mask_id > 0:
                n_barcodes_in_mask[mask_id] += 1
                self.barcodes_in_mask["maskID_" + str(mask_id)].append(i)

        self.n_cells_assigned = np.count_nonzero(n_barcodes_in_mask > 0)
        self.n_cells_unassigned = self.number_masks - self.n_cells_assigned
        self.n_barcodes_in_mask = n_barcodes_in_mask

    def build_distance_matrix(self, mode="mean"):
        """
        

        
        """
        print("building distance matrix")
        barcode_map_roi = self.barcode_map_roi

        # [ builds SCdistanceTable ]
        barcode_map_roi_cell_id = barcode_map_roi.group_by("CellID #")  # ROI data sorted by cellID
        rois, cell_id, n_barcodes, barcode_ids, p = [], [], [], [], []

        for key, group in zip(barcode_map_roi_cell_id.groups.keys, barcode_map_roi_cell_id.groups):

            if key["CellID #"] > 1:  # excludes cellID 0 as this is background
                R = np.column_stack((np.array(group["xcentroid"].data), np.array(group["ycentroid"].data),))
                rois.append(group["ROI #"].data[0])
                cell_id.append(key["CellID #"])
                n_barcodes.append(len(group))
                barcode_ids.append(group["Barcode #"].data)
                p.append(pairwise_distances(R))
                # print("CellID #={}, n_barcodes={}".format(key['CellID #'],len(group)))

        SCdistanceTable = Table()  # [],names=('CellID', 'barcode1', 'barcode2', 'distances'))
        SCdistanceTable["ROI #"] = rois
        SCdistanceTable["CellID #"] = cell_id
        SCdistanceTable["nBarcodes"] = n_barcodes
        SCdistanceTable["Barcode #"] = barcode_ids
        SCdistanceTable["PWDmatrix"] = p

        self.SCdistanceTable = SCdistanceTable

        print("Cells with barcodes found: {}".format(len(SCdistanceTable)))

        # [ builds sc_matrix ]
        number_matrices = len(SCdistanceTable)  # z dimensions of sc_matrix
        unique_barcodes = np.unique(barcode_map_roi["Barcode #"].data)
        number_unique_barcodes = unique_barcodes.shape[0]  # number of unique Barcodes for xy dimensions of sc_matrix
        sc_matrix = np.zeros((number_unique_barcodes, number_unique_barcodes, number_matrices))
        sc_matrix[:] = np.NaN

        for i_cell, sc_pwd_item in zip(range(number_matrices), SCdistanceTable):
            barcodes_to_process = sc_pwd_item["Barcode #"]
            for barcode1, ibarcode1 in zip(barcodes_to_process, range(len(barcodes_to_process))):
                index_barcode_1 = np.nonzero(unique_barcodes == barcode1)[0][0]
                for barcode2, ibarcode2 in zip(barcodes_to_process, range(len(barcodes_to_process))):
                    index_barcode_2 = np.nonzero(unique_barcodes == barcode2)[0][0]
                    if barcode1 != barcode2:
                        newdistance = sc_pwd_item["PWDmatrix"][ibarcode1][ibarcode2]
                        if mode == "last":
                            sc_matrix[index_barcode_1][index_barcode_2][i_cell] = newdistance
                        elif mode == "mean":
                            sc_matrix[index_barcode_1][index_barcode_2][i_cell] = np.nanmean(
                                [newdistance, sc_matrix[index_barcode_1][index_barcode_2][i_cell],]
                            )
                        elif mode == "min":
                            sc_matrix[index_barcode_1][index_barcode_2][i_cell] = np.nanmin(
                                [newdistance, sc_matrix[index_barcode_1][index_barcode_2][i_cell],]
                            )

        self.sc_matrix = sc_matrix
        self.mean_sc_matrix = np.nanmean(sc_matrix, axis=2)
        self.unique_barcodes = unique_barcodes


# =============================================================================
# FUNCTIONS
# =============================================================================

# loads coordinate file

root_folder = "/home/marcnol/data/Experiment_20/Embryo_1"
fullFolder = root_folder + "/rawData/segmentedObjects/"
filename_barcode_coordinates = fullFolder + "segmentedObjects_barcode.dat"

"""
root_folder='/home/marcnol/data/Experiment_15/Embryo_006_ROI18'
fullFolder=root_folder +'/rawData/segmentedObjects/'
filename_barcode_coordinates =  fullFolder+'segmentedObjects_barcode.dat'
filename_roi_masks = fullFolder +'scan_001_DAPI_018_ROI_converted_decon_ch00_Masks.npy'

root_folder='/home/marcnol/Documents/Images/Embryo_debug_dataset'
fullFolder=root_folder +'/raw_images/segmentedObjects/'
filename_barcode_coordinates =  fullFolder+'segmentedObjects_barcode.dat'
filename_roi_masks = fullFolder +'scan_001_DAPI_017_ROI_converted_decon_ch00_Masks.npy'
"""

# Processes Tables
barcode_map = Table.read(filename_barcode_coordinates, format="ascii.ecsv")
barcode_map_roi = barcode_map.group_by("ROI #")

sc_matrix_collated = []
for ROI in range(len(barcode_map_roi.groups.keys)):
    n_roi = barcode_map_roi.groups.keys[ROI][0]  # need to iterate over the first index

    print("rois detected: {}".format(barcode_map_roi.groups.keys))

    barcode_map_single_roi = barcode_map.group_by("ROI #").groups[ROI]

    # finds file for masks
    files_folder = glob.glob(root_folder + "/rawData/" + "*.tif")

    files_to_process = [
        file
        for file in files_folder
        if file.split("_")[-1].split(".")[0] == "ch00"
        and "DAPI" in file.split("_")
        and int(os.path.basename(file).split("_")[3]) == n_roi
    ]

    if len(files_to_process) > 0:

        # loads masks
        filename_roi_masks = os.path.basename(files_to_process[0]).split(".")[0] + "_Masks.npy"
        masks = np.load(fullFolder + filename_roi_masks)

        # Assigns barcodes to masks for a given ROI
        cell_roi = CellID(barcode_map_single_roi, masks, ROI)

        cell_roi.align_by_masking()

        cell_roi.build_distance_matrix("min")  # mean min last

        print("ROI: {}, N cells assigned: {} out of {}".format(ROI, cell_roi.n_cells_assigned, cell_roi.number_masks))

        unique_barcodes = cell_roi.unique_barcodes

        if len(sc_matrix_collated) > 0:
            sc_matrix_collated = np.concatenate((sc_matrix_collated, cell_roi.sc_matrix), axis=2)
        else:
            sc_matrix_collated = cell_roi.sc_matrix
        del cell_roi

#%%
pixel_size = 0.1
mean_sc_matrix = pixel_size * np.nanmedian(sc_matrix_collated, axis=2)
print("rois detected: {}".format(barcode_map_roi.groups.keys))

fig = plt.figure()
pos = plt.imshow(mean_sc_matrix, cmap="seismic")
plt.xlabel("barcode #")
plt.ylabel("barcode #")
plt.title("PWD matrix" + " | n=" + str(sc_matrix_collated.shape[2]) + " | rois=" + str(len(barcode_map_roi.groups.keys)))
plt.xticks(np.arange(sc_matrix_collated.shape[0]), unique_barcodes)
plt.yticks(np.arange(sc_matrix_collated.shape[0]), unique_barcodes)
cbar = plt.colorbar(pos)
cbar.minorticks_on()
cbar.set_label("distance, um")
plt.clim(0, 1.4)

#%%

n_plots_x = 3
size_x, size_y = n_plots_x * 4, 4
fig, (ax1, ax2, ax3) = plt.subplots(figsize=(size_x, size_y), ncols=n_plots_x)

pos1 = ax1.hist(pixel_size * sc_matrix_collated[0, 1, :])
pos2 = ax2.hist(pixel_size * sc_matrix_collated[0, 2, :])
pos3 = ax3.hist(pixel_size * sc_matrix_collated[1, 2, :])
plt.xlabel("distance, um")
plt.ylabel("counts")


#%%
"""
root_folder='/home/marcnol/Documents/Images/Embryo_debug_dataset'
fullFolder=root_folder +'/raw_images/segmentedObjects/'
filename_barcode_coordinates =  fullFolder+'segmentedObjects_barcode.dat'
filename_roi_masks = fullFolder +'scan_001_DAPI_017_ROI_converted_decon_ch00_Masks.npy'

barcode_map = Table.read(filename_barcode_coordinates,format='ascii.ecsv')
masks=np.load(filename_roi_masks)
# loads masks

# Processes Tables 
ROI=0
barcode_map_roi=barcode_map.group_by('ROI #').groups[ROI]

print('rois detected: {}'.format(barcode_map_roi.groups.keys))

# Assigns barcodes to Masks for a given ROI
cellROI2 = cellID(barcode_map_roi,Masks)

cellROI2.align_by_masking()
print('N cells assigned: {} out of {}'.format(cell_roi.n_cells_assigned,cell_roi.number_masks))

cellROI2.build_distance_matrix('min') # mean min last

# collates SC matrices
sc_matrix_collated=np.concatenate((cell_roi.sc_matrix,cellROI2.sc_matrix),axis=2)
"""
