#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 15:35:52 2020

@author: marcnol

"""


import os
import bigfish
import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.plot as plot
from data import check_input_data
import numpy as np
print("Big-FISH version: {0}".format(bigfish.__version__))

#%% do not run

path_input = "/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug/bigfish/input"
path_output = "/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug/output"

# check input images are loaded
check_input_data(path_input)



#%%



path_input = "/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug/input"
path_output = "/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug/output"

# check input images are loaded
# check_input_data(path_input)


#%% 

recipe = {
    "fov": "fov_1",
    "c": ["dapi", "smfish"],
    "opt": "experience_1",
    "ext": "tif",
    "pattern": "opt_c_fov.ext"}

image = stack.build_stack(recipe, input_folder=path_input) 
print("\r shape: {0}".format(image.shape))
print("\r dtype: {0}".format(image.dtype))

rna = image[0, 1, ...]
print("smfish channel")
print("\r shape: {0}".format(rna.shape))
print("\r dtype: {0}".format(rna.dtype))

#%%
# parameters
voxel_size_z = 250
voxel_size_yx = 105
psf_z = 350
psf_yx = 100

# sigma
sigma_z, sigma_yx, sigma_yx = detection.get_sigma(voxel_size_z, voxel_size_yx, psf_z, psf_yx)
print("standard deviation of the PSF (z axis): {:0.3f} pixels".format(sigma_z))
print("standard deviation of the PSF (yx axis): {:0.3f} pixels".format(sigma_yx))

#%%

threshold = 150
spots = detection.detect_spots(rna, threshold, voxel_size_z, voxel_size_yx, psf_z, psf_yx)
print("detected spots")
print("\r shape: {0}".format(spots.shape))
print("\r dtype: {0}".format(spots.dtype))


# sigma
sigma_z, sigma_yx, sigma_yx = detection.get_sigma(voxel_size_z, voxel_size_yx, psf_z, psf_yx)
sigma = (sigma_z, sigma_yx, sigma_yx)

# LoG filter
rna_log = stack.log_filter(rna, sigma)

# local maximum detection
mask = detection.local_maximum_detection(rna_log, sigma)

# thresholding
spots, _ = detection.spots_thresholding(rna_log, mask, threshold)
print("detected spots")
print("\r shape: {0}".format(spots.shape))
print("\r dtype: {0}".format(spots.dtype))

#%%

(radius_z, radius_yx, radius_yx) = detection.get_radius(voxel_size_z, voxel_size_yx, psf_z, psf_yx)
print("radius z axis: {0:0.3f}".format(radius_z))
print("radius yx axis: {0:0.3f}".format(radius_yx))

image_contrasted = stack.rescale(rna, channel_to_stretch=0)
image_contrasted = stack.maximum_projection(image_contrasted)
plot.plot_detection(image_contrasted, spots, radius=radius_yx, framesize=(10, 8), remove_frame=True)


#%% 
spots_post_decomposition, clusters, reference_spot = detection.decompose_cluster(
    rna, spots, 
    voxel_size_z, voxel_size_yx, psf_z, psf_yx)

print("detected spots before decomposition")
print("\r shape: {0}".format(spots.shape))
print("\r dtype: {0}".format(spots.dtype))
print("detected spots after decomposition")
print("\r shape: {0}".format(spots_post_decomposition.shape))
print("\r dtype: {0}".format(spots_post_decomposition.dtype))

#%%

# sigma
sigma = detection.get_sigma(voxel_size_z, voxel_size_yx, psf_z, psf_yx)
large_sigma = tuple([sigma_ * 5 for sigma_ in sigma])

# denoising
rna_denoised = stack.remove_background_gaussian(rna, large_sigma)

# reference spot
reference_spot = detection.build_reference_spot(
    rna_denoised,
    spots,
    voxel_size_z, voxel_size_yx, psf_z, psf_yx)
threshold_cluster = int(reference_spot.max())

# fit a gaussian function on the reference spot
sigma_z, sigma_yx, amplitude, background = detection.modelize_spot(
    reference_spot, voxel_size_z, voxel_size_yx, psf_z, psf_yx)

# detect potential cluster regions
cluster_regions, spots_out_cluster, cluster_size = detection.get_clustered_region(
    rna_denoised, spots, threshold_cluster)

# precompute gaussian function values
max_grid = max(200, cluster_size + 1)
precomputed_gaussian = detection.precompute_erf(
    voxel_size_z, voxel_size_yx, sigma_z, sigma_yx, max_grid)

# gaussian mixtures
spots_in_cluster, _ = detection.fit_gaussian_mixture(
    rna_denoised,
    cluster_regions,
    voxel_size_z,
    voxel_size_yx,
    sigma_z,
    sigma_yx,
    amplitude,
    background,
    precomputed_gaussian)

spots_post_decomposition = np.concatenate((spots_out_cluster, spots_in_cluster[:, :3]), axis=0)
print("detected spots before decomposition")
print("\r shape: {0}".format(spots.shape))
print("\r dtype: {0}".format(spots.dtype))
print("detected spots after decomposition")
print("\r shape: {0}".format(spots_post_decomposition.shape))
print("\r dtype: {0}".format(spots_post_decomposition.dtype))


plot.plot_detection(image_contrasted, spots_post_decomposition, radius=radius_yx, 
                    framesize=(10, 8), remove_frame=True)
plot.plot_reference_spot(reference_spot, rescale=True, remove_frame=True)

