#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 14:55:21 2022

@author: Olivier Messina
"""

import os                             
from IPython.display import Image     
from IPython.core.display import HTML 
import matplotlib.pyplot as plt       
import matplotlib.image as mpimg      
from matplotlib import rcParams 
from tifffile import TiffWriter 
import numpy as np
from glob import glob
import glob 
import shutil

def copy_localization_table(Input_folder):
    source = Input_folder+"/segmentedObjects/segmentedObjects_3D_barcode_original.dat"
    dest = Input_folder+"/segmentedObjects/segmentedObjects_3D_barcode.dat"
    shutil.copyfile(source,dest)

def plot_zprojection(Input_folder,RTs_references,titles,datatype):

    # Figure size in inches optional
    rcParams['figure.figsize'] = 15, 15

    # Enter in the folder contaning the zProjected output images
    os.chdir(Input_folder + 'zProject')

    if (datatype =='DAPI'):

        # Create list of images
        img_A = glob.glob('*DAPI*'+'*ch01*'+'.png')
        img_B = glob.glob('*DAPI*'+'*ch00*'+'.png')
        Concat_images = [img_A,img_B]

    if (datatype =='RT'):

        # Create list of images
        img_A = glob.glob('*'+RTs_references+'*'+'*ch00*'+'.png')
        img_B = glob.glob('*'+RTs_references+'*'+'*ch01*'+'.png')
        Concat_images = [img_A,img_B]
        
    # plots images
    fig, ax = plt.subplots(1,2)
    for x, file in enumerate(Concat_images):
        # Read images
        Image = mpimg.imread(Input_folder + 'zProject/' + file[0])[500:800,0:300]
        # Display images
        ax[x].set_title(titles[x])
        ax[x].imshow(Image)
	     
def plot_alignment(Input_folder,RTs_references,titles):

    # Figure size in inches optional
    rcParams['figure.figsize'] = 15 ,10
    
    # Create list of images
    jpgFilenamesList_RTs_alignement_Difference = glob.glob(Input_folder + 'alignImages/'+'*'+RTs_references+'*'+'_referenceDifference.png')
    jpgFilenamesList_RTs_alignement_overlay = glob.glob(Input_folder + 'alignImages/'+'*'+RTs_references+'*'+'_overlay'+'*'+'.png')

    img_A = mpimg.imread(jpgFilenamesList_RTs_alignement_Difference[0])[1000:2000,1000:2000] # Zoom in the region of interest
    img_B = mpimg.imread(jpgFilenamesList_RTs_alignement_Difference[0])[1000:2000,3500:4500] # Zoom in the region of interest
    img_C = mpimg.imread(jpgFilenamesList_RTs_alignement_overlay[0])[1000:2000,1000:2000]*5 # Zoom in the region of interest


    # Concatenates images
    Concat_images = [img_A,img_B,img_C]

    # Plot 
    fig, ax = plt.subplots(1,3)
    for x, file in enumerate(titles):
        # Read images
        Image = Concat_images[x]
        # Display images
        ax[x].set_title(titles[x])
        ax[x].imshow(Image)    
    
def show_plot(titles,imgs,files):
    if len(imgs)>1:
        fig, ax = plt.subplots(1,2)
        for title, img, axis in zip(titles,imgs, ax):
            axis.set_title(title)
            axis.imshow(img)
    else:
        plt.imshow(imgs[0])
        plt.axis('off')

def plot_segment_object(Input_folder,RTs_references,titles,datatype):

    # Figure size in inches optional
    rcParams['figure.figsize'] = 15 ,10

    if (datatype =='DAPI' or datatype =='RT'):
        if (datatype == 'RT'):
            # Create list of images
            files = glob.glob(Input_folder + 'segmentedObjects/'+'*'+RTs_references+'*'+'_3DimageNlocalizations.png')
            imgs = [mpimg.imread(files[0])[500:4500,500:4500]]
            
        if (datatype == 'DAPI'):
            # Create list of images
            files = glob.glob(Input_folder + 'segmentedObjects/'+'*'+'DAPI'+'*'+'_3Dmasks.png')
            imgs = [mpimg.imread(files[0])[1500:3600,500:4500]]

    elif (datatype =='TRACES'):

        # Create list of images
        files = glob.glob(Input_folder + 'buildsPWDmatrix/'+'*'+'KDtree'+'*'+'_XYZ'+'*'+'.png')
        imgs = [mpimg.imread(x) for x in files]
        
        # crops png to get full and zoomed images
        imgs = [x[500:5500,500:5500] for x in imgs]
        imgs.append(imgs[0][4000:5000,4000:5000])

    print(f"$ Will plot: {files}")

    # Plots images 
    show_plot(titles,imgs,files)
        
def plot_matrix(input_folder, data_type='proximity'):

    # figure size in inches optional
    rcParams['figure.figsize'] = 10 ,10
    
    # Create list of images
    if data_type == 'PWD':
        files = glob.glob(input_folder + 'buildsPWDmatrix/'+'*PWDmatrixKDE.png')
    elif data_type == 'proximity':
        files = glob.glob(input_folder + 'buildsPWDmatrix/'+'*HiMmatrix.png')
       
    titles = [data_type]
    
    img = mpimg.imread(files[0])

    show_plot(titles,[img],"-")


