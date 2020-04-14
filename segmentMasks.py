#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 14:59:43 2020

@author: marcnol
"""



import glob,os
import matplotlib.pylab as plt
import numpy as np
#import cv2
#import matplotlib.pyplot as plt
#from matplotlib import cm
#from skimage import io
from astropy.stats import sigma_clipped_stats,SigmaClip,gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
from astropy.visualization import SqrtStretch,simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.table import Table, vstack, Column
from photutils import DAOStarFinder,CircularAperture,find_peaks,detect_sources
from photutils import detect_threshold,deblend_sources
from photutils import background, Background2D, MedianBackground

from imageProcessing import Image,saveImage2Dcmd
from fileManagement import folders
from fileManagement import session,writeString2File


def showsImageSources(im,im1_bkg_substracted,log1,sources,outputFileName):
    # show results
    fig = plt.figure()
    fig.set_size_inches((30,30))
    
    positions = np.transpose((sources['xcentroid']+.5, sources['ycentroid']+.5)) # for some reason sources are always displays 1/2 px from center of spot
    apertures = CircularAperture(positions, r=4.)
    norm = simple_norm(im, 'sqrt', percent=99.9)
    plt.imshow(im1_bkg_substracted, cmap='Greys', origin='lower', norm=norm)
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.xlim(0, im.shape[1]-1)
    plt.ylim(0, im.shape[0]-1)
    plt.savefig(outputFileName+'_segmentedSources.png')
    plt.close()
    writeString2File(log1.fileNameMD,"{}\n ![]({})\n".format(os.path.basename(outputFileName),outputFileName+'_segmentedSources.png'),'a')
    
def showsImageMasks(im,log1,segm_deblend,outputFileName):

    norm = ImageNormalize(stretch=SqrtStretch())
    cmap = segm_deblend.make_cmap(random_state=12345)

    fig = plt.figure()
    fig.set_size_inches((30,30))
    plt.imshow(im, cmap='Greys_r', origin='lower', norm=norm)
    plt.imshow(segm_deblend , origin='lower', cmap=cmap)
    
    '''
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
    ax1.imshow(im, origin='lower', cmap='Greys_r', norm=norm)
    ax1.set_title('Data')
    ax2.imshow(segm_deblend , origin='lower', cmap=cmap)
    ax2.set_title('Segmentation Image')
    '''
    plt.savefig(outputFileName+'_segmentedMasks.png')
    plt.close()
    writeString2File(log1.fileNameMD,"{}\n ![]({})\n".format(os.path.basename(outputFileName),outputFileName+'_segmentedMasks.png'),'a')


def segmentSourceInhomogBackground(im,param,log1,outputFileName):
    threshold_over_std=param.param['segmentedObjects']['threshold_over_std']
    fwhm=param.param['segmentedObjects']['fwhm']
    brightest=1000 # keeps brightest sources

    # estimates inhomogeneous background
    sigma_clip = SigmaClip(sigma=3.)
    bkg_estimator = MedianBackground()
    bkg = Background2D(im, (50, 50), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    
    im1_bkg_substracted=im-bkg.background
    mean, median, std = sigma_clipped_stats(im1_bkg_substracted, sigma=3.0) 
    
    # estimates sources
    daofind = DAOStarFinder(fwhm=fwhm, 
                            threshold=threshold_over_std*std, 
                            brightest=brightest,
                            exclude_border=True)  
    sources = daofind(im1_bkg_substracted) 
    
    # show results
    showsImageSources(im,im1_bkg_substracted,log1,sources,outputFileName)

    return sources

def segmentSourceFlatBackground(im,param,log1,outputFileName):
    threshold_over_std=param.param['segmentedObjects']['threshold_over_std']
    fwhm=param.param['segmentedObjects']['fwhm']

    # removes background
    mean, median, std = sigma_clipped_stats(im, param.param['segmentedObjects']['background_sigma']) 
    im1_bkg_substracted = im - median
    
    # estimates sources
    daofind = DAOStarFinder(fwhm=fwhm, 
                            threshold=threshold_over_std*std, 
                            exclude_border=True)  
    sources = daofind(im - median)  
    
    # show results
    showsImageSources(im,im1_bkg_substracted,log1,sources,outputFileName)

    return sources

def segmentMaskInhomogBackground(im,param,log1,outputFileName):

    # removes background
    threshold = detect_threshold(im, nsigma=2.)
    sigma_clip = SigmaClip(sigma=param.param['segmentedObjects']['background_sigma'])
    
    bkg_estimator = MedianBackground()
    bkg = Background2D(im, (50, 50), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    threshold = bkg.background + (param.param['segmentedObjects']['threshold_over_std'] * bkg.background_rms)  # background-only error image, typically 1.0
    
    sigma = param.param['segmentedObjects']['fwhm'] * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()

    # estimates masks and deblends
    segm = detect_sources(im, threshold, npixels=5, filter_kernel=kernel)
    
    segm_deblend = deblend_sources(im, segm, npixels=param.param['segmentedObjects']['area_min'], # typically 50 for DAPI
                                   filter_kernel=kernel, nlevels=32,
                                   contrast=0.001)
    # show results
    showsImageMasks(im,log1,segm_deblend,outputFileName)
    
    return segm_deblend

def makesSegmentations(fileName,param,log1,session1,dataFolder):

    rootFileName = os.path.basename(fileName).split('.')[0]
    outputFileName = dataFolder.outputFolders['segmentedObjects']+os.sep+rootFileName
    fileName_2d_aligned = dataFolder.outputFolders['alignImages']+os.sep+rootFileName+'_2d_registered.npy'

    if (fileName in session1.data) and \
        param.param['segmentedObjects']['operation']=='overwrite' and \
        os.path.exists(fileName_2d_aligned): # file exists
        
        ROI = os.path.basename(fileName).split('_')[param.param['acquisition']['positionROIinformation']]
        label=param.param['acquisition']['label']

        # loading registered 2D projection
        Im = Image()
        Im.loadImage2D(fileName,log1,dataFolder.outputFolders['alignImages'],tag='_2d_registered')
        im=Im.data_2D
        log1.report("[{}] Loaded 2D registered file: {}".format(label,os.path.basename(fileName)))        
            
        if label=='barcode' and len([i for i in rootFileName.split('_') if 'RT' in i])>0:        
            if param.param['segmentedObjects']['background_method']=='flat':
                output=segmentSourceFlatBackground(im,param,log1,outputFileName)
            elif param.param['segmentedObjects']['background_method']=='inhomogeneous':
                output=segmentSourceInhomogBackground(im,param,log1,outputFileName)
            else:
                log1.report('segmentedObjects/background_method not specified in json file','ERROR')
                output=Table()            

            barcodeID=os.path.basename(fileName).split('_')[2].split('RT')[1]
            colROI=Column(int(ROI)*np.ones(len(output)),name='ROI #')
            colBarcode=Column(int(barcodeID)*np.ones(len(output)),name='Barcode #')
            output.add_column(colBarcode,index=0)
            output.add_column(colROI,index=0)
               
            for col in output.colnames:  
                output[col].info.format = '%.8g'  # for consistent table output

        elif label=='DAPI' and rootFileName.split('_')[2]=='DAPI':
            if param.param['segmentedObjects']['background_method']=='flat':
                output = segmentMaskInhomogBackground(im,param,log1,outputFileName)
            elif param.param['segmentedObjects']['background_method']=='inhomogeneous':
                output = segmentMaskInhomogBackground(im,param,log1,outputFileName)
            else:
                log1.report('segmentedObjects/background_method not specified in json file','ERROR')
                output=np.zeros(1)
            
            # saves output 2d zProjection as matrix
            Im.saveImage2D(log1,dataFolder.outputFolders['zProject'])
            saveImage2Dcmd(output,outputFileName+'_Masks',log1)
        else:
            output=[]
        del Im

        return output
    else:
        log1.report("2D aligned file does not exist:{}".format(fileName_2d_aligned),'Error')
        return []
    
def segmentMasks(param,log1,session1):
    sessionName='segmentMasks'
 
    # processes folders and files 
    dataFolder=folders(param.param['rootFolder'])
    dataFolder.setsFolders()
    log1.addSimpleText("\n===================={}:{}====================\n".format(sessionName,
                                                          param.param['acquisition']['label']))
    log1.report('folders read: {}'.format(len(dataFolder.listFolders)))
    writeString2File(log1.fileNameMD,"## {}: {}\n".format(sessionName,
                                                          param.param['acquisition']['label']),'a') 
    barcodesCoordinates=Table()
    
    for currentFolder in dataFolder.listFolders:
        #currentFolder=dataFolder.listFolders[0]
        filesFolder=glob.glob(currentFolder+os.sep+'*.tif')
        dataFolder.createsFolders(currentFolder,param)
    
        # generates lists of files to process    
        param.files2Process(filesFolder)
        log1.report('About to read {} files\n'.format(len(param.fileList2Process)))
        
        for fileName in param.fileList2Process:
            label=param.param['acquisition']['label']
            if label!='fiducial':
                output = makesSegmentations(fileName,param,log1,session1,dataFolder)
                if label=='barcode':
                    outputFile = dataFolder.outputFiles['segmentedObjects']+'_'+label+'.dat'
                    barcodesCoordinates = vstack([barcodesCoordinates,output])
                    barcodesCoordinates.write(outputFile,format='ascii.ecsv',overwrite=True)
                    log1.report('File {} written to file.'.format(outputFile),'info')
                session1.add(fileName,sessionName)


    
    



