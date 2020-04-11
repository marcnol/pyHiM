#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 09:22:18 2020

@author: marcnol
"""
import numpy as np
import matplotlib.pyplot as plt
import os,glob
from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft
from skimage import exposure
#from scipy.ndimage import fourier_shift
from scipy.ndimage import shift as shiftImage
from imageProcessing import Image, imageAdjust, save2imagesRGB, saveImage2Dcmd
from fileManagement import folders,writeString2File,saveJSON,loadJSON
from collections import defaultdict

 
def RT2fileName(fileList2Process,referenceBarcode,positionROIinformation=3):    
    fileNameReferenceList = []
    ROIList = {}
    
    for file in fileList2Process:
        if referenceBarcode in file.split('_'):
            fileNameReferenceList.append(file)
            ROIList[file]= os.path.basename(file).split('_')[positionROIinformation]
                    
    return fileNameReferenceList, ROIList

def displaysEqualizationHistograms(I_histogram,lower_threshold,outputFileName,log1,verbose=False):
        # hist1_before, hist1_after,hist2_before, hist2_after, , vebose=False, fileName='test'):
    fig= plt.figure(figsize=(6, 3))
    ax1 = plt.subplot(2, 2, 1)
    ax2 = plt.subplot(2, 2, 2)
    ax3 = plt.subplot(2, 2, 3)
    ax4 = plt.subplot(2, 2, 4)
    
    ax1.plot(I_histogram['Im1'][0][1],I_histogram['Im1'][0][0])
    ax2.plot(I_histogram['Im2'][0][1],I_histogram['Im2'][0][0])
    ax3.plot(I_histogram['Im1'][1][1],I_histogram['Im1'][1][0])
    ax4.plot(I_histogram['Im2'][1][1],I_histogram['Im2'][1][0])
    ax3.set_yscale('log')
    ax4.set_yscale('log')
    ax1.vlines(lower_threshold['Im1'],0,I_histogram['Im1'][0][0].max(),colors='r')
    ax2.vlines(lower_threshold['Im2'],0,I_histogram['Im2'][0][0].max(),colors='r')
    plt.savefig(outputFileName+'_intensityHist.png')
    writeString2File(log1.fileNameMD,"{}\n ![]({})\n".format(os.path.basename(outputFileName),outputFileName+'_intensityHist.png'),'a')

    if not verbose:
        plt.close(fig)

def showCCimage(image1_uncorrected,image2_uncorrected,outputFileName,shift,verbose=False):
    image_product = np.fft.fft2(image1_uncorrected) * np.fft.fft2(image2_uncorrected).conj()
    cc_image = _upsampled_dft(image_product, 150, 100, (shift*100)+75).conj()
    if verbose:
        plt.figure(figsize=(60, 30))
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
        RGB_falsecolor_image=np.dstack([image1_uncorrected,image2_uncorrected,np.zeros([2048,2048])])
        ax1.imshow(RGB_falsecolor_image,origin='lower',interpolation='nearest')
        ax1.set_axis_off()
        ax1.set_title("Super-imposed images")
    
        ax2.imshow(cc_image.real)
        ax2.set_axis_off()
        ax2.set_title("Supersampled XC sub-area")
    else:
        plt.figure(figsize=(30, 30))
        plt.imsave(outputFileName+'_CC.png',cc_image)
        
def align2Files(fileName,imReference, param,log1,session1,dataFolder,verbose):

    fileName1=imReference.fileName
    fileName2=fileName
       
    outputFileName=dataFolder.outputFolders['alignImages']+os.sep+os.path.basename(fileName2).split('.')[0]

    # loads image
    Im2 = Image()
    Im2.loadImage2D(fileName2,log1,dataFolder)
    
    # adjusts image levels
    image1_uncorrected=imReference.data_2D/imReference.data_2D.max()
    image2_uncorrected=Im2.data_2D/Im2.data_2D.max()
    image1_adjusted,hist1_before,hist1_after, lower_cutoff1, higher_cutoff1 = imageAdjust(image1_uncorrected,
                                                  log1,
                                                  outputFileName+'_ref',
                                                  lower_threshold=0.999,
                                                  higher_threshold=.9999999, 
                                                  display=verbose) 
    image2_adjusted,hist2_before,hist2_after, lower_cutoff2, higher_cutoff2 = imageAdjust(image2_uncorrected,
                                                  log1,
                                                  outputFileName,
                                                  lower_threshold=0.999,
                                                  higher_threshold=.9999999, 
                                                  display=verbose) 
    # displays intensity histograms
    lower_threshold={'Im1':lower_cutoff1, 'Im2':lower_cutoff2}
    I_histogram={'Im1':(hist1_before,hist1_after),'Im2':(hist2_before,hist2_after)} 
    displaysEqualizationHistograms(I_histogram,lower_threshold,outputFileName,log1,verbose)
    
    # calculates shift
    shift, error, diffphase = register_translation(image1_adjusted, 
                                                   image2_adjusted, 100)
    
    # corrects image
    # The shift corresponds to the pixel offset relative to the reference image
    image2_corrected = shiftImage(image2_adjusted,shift)
    image2_corrected_raw = shiftImage(Im2.data_2D,shift)
    image2_corrected = exposure.rescale_intensity(image2_corrected ,out_range=(0,1))
    
    log1.report(f"Detected subpixel offset (y, x): {shift} px\n")
    
   # saves uncrrected images to file
    save2imagesRGB(image1_uncorrected,image2_uncorrected,
                   outputFileName+'_overlay_uncorrected.png')
    writeString2File(log1.fileNameMD,"{}\n ![]({})\n".format(os.path.basename(outputFileName),outputFileName+'_overlay_uncorrected.png'),'a')

    
    # thresholds corrected images for better display and saves
    image1_corrected=image1_adjusted>0.1
    image2_corrected=image2_corrected>0.1
    save2imagesRGB(image1_corrected,image2_corrected,
                   outputFileName+'_overlay_corrected.png')

    alignmentOutput = dataFolder.outputFiles['alignImages']
    list2output = "{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(os.path.basename(fileName2),
                                              os.path.basename(fileName1),
                                              shift[0],
                                              shift[1],
                                              error,
                                              diffphase)
    writeString2File(alignmentOutput+'.bed', list2output,'a')

    # saves raw registered image of fiducial
    saveImage2Dcmd(image2_corrected_raw,outputFileName+'_2d_registered',log1)
    
    del Im2
    return shift

def alignImages(param,log1,session1):

    if param.param['alignImages']['operation']=='overwrite':
        
        verbose=False
        
        # session
        sessionName='alignImages'
    
        # processes folders and files 
        dataFolder=folders(param.param['rootFolder'])
        dataFolder.setsFolders()
        log1.addSimpleText("\n===================={}====================\n".format(sessionName))
        log1.report('folders read: {}'.format(len(dataFolder.listFolders)))
        writeString2File(log1.fileNameMD,"## {}: {}\n".format(sessionName,param.param['acquisition']['label']),'a') # initialises MD file

        for currentFolder in dataFolder.listFolders:
            #currentFolder=dataFolder.listFolders[0]
            filesFolder=glob.glob(currentFolder+os.sep+'*.tif')
            dataFolder.createsFolders(currentFolder,param)
            dictShifts = {} #defaultdict(dict) # contains dictionary of shifts for each folder      
            
            # generates lists of files to process    
            param.files2Process(filesFolder)
            log1.report('About to process {} files\n'.format(len(param.fileList2Process)))
            writeString2File(dataFolder.outputFiles['alignImages']+'.bed', 
                          "File1 \t File_reference \t shift_y \t shift_x \t error \t diffphase",
                          'w')
        
            # Finds and loads Reference fiducial
            positionROIinformation=param.param['acquisition']['positionROIinformation']
            referenceBarcode=param.param['alignImages']['referenceFiducial']
            fileNameReferenceList, ROIList = RT2fileName(param.fileList2Process,referenceBarcode,positionROIinformation)    
            if len(fileNameReferenceList)>0:                        
                # loops over reference fiducials one ROI at a time
                for fileNameReference in fileNameReferenceList:
                    ROI = ROIList[fileNameReference]
                    imReference = Image()
                    imReference.loadImage2D(fileNameReference,log1,dataFolder)
                    dictShiftROI={}
                    
                    # loops over files in file list
                    for fileName in param.fileList2Process:
                        # excludes the reference fiducial and processes files in the same ROI
                        label=os.path.basename(fileName).split('_')[2]
                        if (fileName not in fileNameReference) and \
                                os.path.basename(fileName).split('_')[positionROIinformation]==ROI:
                            shift = align2Files(fileName,imReference,param,log1,session1,dataFolder,verbose)
                            session1.add(fileName,sessionName)
                            dictShiftROI[label]=shift.tolist()
                    dictShifts['ROI:'+ROI]=dictShiftROI
                    del imReference
                    
                #saves dicShifts
                saveJSON(dataFolder.outputFiles['dictShifts']+'.json',dictShifts)
            else:
                log1.report("Reference Barcode file does not exist: {}",format(referenceBarcode))

        del dataFolder, dictShifts

def appliesRegistrations(param,log1,session1):
    '''This function will 
    - load DAPI and barcode 2D projected images, 
    - apply registrations
    - save registered images as npy arrays 
    '''
    
    verbose=False
    if param.param['alignImages']['operation']=='overwrite':

        # session
        sessionName='registersImages'
    
        # processes folders and files 
        dataFolder=folders(param.param['rootFolder'])
        dataFolder.setsFolders()
        log1.addSimpleText("\n===================={}====================\n".format(sessionName))
        log1.report('folders read: {}'.format(len(dataFolder.listFolders)))
        #writeString2File(log1.fileNameMD,"## {}: {}\n".format(sessionName,param.param['acquisition']['label']),'a') # initialises MD file

        positionROIinformation=param.param['acquisition']['positionROIinformation']

        for currentFolder in dataFolder.listFolders:
            #currentFolder=dataFolder.listFolders[0] # only one folder processed so far...
            filesFolder=glob.glob(currentFolder+os.sep+'*.tif')
            dataFolder.createsFolders(currentFolder,param)

            #loads dicShifts with shifts for all ROIs and all labels
            dictShifts = loadJSON(dataFolder.outputFiles['dictShifts']+'.json')
        
            # generates lists of files to process    
            param.files2Process(filesFolder)
            log1.report('About to process {} files\n'.format(len(param.fileList2Process)))

            if len(param.fileList2Process)>0:                        
                # loops over files in file list
                 for fileName in param.fileList2Process:

                    # gets shift from dictionary
                    ROI = os.path.basename(fileName).split('_')[positionROIinformation]
                    label=os.path.basename(fileName).split('_')[2]
                    try:
                         shiftArray=dictShifts['ROI:'+ROI][label]
                    except KeyError:
                         shiftArray=None
                    
                    if shiftArray!=None:
                        
                        shift = np.asarray(shiftArray)
                        # loads 2D image and applies registration
                        Im = Image()
                        Im.loadImage2D(fileName,log1,dataFolder)
                        Im.data_2D = shiftImage(Im.data_2D,shift)
                        log1.report('Image registered using ROI:{}, label:{}, shift={}'.format(ROI,label,shift),'info')
                        
                        # saves registered 2D image
                        Im.saveImage2D(log1,dataFolder.outputFolders['alignImages'],tag='_2d_registered')

                        # logs output
                        session1.add(fileName,sessionName)
                    else:
                        log1.report('No shift found in dictionary for ROI:{}, label:{}'.format(ROI,label),'Warning')
                        
            del dataFolder





























