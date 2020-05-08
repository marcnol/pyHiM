#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:23:36 2020

@author: marcnol

test fitting barcode spots to masks


TO SOLVE:
    - I need to find a simple way of specifying the genomic coordinates of 
    barcodes for the production of the Hi-M matrix. At the moment it is just 
    using barcodeID from the file name.
    

"""

# =============================================================================
# IMPORTS
# =============================================================================


import glob,os
import argparse
import uuid
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from imageProcessing import Image
from fileManagement import folders,isnotebook
from fileManagement import session,writeString2File

from astropy.table import Table, vstack, Column
from astropy.visualization import SqrtStretch,simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import sigma_clip,sigma_clipped_stats

from photutils.segmentation import SegmentationImage
from sklearn.metrics import pairwise_distances

# =============================================================================
# CLASSES
# =============================================================================

class cellID():
    def __init__(self,barcodeMapROI,Masks,ROI):
        self.barcodeMapROI=barcodeMapROI
        self.Masks=Masks
        self.NcellsAssigned = 0
        self.NcellsUnAssigned = 0
        self.NbarcodesinMask=0

        self.SegmentationMask=SegmentationImage(self.Masks)
        self.numberMasks=self.SegmentationMask.nlabels
        self.ROI=ROI
        
        self.barcodesinMask=dict()
        for mask in range(self.numberMasks+1):
            self.barcodesinMask['maskID_'+str(mask)]=[]        

    def visualize(self):

        imageBarcodes=np.zeros([2048,2048])
        MasksBarcodes=Masks
        R=[]

        for i in range(len(self.barcodeMapROI.groups[0])):
            y_int=int(self.barcodeMapROI.groups[0]['xcentroid'][i])
            x_int=int(self.barcodeMapROI.groups[0]['ycentroid'][i])
            barcodeID = self.barcodeMapROI.groups[0]['Barcode #'][i]
            imageBarcodes[x_int][y_int] = barcodeID
            MasksBarcodes[x_int][y_int] += 20*barcodeID
            R.append([y_int,x_int,barcodeID])

        # Shows results
        Ra=np.array(R)
        plt.imshow(Masks, origin='lower', cmap='jet')
        plt.scatter(Ra[:,0],Ra[:,1],s=5,c=Ra[:,2],alpha=0.5)
      
    def alignByMasking(self):
        # [ Assigns barcodes to masks and creates <NbarcodesinMask> ]
        NbarcodesinMask=np.zeros(self.numberMasks+2)
        print('ROI:{}'.format(self.ROI))
        for i in range(len(self.barcodeMapROI.groups[0])):
            y_int=int(self.barcodeMapROI.groups[0]['xcentroid'][i])
            x_int=int(self.barcodeMapROI.groups[0]['ycentroid'][i])
            #barcodeID = self.barcodeMapROI.groups[ROI]['Barcode #'][i]
            maskID = self.Masks[x_int][y_int]
            self.barcodeMapROI['CellID #'][i]=maskID
            if maskID>0:
                NbarcodesinMask[maskID]+=1
                self.barcodesinMask['maskID_'+str(maskID)].append(i)
            
        # Total number of masks assigned and not assigned    
        self.NcellsAssigned = np.count_nonzero(NbarcodesinMask>0)
        self.NcellsUnAssigned = self.numberMasks - self.NcellsAssigned
        
        # this list contains which barcodes are allocated to which masks
        self.NbarcodesinMask=NbarcodesinMask
        
    def buildsdistanceMatrix(self,mode='mean'):
        '''
        

        
        '''
        print('building distance matrix')
        barcodeMapROI=self.barcodeMapROI

        # [ builds SCdistanceTable ]
        
        # sorts Table by cellID
        barcodeMapROI_cellID=barcodeMapROI.group_by('CellID #') # ROI data sorted by cellID
        ROIs,cellID,nBarcodes,barcodeIDs,p,cuid,buid=[],[],[],[],[],[],[]
        
        # iterates over all cell masks in an ROI  
        for key, group in zip(barcodeMapROI_cellID.groups.keys, barcodeMapROI_cellID.groups):
            
            if key['CellID #']>1: # excludes cellID 0 as this is background
                R=np.column_stack((np.array(group['xcentroid'].data),np.array(group['ycentroid'].data)))
                ROIs.append(group['ROI #'].data[0])
                cellID.append(key['CellID #'])
                nBarcodes.append(len(group))
                barcodeIDs.append(group['Barcode #'].data)
                buid.append(group['Buid'].data)
                p.append(pairwise_distances(R))
                cuid.append(str(uuid.uuid4())) # creates cell unique identifier
                #print("CellID #={}, nBarcodes={}".format(key['CellID #'],len(group)))

        SCdistanceTable=Table() #[],names=('CellID', 'barcode1', 'barcode2', 'distances'))
        SCdistanceTable['Cuid']=cuid
        SCdistanceTable['ROI #']=ROIs
        SCdistanceTable['CellID #']=cellID
        SCdistanceTable['nBarcodes']=nBarcodes
        SCdistanceTable['Barcode #']=barcodeIDs
        SCdistanceTable['Buid']=buid
        SCdistanceTable['PWDmatrix']=p
        
        # NEEED TO ADD THE IDENTITIES OF THE BARCODES
        
        self.SCdistanceTable=SCdistanceTable

        print('Cells with barcodes found: {}'.format(len(SCdistanceTable)))
        
        # [ builds SCmatrix ]
        numberMatrices=len(SCdistanceTable) # z dimensions of SCmatrix 
        uniqueBarcodes=np.unique(barcodeMapROI['Barcode #'].data)
        numberUniqueBarcodes=uniqueBarcodes.shape[0] # number of unique Barcodes for xy dimensions of SCmatrix
        SCmatrix=np.zeros((numberUniqueBarcodes,numberUniqueBarcodes,numberMatrices))
        SCmatrix[:]=np.NaN
        
        for iCell,scPWDitem in zip(range(numberMatrices),SCdistanceTable):
            barcodes2Process=scPWDitem['Barcode #']
            for barcode1,ibarcode1 in zip(barcodes2Process,range(len(barcodes2Process))):
                indexBarcode1=np.nonzero(uniqueBarcodes==barcode1)[0][0]
                for barcode2,ibarcode2 in zip(barcodes2Process,range(len(barcodes2Process))):
                    indexBarcode2=np.nonzero(uniqueBarcodes==barcode2)[0][0]
                    if barcode1!=barcode2:
                        newdistance=scPWDitem['PWDmatrix'][ibarcode1][ibarcode2]
                        if mode=='last':
                            SCmatrix[indexBarcode1][indexBarcode2][iCell]=newdistance
                        elif mode=='mean':
                            SCmatrix[indexBarcode1][indexBarcode2][iCell]=np.nanmean([newdistance,SCmatrix[indexBarcode1][indexBarcode2][iCell]])
                        elif mode=='min':
                            SCmatrix[indexBarcode1][indexBarcode2][iCell]=np.nanmin([newdistance,SCmatrix[indexBarcode1][indexBarcode2][iCell]])
       
        self.SCmatrix=SCmatrix
        self.meanSCmatrix=np.nanmean(SCmatrix,axis=2)
        self.uniqueBarcodes=uniqueBarcodes
        
# =============================================================================
# FUNCTIONS
# =============================================================================

def plotMatrix(SCmatrixCollated,uniqueBarcodes, pixelSize, numberROIs=1, outputFileName='test',logNameMD='log.md',clim=1.4,cm='seismic',figtitle='PWD matrix',cmtitle='distance, um',nCells=0):

    # projects matrix by calculating median in the nCell direction
    if len(SCmatrixCollated.shape)==3:
        meanSCmatrix = pixelSize*np.nanmedian(SCmatrixCollated,axis=2)
        nCells=SCmatrixCollated.shape[2]
    else:
        meanSCmatrix = pixelSize*SCmatrixCollated

    # plots figure   
    fig=plt.figure(figsize=(10, 10))
    pos = plt.imshow(meanSCmatrix,cmap=cm) # colormaps RdBu seismic
    plt.xlabel('barcode #')
    plt.ylabel('barcode #')
    plt.title(figtitle+ ' | ' + str(meanSCmatrix.shape[0]) + ' barcodes | n='+str(nCells)+' | ROIs='+str(numberROIs))
    plt.xticks(np.arange(SCmatrixCollated.shape[0]), uniqueBarcodes)
    plt.yticks(np.arange(SCmatrixCollated.shape[0]), uniqueBarcodes)
    cbar=plt.colorbar(pos,fraction=0.046, pad=0.04)
    cbar.minorticks_on()
    cbar.set_label(cmtitle)
    plt.clim(0,clim)

    plt.savefig(outputFileName+'_HiMmatrix.png')
    
    if not isnotebook():
        plt.close()
    
    writeString2File(logNameMD,"![]({})\n".format(outputFileName+'_HiMmatrix.png'),'a')

def plotDistanceHistograms(SCmatrixCollated,pixelSize,outputFileName='test',logNameMD='log.md'):
    
    if not isnotebook():
        NplotsX = NplotsY =  SCmatrixCollated.shape[0]
    else:
        NplotsX = NplotsY =  min([10,SCmatrixCollated.shape[0]]) # sets a max of subplots if you are outputing to screen!
            
    bins=np.arange(0,4,0.25)
    
    sizeX,sizeY=NplotsX*4, NplotsY*4
    
    fig, axs = plt.subplots(figsize=(sizeX,sizeY), ncols=NplotsX, nrows=NplotsY, sharex=True)
    
    for i in range(NplotsX):
        for j in range(NplotsY):
            if i!=j:
                #print('Printing [{}:{}]'.format(i,j))
                axs[i,j].hist(pixelSize*SCmatrixCollated[i,j,:],bins=bins)
                    
    plt.xlabel('distance, um')
    plt.ylabel('counts')
    plt.savefig(outputFileName+'_PWDhistograms.png')
   
    if not isnotebook():
       plt.close()

    writeString2File(logNameMD,"![]({})\n".format(outputFileName+'_PWDhistograms.png'),'a')

def buildsPWDmatrix(currentFolder, fileNameBarcodeCoordinates, outputFileName, pixelSize=0.1, logNameMD='log.md'):     
    
    # Processes Tables 
    barcodeMap = Table.read(fileNameBarcodeCoordinates,format='ascii.ecsv')
    barcodeMapROI=barcodeMap.group_by('ROI #')
    
    SCmatrixCollated, uniqueBarcodes = [], []
    numberROIs=len(barcodeMapROI.groups.keys)
    filesinFolder=glob.glob(currentFolder+os.sep+'*.tif')

    for ROI in range(numberROIs):
        nROI=barcodeMapROI.groups.keys[ROI][0] # need to iterate over the first index
    
        print('\nROIs detected: {}'.format(barcodeMapROI.groups.keys))
        
        barcodeMapSingleROI=barcodeMap.group_by('ROI #').groups[ROI]
        
        # finds file for masks
        fileList2Process=[file for file in filesinFolder 
                  if file.split('_')[-1].split('.')[0]=='ch00' and
                  'DAPI' in file.split('_')
                  and int(os.path.basename(file).split('_')[3])==nROI]
        
        if len(fileList2Process)>0:
        
            # loads Masks
            fileNameROImasks = os.path.basename(fileList2Process[0]).split('.')[0]+'_Masks.npy'
            fullFileNameROImasks = os.path.dirname(fileNameBarcodeCoordinates)+os.sep+fileNameROImasks 
            if os.path.exists(fullFileNameROImasks):
                Masks=np.load(fullFileNameROImasks)
                
                # Assigns barcodes to Masks for a given ROI
                cellROI = cellID(barcodeMapSingleROI,Masks,ROI)
                
                cellROI.alignByMasking()
                
                cellROI.buildsdistanceMatrix('min') # mean min last
                
                print('ROI: {}, N cells assigned: {} out of {}'.format(ROI,cellROI.NcellsAssigned,cellROI.numberMasks))
            
                uniqueBarcodes = cellROI.uniqueBarcodes    
                
                # saves Table with results per ROI
                
                cellROI.SCdistanceTable.write(outputFileName+'_ROI'+str(nROI)+'.ecsv',format='ascii.ecsv',overwrite=True)
    
                if len(SCmatrixCollated)>0:
                    SCmatrixCollated=np.concatenate((SCmatrixCollated,cellROI.SCmatrix),axis=2)
                else:
                    SCmatrixCollated=cellROI.SCmatrix
                del cellROI
            else:
                print('Error, no DAPI mask file found for {}\n Expected: {}',format(fileNameBarcodeCoordinates, fullFileNameROImasks))
                
    # saves output
    np.save(outputFileName+'_HiMscMatrix.npy',SCmatrixCollated)
    np.savetxt(outputFileName+'_uniqueBarcodes.ecsv', uniqueBarcodes, delimiter=" ",fmt='%d')
   
    # plots outputs
    plotMatrix(SCmatrixCollated,uniqueBarcodes, pixelSize, numberROIs, outputFileName,logNameMD)
    plotDistanceHistograms(SCmatrixCollated, pixelSize, outputFileName,logNameMD)
    

def processesPWDmatrices(param,log1,session1):
    sessionName='buildsPWDmatrix'
 
    # processes folders and files 
    dataFolder=folders(param.param['rootFolder'])
    log1.addSimpleText("\n===================={}====================\n".format(sessionName))
    log1.report('folders read: {}'.format(len(dataFolder.listFolders)))
    writeString2File(log1.fileNameMD,"## {}\n".format(sessionName),'a') 
 
    for currentFolder in dataFolder.listFolders:
        #filesFolder=glob.glob(currentFolder+os.sep+'*.tif')
        dataFolder.createsFolders(currentFolder,param)
        log1.report("-------> Processing Folder: {}".format(currentFolder))

        fileNameBarcodeCoordinates = dataFolder.outputFiles['segmentedObjects']+'_barcode.dat'
        #segmentedMasksFolder=dataFolder.outputFolders['segmentedObjects']
        outputFileName = dataFolder.outputFiles['buildsPWDmatrix']
        pixelSize=0.1 

        buildsPWDmatrix(currentFolder, fileNameBarcodeCoordinates,outputFileName,pixelSize,log1.fileNameMD)
        session1.add(currentFolder,sessionName)
     
        log1.report('HiM matrix in {} processed'.format(currentFolder),'info')

# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
        
    parser = argparse.ArgumentParser()
    parser.add_argument("-F","--rootFolder", help="Folder with images")
    args = parser.parse_args()
    
    print("\n--------------------------------------------------------------------------")
    
    if args.rootFolder:
        rootFolder=args.rootFolder
    else:
        rootFolder='/home/marcnol/data/Experiment_20/Embryo_1'
        #rootFolder='/home/marcnol/data/Experiment_15/Embryo_006_ROI18'
        #rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset'
    print("parameters> rootFolder: {}".format(rootFolder))

    # sets name for coordinate file and for masks
    segmentedMasksFolder = rootFolder +'/rawData/segmentedObjects/'
    fileNameBarcodeCoordinates = segmentedMasksFolder + 'segmentedObjects_barcode.dat'
    #currentFolder=glob.glob(rootFolder+'/rawData/'+'*.tif')
    currentFolder=rootFolder+'/rawData'

    outputFileName = rootFolder.split('/')[-1]
    pixelSize=0.1
    
    buildsPWDmatrix(currentFolder, fileNameBarcodeCoordinates,outputFileName,pixelSize)
           