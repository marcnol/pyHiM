#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 15:18:30 2022
@author: Julian Gurgo
- This script takes trace tables and decodes using masks and/or combinatorial 
designs.
"""

# =============================================================================
# IMPORTS
# =============================================================================

from operator import length_hint
import string
from scipy.spatial import KDTree
import os
import numpy as np
import uuid
from astropy.table import Table
import ast
import glob
import re
from matrixOperations.chromatin_trace_table import ChromatinTraceTable
import copy


# =============================================================================
# CLASSES
# =============================================================================

class Decode:
    
    def __init__(self,dicPath=None,lenwords=None):
        
        self.lenwords=lenwords
        self.sourceDict = dicPath # full path: data_folder + os.sep + 'dictionary.dat'
        self.GroupSpotsID = None

    def DecodeCombinatorial(self,Tracemethod):

        if self.sourceDict == None:
            stringExit = 'Error: Provide a dictionary path to decode.'
            return stringExit

        if self.lenwords == None:
            stringExit = 'Error: Provide a valid word length.'
            return stringExit
        
        if Tracemethod == 'KDE':
            Traces = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_KDtree_ROI:+\d.ecsv', x)]### check name: Trace_combined_3D_method:KDTree_...
        elif Tracemethod == '3D_mask':
            Traces = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_mask:mask\d_ROI:+\d.ecsv', x)]### check name: Trace_combined_3D_method:mask_...
        elif Tracemethod == '3D_DAPI':
            Traces = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_mask:DAPI\d_ROI:+\d.ecsv', x)]### check name: Trace_combined_3D_method:DAPI_...

        # Get barcodes of grouped detections by the traces here (same trace ID). This is done per ROI
        # This should be later separated into class methods.
        for ROI in range(len(Traces)):
            
            TracesData = Table.read(Traces[ROI])
            # Get unique traces
            groups = list(set(TracesData['Trace_ID']))
            # For each unique trace, get barcode list. KDTree clustered barcodes share the same traceID
            self.BarcodeList = []
            self.GroupSpotsID = []
            for Trace in groups:
                self.BarcodeList.append(list(TracesData[TracesData['Trace_ID']==Trace]['Barcode #'])) # gets Barcode # for the same trace ID and puts them into a list 
                self.GroupSpotsID.append(list(TracesData[TracesData['Trace_ID']==Trace]['Spot_ID']))
             # Decode Elements in BarcodeList using the dictionary, and append the codeword to a new column or label column
            
            # All this has to be done by ROI
            self.wordsDetected = []

            # if self.sourceDict is not None:
            if self.sourceDict is not None:   
                
                # Assigns a word to each of the grouped detections
                for detection in self.BarcodeList:
                    word = [0]* self.lenwords
                    for item in detection:
                        word[item-1] = 1
                    self.wordsDetected.append(word)

                # Find words detected in dictionary and assign them a label
                print('\n >>> Assigns ID to detected words')
                file = open(self.sourceDict,'r')
                wordFile = file.readlines()
                words = [ast.literal_eval((line.split('\n')[0])) for line in wordFile]
                RTsID = ['RT' + str(idx+1) for idx in range(len(words))]
                dictionary = {k: v for k,v in zip(RTsID,words)}
                self.WordsID = [k for word in self.wordsDetected for k,v in dictionary.items() if v == word ]
                print('\n ID assignment for detected words done') 

                ##########################################################################################################################
                # Warning: Adding a col to matrix, maybe just modify label?
                ##########################################################################################################################
                # Adds a column with the RT ID of the word detected. 
                dataToModify = ['     ']*len(TracesData)
                TracesData.add_column(dataToModify,name='Decoded Word')
                for group in self.GroupSpotsID:
                    word = self.WordsID[self.GroupSpotsID.index(group)] ############################ MODIFY TO WRITE BINARY WORD
                    for spotID in group:
                        TracesData['Decoded Word'][np.where(TracesData['Spot_ID']==spotID)] = "{:<5}".format(word)

                TracesData.write(Traces[ROI].split('.')[0]+'_decoded.ecsv', format="ascii.ecsv")
                
        return self.wordsDetected, self.BarcodeList, self.GroupSpotsID
    
    ### First add to table the detected words and groups from previous step, then save it. And continue here
    
    def DecodeMasks(self,Tracemethod):
        
        # Case 1: We have a self.GroupSpotsID that comes from combinatorial labeling
        if self.GroupSpotsID is not None:

            print('Found combinatorial decoded words, decoding masks')
        
            maskList = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_mask', x)] #get masks files in dir
            masksNames = [x.split('mask:')[1].split('_')[0] for x in maskList]  # len(masksNames) is the number of zeros to add to the codeword. 

            '''
            The idea is to take the localizations that are saved in the KDTree file. 
            Then loop over the spots that are contained in a particular group of detections, and see if they are found in a particular mask.
            If they are, then make the mask ID 0 --> 1. 
            '''
            if Tracemethod == 'KDE':
                Traces = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_KDtree_ROI:+\d_decoded.ecsv', x)]### check name: Trace_combined_3D_method:KDTree_...
            elif Tracemethod == '3D_mask':
                Traces = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_mask:mask\d_ROI:+\d_decoded.ecsv', x)]### check name: Trace_combined_3D_method:mask_...
            elif Tracemethod == '3D_DAPI':
                Traces = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_mask:DAPI\d_ROI:+\d_decoded.ecsv', x)]### check name: Trace_combined_3D_method:DAPI_...

            for ROI in range(len(Traces)):

                TableKDTree = Table.read(Traces[ROI]) 
                TableWithMasks = TableKDTree.copy() # To add the detected mask ID

                MasksToAdd = [[None]*(self.lenwords+len(masksNames))]*len(TableWithMasks)
                TableWithMasks.add_column(MasksToAdd,name='Decoded Mask')

                self.wordsDetectedToExtend = copy.deepcopy(self.wordsDetected)
                #masks = []
                for group in self.GroupSpotsID: # select each group of detections, that has a combinatory word assigned
                    self.wordsDetectedToExtend[self.GroupSpotsID.index(group)].extend([0]*len(masksNames))
                    codeword = self.wordsDetectedToExtend[self.GroupSpotsID.index(group)].append([0]*len(masksNames)) # append as many 0's as masks detected to the codeword
                    for item in range(len(masksNames)): # iterate through masks to find if spots are there
                        #masksGroup = []# not necessary
                        inputMask = maskList[item]
                        maskTable = Table.read(inputMask,format='ascii') # read mask element
                        spotsInMask = list(maskTable['Spot_ID'])
                        for spot in range(len(group)): #check if given spots in group are located in a mask
                            if group[spot] in spotsInMask:
                                #TableWithMasks['Decoded Mask'][np.where(TableWithMasks['Spot_ID'] == group[spot])] = masksNames[item]
                                #masksGroup.append(masksNames[item])# not necessary
                                codeword[self.lenwords+masksNames.index(masksNames[item])] = 1 
                                TableWithMasks['Decoded Mask'][np.where(TableWithMasks['Spot_ID'] == group[spot])] = codeword
                        
                        #now the ID has been added to the codeword, and this can be generalized to many masks, where the 
                        # same spot can be located in more than one mask
                        # we do not need anymore to write the mask and RT in the table, rather this new codeword can be directly assigned to a locus 
                        # identifier, that has to be given to the program. 

                            else:
                                continue


        # Case 2: We do not have a self.GroupSpotsID from combinatorial labeling, proceed with mask decoding
        else:
            
            print('No combinatorial decoded words, decoding masks')

            maskList = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_mask', x)] #get masks files in dir
            masksNames = [x.split('mask:')[1].split('_')[0] for x in maskList]  # len(masksNames) is the number of zeros to add to the codeword. 

            #Traces = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_KDtree_ROI:+\d.ecsv', x)]

            if Tracemethod == 'KDE':
                Traces = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_KDtree_ROI:+\d_decoded.ecsv', x)]### check name: Trace_combined_3D_method:KDTree_...
            elif Tracemethod == '3D_mask':
                Traces = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_mask:mask\d_ROI:+\d_decoded.ecsv', x)]### check name: Trace_combined_3D_method:mask_...
            elif Tracemethod == '3D_DAPI':
                Traces = [x for x in os.listdir('.') if re.match('Trace_3D_barcode_mask:DAPI\d_ROI:+\d_decoded.ecsv', x)]### check name: Trace_combined_3D_method:DAPI_...

            for ROI in range(len(Traces)):
                
                TableKDTree = Table.read(Traces[ROI]) # or table obtained via any other method
                TableWithMasks = TableKDTree.copy() # To add the detected mask ID

                MasksToAdd = [[None]*len(masksNames)]*len(TableWithMasks)
                TableWithMasks.add_column(MasksToAdd,name='Decoded Mask')

                for item in range(len(masksNames)): # iterate through masks to find if spots are there
                    
                    codeword = list([0]*len(masksNames))
                    inputMask = maskList[item]
                    maskTable = Table.read(inputMask,format='ascii') # read mask element
                    spotsInMask = list(maskTable['Spot_ID'])

                    for spot in spotsInMask: #check if given spots in group are located in a mask

                        codeword[masksNames.index(masksNames[item])] = 1 
                        TableWithMasks['Decoded Mask'][np.where(TableWithMasks['Spot_ID'] == spot)] = codeword

        TableWithMasks.write(Traces[ROI].split('.')[0]+'_decoded_mask.ecsv', format="ascii.ecsv") #writes for all

                
        

        
        # def ()  to assign the correct locus to the word saved in 'Decoded Mask' column