 

###################################################################################
#IMPORTS
###################################################################################
import numpy as np
import pandas as pd
import matplotlib as plt
import random
import ast
from astropy.table import Table
from astropy.table import vstack
import uuid
import os

###################################################################################
#CLASSES
###################################################################################
class ChromatinTraceTable:

    def __init__(self, XYZ_unit='micron', genome_assembly='dm3'):
        self.a = 1
        self.XYZ_unit=XYZ_unit
        self.genome_assembly=genome_assembly

    def initialize(self):
        self.data = Table(
                names=(
                    "Spot_ID",
                    "Trace_ID",
                    "x",
                    "y",
                    "z",
                    "Chrom",
                    "Chrom_Start",
                    "Chrom_End",
                    "ROI #",
                    "Mask_id",
                    "Barcode #",
                    "label",
                ),
                dtype=("S2",
                       "S2",
                       "f4",
                       "f4",
                       "f4",
                       "S2",
                       "int",
                       "int",
                       "int",
                       "int",
                       "int",
                       "S2",
                       ),
            )

        self.data.meta['comments']=["XYZ_unit={}".format(self.XYZ_unit),
                                    "genome_assembly={}".format(self.genome_assembly),
                                    ]
    def add_row(self,entry):
        self.data.add_row(entry)
    
    def append(self, table):
        """
        appends <table> to self.data
        Parameters
        ----------
        table : astropy table
            table to append to existing self.data table.
        Returns
        -------
        None.
        """

        self.data = vstack([self.data, table])
        
class CreateTraceTable:
    
    def __init__(self, dictionary, wordsToUse=6, sizeFOV=2000):
        self.wordsToUse = wordsToUse
        self.sizeFOV = sizeFOV
        self.dictionary = dictionary
        
    def GetWordsFromDic(self):
        
        self.experimentWords = [self.dictionary[np.random.randint(0,len(self.dictionary))] for idx in range(self.wordsToUse) ]
        
        return self.experimentWords
        
    def GenerateSpotsPosition(self):
        self.coordinates = []
        for spot in range(self.wordsToUse):
            detection = list(np.random.uniform(0,self.sizeFOV,3))
            self.coordinates.append(detection)
        # return self.coordinates
            
    def GenerateTable(self):
        
        self.outputTable = ChromatinTraceTable()
        self.outputTable.initialize()
        counts = 0
        for word in self.experimentWords:
            traceID = str(uuid.uuid4())
            maskID = self.experimentWords.index(word)
            for HybRound in range(len(word)):
                if word[HybRound] == 1:
                    
                    roi = 1
                    barcode = HybRound + 1
                    # z,x,y = spots_subpixel[i,:]
                    Table_entry = [str(uuid.uuid4()),
                                    traceID,
                                    self.coordinates[self.experimentWords.index(word)][0],#x
                                    self.coordinates[self.experimentWords.index(word)][1],#y
                                    self.coordinates[self.experimentWords.index(word)][2],#z
                                    'Chr3R',
                                    129999, # chr start
                                    149999, # chr end
                                    roi, # ROI #
                                    maskID, #mask id
                                    barcode, #barcde #
                                    '    ' #label
                                    ]
                    counts = counts + 1
                    self.outputTable.add_row(Table_entry)
                    # outputTable.append(Table_entry)
        return self.outputTable
    
    def SaveTable(self,outputPath,outputName):
        
        outputFileName = outputPath + os.sep + outputName
        self.outputTable.data.write(outputFileName,format="ascii.ecsv")
        
    def SaveChosenWords(self,wordsPath):
        
        f = open(wordsPath + os.sep + 'codewords.dat','w')
        for element in self.experimentWords:
            f.write(str(element) + '\n')
        f.close()
        
###################################################################################
#FUNCTIONS
###################################################################################

# # Construct and dictionary with RT names ##############
# lenWords = 16                                         #
# keys = ['RT' + str(idx+1) for idx in range(lenWords)] # Not really useful, appart from keys
# words = {k:[] for k in keys}                          #
# #######################################################

def ProcessDictionary(dicpath):
    
    new_df = pd.read_excel(dicpath,header=1)
    dictionaryString = new_df['Codewords'].tolist()
    dictstring=[z0.split('  ') for z0 in dictionaryString]
    dictionary = [list(map(int, ze)) for ze in dictstring]
    
    return dictionary

#%%
# in a first approximation, give all spots the same coordinate
dicpath = '/mnt/PALM_dataserv/DATA/gurgo/Projects_2022/Multiplexed_image_analysis/dictionary_chen_moffit/dictionary_moffit.xlsx'
dictionary = ProcessDictionary(dicpath)

TraceTable = CreateTraceTable(dictionary)
TraceTable.GetWordsFromDic()
TraceTable.GenerateSpotsPosition()
# outputTable = TraceTable.GenerateTable()
TraceTable.GenerateTable()

outputPath = '/mnt/PALM_dataserv/DATA/gurgo/Projects_2022/Multiplexed_image_analysis/simulated_dataset'
outputName = "Trace_3D_barcode_KDtree_ROI:2.ecsv"
TraceTable.SaveTable(outputPath,outputName)
TraceTable.SaveChosenWords(outputPath)
#%% Create fake mask files, assigning words from previous step to 2 different masks

# Read KDTree data

TableKDTree = Table.read('Trace_3D_barcode_KDtree_ROI:2_decoded.ecsv')

name1 = 'Trace_3D_barcode_mask:mask1_ROI:2.ecsv'
name2 = 'Trace_3D_barcode_mask:mask2_ROI:2.ecsv'

detectedRTs = set(TableKDTree['Decoded Word'])
mask1 =  set(random.sample(detectedRTs, 3))
mask2 = detectedRTs - mask1

mask1Members = list(mask1)
mask2Members = list(mask2)

for element in mask1Members:
    TableKDTree['Mask_id'][np.where(TableKDTree['Decoded Word']==element)] = 1

for element in mask2Members:
    TableKDTree['Mask_id'][np.where(TableKDTree['Decoded Word']==element)] = 2

# print(mask1Members)
# print(mask2Members)
    
# Drop elements of mask 2 to save a file for mask 1, and viceversa
TableMask1 = TableKDTree.copy()
TableMask1.remove_rows(tuple(np.where(TableMask1['Mask_id'] == 2)))
TableMask1.write(name1, format="ascii.ecsv")

TableMask2 = TableKDTree.copy()
TableMask2.remove_rows(tuple(np.where(TableMask2['Mask_id'] == 1)))
TableMask2.write(name2, format="ascii.ecsv")