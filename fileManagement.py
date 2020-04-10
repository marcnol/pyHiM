#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:16:58 2020

@author: marcnol
"""

from datetime import datetime
import glob
import os
#import sys
from os import path
import json

class log:
    def __init__(self,fileName='test.log'):
        self.fileName=fileName
        
    def eraseFile(self):
        with open(self.fileName, 'w') as file:
            file.write("")
        
    # cmd line output only
    def info(self,text):
        print("INFO:{}".format(text))

   # saves to logfile, no display to cmd line
    def save(self,text='',status='info'):
        with open(self.fileName, 'a') as file:
            file.write(self.getFullString(text,status)+'\n')

    # thisfunction will output to cmd line and save in logfile
    def report(self,text,status='info'):
        print(self.getFullString(text,status))
        self.save("\n"+text,status)

    # returns formatted line to be outputed
    def getFullString(self,text='',status='info'):
        now=datetime.now()
        return "{}|{}>{}".format(now.strftime("%d/%m/%Y %H:%M:%S"),status,text)
    
    def addSimpleText(self,title):
        print("{}".format(title))
        with open(self.fileName, 'a') as file:
            file.write(title)
        
        
class folders():

    def __init__(self,masterFolder=r'/home/marcnol/Documents/Images'):
        self.masterFolder=masterFolder
        self.listFolders=[]; # list of sub-folders in rootFilder with images
        self.zProjectFolder=''
        self.outputFolders={}
        self.outputFiles={}
       
    # returns list of directories with given extensions
    def setsFolders(self,extension='tif'):
        
        # checks if there are files with the required extension in the root folder provided
        if os.path.isdir(self.masterFolder) and len(glob.glob(self.masterFolder+os.sep+'*.'+extension))>0:
            self.listFolders=self.masterFolder
        else:
            # finds more folders inside the given folder 
            hfolders = [folder for folder in glob.glob(self.masterFolder+os.sep+'*')
                       if os.path.isdir(folder) and 
                       len(glob.glob(folder+os.sep+'*.'+extension))>0] 
                       #os.path.name(folder)[0]!='F']
            self.listFolders=hfolders
        
    # creates folders for outputs
    def createsFolders(self,filesFolder,param):
        self.zProjectFolder=filesFolder+os.sep+param.param['zProject']['folder']
        self.createSingleFolder(self.zProjectFolder)            

        self.outputFolders['alignImages']=filesFolder+os.sep+param.param['alignImages']['folder']
        self.createSingleFolder(self.outputFolders['alignImages'])            

        self.outputFiles['alignImages']=self.outputFolders['alignImages']+os.sep+param.param['alignImages']['outputFile']
        self.outputFiles['dictShifts']=self.masterFolder+os.sep+param.param['alignImages']['outputFile']

    def createSingleFolder(self,folder):
        if not path.exists(folder):
            os.mkdir(folder)            
            print("Folder created: {}".format(folder))

class session:
    def __init__(self,name='dummy',fileName='session.json'):
        self.name=name
        self.fileName=fileName
        self.data={}   
        
    # loads existing session    
    def load(self):
        self.data = loadJSON(self.fileName)
        print("Session information read: {}".format(self.fileName))
    
    # saves session to file        
    def save(self,log):
        saveJSON(self.fileName,self.data)
        log.info("Saved json session file to {}".format(self.fileName))

   # add new task to session
    def add(self,key,value):
        if key not in self.data:
            self.data[key] = value
        else:
            self.data[key] = [self.data[key],value]
        
    def clearData(self):
        self.data={}
        
class Parameters:
    def __init__(self,label=''):
        self.label= label
        self.paramFile = "infoList.json"
        self.param = {
                        'image': {
                                    'currentPlane': 1,
                                    # 'contrastMin': 0.1,
                                    # 'contrastMax': 0.9,
                                    'claheGridH': 8,
                                    'claheGridW': 8
                        },
                        'acquisition': {
                                    'label': 'DAPI', # barcode, fiducial
                                    'positionROIinformation': 3
                        },
                        'zProject': {
                                    'folder':'zProject', # output folder
                                    'operation': 'skip', # overwrite, skip
                                    'mode': 'full', # full, manual, automatic
                                    'display': True,
                                    'saveImage': True,
                                    'zmin': 1,
                                    'zmax': 59,
                                    'zwindows': 10,
                                    'windowSecurity': 2,
                                    'zProjectOption': 'sum'  # sum or MIP
                        },
                       'alignImages': {
                                    'folder':'alignImages', # output folder
                                    'operation': 'overwrite', # overwrite, skip
                                    'outputFile': 'alignImages',
                                    'referenceFiducial': 'RT18'
                        }
                          
        }


    def get_param(self, param=False):
        if not param:
            return self.param
        else:
            return self.param[param]

    def initializeStandardParameters(self):
        with open(self.paramFile, 'w') as f:
            #json.dump(json.dumps(self.param), f, ensure_ascii=False, indent=4)
            json.dump(self.param, f, ensure_ascii=False, sort_keys=True, indent=4)
        print("Model parameters file saved to: {}".format(os.getcwd()+os.sep+self.paramFile))     
        
    def loadParametersFile(self, fileName):
        if path.exists(fileName):
            with open(fileName) as json_file:
                self.param = json.load(json_file)
                
            print("Parameters file read: {}".format(fileName))
    
    # method returns label specific filenames from filename list
    def files2Process(self, filesFolder):

        # selects DAPI files
        if self.param['acquisition']['label']=='DAPI':
            self.fileList2Process=[file for file in filesFolder 
                      if file.split('_')[-1].split('.')[0]=='ch00' and 'DAPI' in file.split('_')]

        # selects barcode files
        elif self.param['acquisition']['label']=='barcode':
            self.fileList2Process=[file for file in filesFolder 
               if len([i for i in file.split('_') if 'RT' in i])>0 and file.split('_')[-1].split('.')[0]=='ch01'] 

        # selects fiducial files
        elif self.param['acquisition']['label']=='fiducial':
            self.fileList2Process=[file for file in filesFolder 
               if (len([i for i in file.split('_') if 'RT' in i])>0 and file.split('_')[-1].split('.')[0]=='ch00') or 
                   ('DAPI' in file.split('_') and file.split('_')[-1].split('.')[0]=='ch01')] 


#################################################################################
# Functions
#################################################################################

        
def writeString2File(fileName,list2output,attribute='a'):
   with open(fileName,attribute) as fileHandle:
        fileHandle.write("{}\n".format(list2output))

def saveJSON(fileName,data):
    with open(fileName, 'w') as f:
        json.dump(data, f, ensure_ascii=False, sort_keys=True, indent=4)
   
def loadJSON(fileName):
    if path.exists(fileName):
        with open(fileName) as json_file:
            data = json.load(json_file)
    else:
        data = {}
    return data            
