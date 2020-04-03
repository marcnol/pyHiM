#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:16:58 2020

@author: marcnol
"""

from datetime import datetime

import glob,os,sys

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
        self.save(text,status)

    # returns formatted line to be outputed
    def getFullString(self,text='',status='info'):
        now=datetime.now()
        return "{}|{}>{}".format(now.strftime("%d/%m/%Y %H:%M:%S"),status,text)
    
    
    
class folders():

    def __init__(self,masterFolder=r'/home/marcnol/Documents/Images'):
        self.masterFolder=masterFolder
        self.listFolders=[];
        
    # returns list of directories with given extensions
    def setsFolders(self,extension='tif'):
        
        hfolders = [folder for folder in glob.glob(self.masterFolder+os.sep+'*')
                   if os.path.isdir(folder) and len(glob.glob(folder+os.sep+'*.'+extension))>0 and os.path.basename(folder)[0]!='F']
    
        self.listFolders=hfolders
