#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 13:50:54 2021

@author: marcnol

unzips HiM_run.tar.gz recursively

In the command line, run as 
$ unzipHiM_run.py 

to unzip all the directories recursively

"""
import os
import glob
import argparse

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images, default: .")
    parser.add_argument("-R", "--recursive", help="One more depth of folders will be explored and zipped", action="store_true")
    
    args = parser.parse_args()

    if args.rootFolder:
        rootFolder0 = args.rootFolder
    else:
        rootFolder0 = os.getcwd()

    if args.recursive:
        recursive= args.recursive
    else:
        recursive= False

    if recursive:
        allFiles = glob.glob(rootFolder0+"/*")
        rootFolders = [x for x in allFiles if os.path.isdir(x)]
    else:
        rootFolders  = [rootFolder0]
        
    print("RootFolders: {}".format(rootFolders))        
    
    for rootFolder in rootFolders:
        
        # opens tarfile
        os.chdir(rootFolder)
        tarFileName = "HiMrun.tar.gz"
    
        if os.path.exists(rootFolder+os.sep+tarFileName):
            
            print("Unzipping archive: {}".format(tarFileName))
        
            
            # tar files in rootFolder
            filesMD = [os.path.basename(f) for f in glob.glob(rootFolder + os.sep + "HiM_analysis*.md", recursive=True)]
            # filesLOGMD = [os.path.basename(f) for f in glob.glob(rootFolder + os.sep + "log*.txt", recursive=True)]
            filesLOG = [os.path.basename(f) for f in glob.glob(rootFolder + os.sep + "HiM_analysis*.log", recursive=True)]
            filesSession = [os.path.basename(f) for f in glob.glob(rootFolder + os.sep + "Session*.json", recursive=True)]
        
            tarcmd = "tar -xzvf " + tarFileName 
            print("cmd> {}".format(tarcmd))       
       
            os.system(tarcmd)
        else:
            print("Nothing to unzip in: {}".format(rootFolder+os.sep+tarFileName))
