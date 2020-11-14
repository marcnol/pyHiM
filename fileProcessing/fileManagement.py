#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:16:58 2020

@author: marcnol

Classes and functions for file management

"""
# =============================================================================
# IMPORTS
# =============================================================================

from datetime import datetime
import glob
import os
from os import path
import json
import re
from warnings import warn
import multiprocessing
import numpy as np

from dask.distributed import Client, LocalCluster


# =============================================================================
# CLASSES
# =============================================================================


class log:
    def __init__(self, rootFolder="./", fileNameRoot="HiM_analysis", parallel=False):
        now = datetime.now()
        dateTime = now.strftime("%Y%m%d_%H%M%S")
        self.fileName = rootFolder + os.sep + fileNameRoot + dateTime + ".log"
        self.fileNameMD = self.fileName.split(".")[0] + ".md"
        self.parallel=parallel
        self.eraseFile()
        self.report("Starting to log to: {}".format(self.fileName))
        
    def eraseFile(self):
        # with open(self.fileName, 'w') as file:
        #    file.write("")
        writeString2File(self.fileName, "", "w")

    # cmd line output only
    def info(self, text):
        print("INFO:{}".format(text))

    # saves to logfile, no display to cmd line
    def save(self, text="", status="info"):
        # with open(self.fileName, 'a') as file:
        #    file.write(self.getFullString(text,status)+'\n')
        writeString2File(self.fileName, self.getFullString(text, status), "a")

    # thisfunction will output to cmd line and save in logfile
    def report(self, text, status="info"):
        if not self.parallel or status.lower()=="error":
            print(self.getFullString(text, status))
            self.save("\n" + text, status)
        else:
            self.save("\n" + text, status)
       
    # returns formatted line to be outputed
    def getFullString(self, text="", status="info"):
        now = datetime.now()
        return "{}> {}".format(now.strftime("%d/%m/%Y %H:%M:%S"), text)

    def addSimpleText(self, title):
        print("{}".format(title))
        writeString2File(self.fileName, title, "a")


class folders:
    def __init__(self, masterFolder=r"/home/marcnol/Documents/Images"):
        self.masterFolder = masterFolder
        self.listFolders = []
        
        # list of sub-folders in rootFilder with images
        self.zProjectFolder = ""
        self.outputFolders = {}
        self.outputFiles = {}
        self.setsFolders()

    # returns list of directories with given extensions
    def setsFolders(self, extension="tif"):

        # finds more folders inside the given folder
        hfolders = [
            folder
            for folder in glob.glob(self.masterFolder + os.sep + "*")
            if os.path.isdir(folder) and len(glob.glob(folder + os.sep + "*." + extension)) > 0
        ]
        # os.path.name(folder)[0]!='F']
        if len(hfolders) > 0:
            self.listFolders = hfolders
        else:
            self.listFolders = []

        # checks if there are files with the required extension in the root folder provided
        if os.path.isdir(self.masterFolder) and len(glob.glob(self.masterFolder + os.sep + "*." + extension)) > 0:
            # self.listFolders=self.masterFolder
            self.listFolders.append(self.masterFolder)

        print("\n> SetsFolders> Detected {} folders with images".format(len(self.listFolders)))

    # creates folders for outputs
    def createsFolders(self, filesFolder, param):
        """
        this function will create all the folders required for processingPipeline

        Parameters
        ----------
        filesFolder : string 
            rootFolder
        param : Parameters Class 
            with filenames of folders to be created

        Returns
        -------
        None.

        """
        self.outputFolders["zProject"] = filesFolder + os.sep + param.param["zProject"]["folder"]
        self.outputFolders["alignImages"] = filesFolder + os.sep + param.param["alignImages"]["folder"]
        self.outputFolders["segmentedObjects"] = filesFolder + os.sep + param.param["segmentedObjects"]["folder"]
        self.outputFolders["buildsPWDmatrix"] = filesFolder + os.sep + "buildsPWDmatrix"
        self.outputFolders["projectsBarcodes"] = filesFolder + os.sep + param.param["projectsBarcodes"]["folder"]

        self.createSingleFolder(self.outputFolders["zProject"])
        self.createSingleFolder(self.outputFolders["alignImages"])
        self.createSingleFolder(self.outputFolders["segmentedObjects"])
        self.createSingleFolder(self.outputFolders["buildsPWDmatrix"])
        self.createSingleFolder(self.outputFolders["projectsBarcodes"])

        # self.outputFiles['zProject']=self.outputFolders['zProject']+os.sep+param.param['zProject']['outputFile']
        self.outputFiles["alignImages"] = (
            self.outputFolders["alignImages"] + os.sep + param.param["alignImages"]["outputFile"]
        )
        self.outputFiles["dictShifts"] = self.outputFolders["alignImages"] + os.sep + param.param["alignImages"]["outputFile"]
        self.outputFiles["segmentedObjects"] = (
            self.outputFolders["segmentedObjects"] + os.sep + param.param["segmentedObjects"]["outputFile"]
        )
        self.outputFiles["buildsPWDmatrix"] = self.outputFolders["buildsPWDmatrix"] + os.sep + "buildsPWDmatrix"
        self.outputFiles["projectsBarcodes"] = (
            self.outputFolders["projectsBarcodes"] + os.sep + param.param["projectsBarcodes"]["outputFile"]
        )

    def createSingleFolder(self, folder):
        if not path.exists(folder):
            os.mkdir(folder)
            print("Folder created: {}".format(folder))


class session:
    def __init__(self, rootFolder, name="dummy"):
        now = datetime.now()
        sessionRootName = now.strftime("%d%m%Y_%H%M%S")
        self.fileName = rootFolder + os.sep + "Session_" + sessionRootName + ".json"
        self.name = name
        self.data = {}

    # loads existing session
    def load(self):
        self.data = loadJSON(self.fileName)
        print("Session information read: {}".format(self.fileName))

    # saves session to file
    def save(self, log):
        saveJSON(self.fileName, self.data)
        log.info("Saved json session file to {}".format(self.fileName))

    # add new task to session
    def add(self, key, value):
        if key not in self.data:
            self.data[key] = value
        else:
            self.data[key] = [self.data[key], value]

    def clearData(self):
        self.data = {}


class FileHandling:
    def __init__(self, fileName):
        self.fileName = fileName
        self.positionROIinformation = 3

    def getROI(self):
        return os.path.basename(self.fileName).split("_")[self.positionROIinformation]

class Parameters:
    def __init__(self, rootFolder="./", label=""):
        self.label = label
        self.paramFile = "infoList.json"
        self.param = {
            "acquisition": {
                "label": "DAPI",
                "positionROIinformation": 3,
                "fileNameRegExp": "scan_(?P<runNumber>[0-9]+)_(?P<cycle>[\\w|-]+)_(?P<roi>[0-9]+)_ROI_converted_decon_(?P<channel>[\\w|-]+).tif",
                "DAPI_channel": "ch00",
                "fiducialDAPI_channel": "ch01",
                "RNA_channel": "ch02",
                "fiducialBarcode_channel": "ch00",
                "barcode_channel": "ch01",
                "pixelSizeXY": 0.1,
                "pixelSizeZ": 0.25,
                },  # barcode, fiducial
            "zProject": {
                "folder": "zProject",  # output folder
                "operation": "skip",  # overwrite, skip
                "mode": "full",  # full, manual, automatic
                "display": True,
                "saveImage": True,
                "zmin": 1,
                "zmax": 59,
                "zwindows": 10,
                "windowSecurity": 2,
                "zProjectOption": "sum",  # sum or MIP
            },
            "alignImages": { 
                "folder": "alignImages",  # output folder
                "operation": "overwrite",  # overwrite, skip
                "outputFile": "alignImages",
                "referenceFiducial": "RT18",
                "alignByBlock": True, # alignByBlock True will perform block alignment
                "tolerance": 0.1, #Used in blockAlignment to determine the % of error tolerated
                "lower_threshold": 0.999, # lower threshold to adjust image intensity levels before xcorrelation
                "higher_threshold": 0.9999999, # higher threshold to adjust image intensity levels before xcorrelation
                "localShiftTolerance": 1,
                "bezel": 20,                
            },
            "projectsBarcodes": {
                "folder": "projectsBarcodes",  # output folder
                "operation": "overwrite",  # overwrite, skip
                "outputFile": "projectsBarcodes",
            },
            "segmentedObjects": {
                "folder": "segmentedObjects",  # output folder
                "operation": "overwrite",  # overwrite, skip
                "outputFile": "segmentedObjects",
                "background_method": "inhomogeneous",  # flat or inhomogeneous or stardist
                "stardist_network": "stardist_nc14_nrays:64_epochs:40_grid:2",
                "stardist_basename": "/mnt/grey/DATA/users/marcnol/models",
                "tesselation": True,  # tesselates DAPI masks
                "background_sigma": 3.0,  # used to remove inhom background
                "threshold_over_std": 1.0,  # threshold used to detect sources
                "fwhm": 3.0,  # source size in px
                "brightest": 1100,  # max number of objects segmented per FOV
                "intensity_min": 0,  # min int to keep object
                "intensity_max": 59,  # max int to keeep object
                "area_min": 50,  # min area to keeep object
                "area_max": 500,  # max area to keeep object
                "flux_min": 200,  # min flux to keeep object                
                "residual_max": 2.5,  # max residuals to keeep object                
                "sigma_max": 5,  # max sigma 3D fitting to keeep object                
                "centroidDifference_max": 5,  # max diff between Moment and Gaussian z fits to keeep object                
                "3DGaussianfitWindow": 3,  # size of window to extract subVolume, px. 3 means subvolume will be 7x7.
                "toleranceDrift":1, # tolerance used for block drift correction, in px
            },
        }
        self.initializeStandardParameters()
        self.paramFile = rootFolder + os.sep + label
        self.loadParametersFile(self.paramFile)
        self.param["rootFolder"] = rootFolder
        self.fileParts={}
        
    def get_param(self, param=False):
        if not param:
            return self.param
        else:
            return self.param[param]

    def initializeStandardParameters(self):
        with open(self.paramFile, "w") as f:
            # json.dump(json.dumps(self.param), f, ensure_ascii=False, indent=4)
            json.dump(self.param, f, ensure_ascii=False, sort_keys=True, indent=4)
        print("Model parameters file saved to: {}".format(os.getcwd() + os.sep + self.paramFile))

    def loadParametersFile(self, fileName):
        if path.exists(fileName):
            with open(fileName) as json_file:
                self.param = json.load(json_file)

            print("Parameters file read: {}".format(fileName))

    def setsChannel(self, key, default):
        if key in self.param["acquisition"].keys():
            channel = self.param["acquisition"][key]
        else:
            channel = default

        return channel

    # method returns label specific filenames from filename list
    def files2Process(self, filesFolder):

        # defines channel for DAPI, fiducials and barcodes
        channelDAPI = self.setsChannel("DAPI_channel", "ch00")
        channelbarcode = self.setsChannel("barcode_channel", "ch01")
        channelfiducial = self.setsChannel("fiducialBarcode_channel", "ch00")

        # finds if there is 2 or 3 channels for DAPI acquisition
        fileList2Process = [
            file for file in filesFolder if self.decodesFileParts(path.basename(file))["channel"] == "ch02" and "DAPI" in file.split("_")
        ]

        # defines channels for RNA and DAPI-fiducial
        if len(fileList2Process) > 0:
            channelDAPI_fiducial = self.setsChannel("fiducialDAPI_channel", "ch02")
            channelDAPI_RNA = self.setsChannel("fiducialDAPI_channel", "ch01")
        else:
            channelDAPI_fiducial = self.setsChannel("fiducialDAPI_channel", "ch01")
            channelDAPI_RNA = self.setsChannel("fiducialDAPI_channel", "ch04")

        if channelDAPI_fiducial and len(fileList2Process)==0:
            warn("\n\n****You are using ch02 for channelDAPI_fiducial but there are only 2 channels for DAPI!\n\n")
        
        # selects DAPI files
        if self.param["acquisition"]["label"] == "DAPI":
            self.fileList2Process = [
                file
                for file in filesFolder
                if self.decodesFileParts(path.basename(file))["channel"] == channelDAPI and "DAPI" in file.split("_")
            ]

        # selects DAPIch2 files
        if self.param["acquisition"]["label"] == "RNA":
            self.fileList2Process = [
                file
                for file in filesFolder
                if self.decodesFileParts(path.basename(file))["channel"]== channelDAPI_RNA and "DAPI" in file.split("_")
            ]

        # selects barcode files
        elif self.param["acquisition"]["label"] == "barcode":
            self.fileList2Process = [
                file
                for file in filesFolder
                if len([i for i in file.split("_") if "RT" in i]) > 0
                and self.decodesFileParts(path.basename(file))["channel"] == channelbarcode
            ]

        # selects fiducial files
        elif self.param["acquisition"]["label"] == "fiducial":
            self.fileList2Process = [
                file
                for file in filesFolder
                if (
                    len([i for i in file.split("_") if "RT" in i]) > 0
                    and self.decodesFileParts(path.basename(file))["channel"]  == channelfiducial
                )
                or ("DAPI" in file.split("_") and self.decodesFileParts(path.basename(file))["channel"] == channelDAPI_fiducial)
            ]


    def decodesFileParts(self, fileName):
        '''
        decodes variables from an input file. typically, RE takes the form:
            
        "DAPI_(?P<runNumber>[0-9]+)_(?P<cycle>[\w|-]+)_(?P<roi>[0-9]+)_ROI_converted_decon_(?P<channel>[\w|-]+).tif"
        
        thus, by running decodesFileParts(param,fileName) you will get back either an empty dict if the RE were not present 
        in your infoList...json file or a dict as follows if it all worked out fine:
            
        fileParts['runNumber']: runNumber number
        fileParts['cycle']: cycle string
        fileParts['roi']: roi number
        fileParts['channel']: channel string
        
        Parameters
        ----------
        param : Parameters class
        fileName : string
            filename to decode
    
        Returns
        -------
        Dict with fileParts.
    
        '''
        # decodes regular expressions
        if 'fileNameRegExp' in self.param['acquisition'].keys():
            fileParts=re.search(self.param['acquisition']['fileNameRegExp'],fileName)
            return fileParts
        else:
            return {}
    
class daskCluster:
    def __init__(self, requestedNumberNodes,maximumLoad=0.6,memoryPerWorker = 2000):
        self.requestedNumberNodes = requestedNumberNodes
        # self.nThreads will be created after exetution of initializeCluster()
        self.maximumLoad = maximumLoad  # max number of workers that I can take
        self.memoryPerWorker = memoryPerWorker# in Mb       
        self.initializeCluster()

    def initializeCluster(self):
        
        numberCoresAvailable = multiprocessing.cpu_count()

        # we want at least 2 GB per worker
        _, _, free_m = map(int, os.popen("free -t -m").readlines()[-1].split()[1:])
 
        maxNumberThreads = int(np.min([numberCoresAvailable*self.maximumLoad, free_m / self.memoryPerWorker]))
        
        self.nThreads = int(np.min([maxNumberThreads, self.requestedNumberNodes]))

        print("Cluster with {} workers started ({} requested)".format(self.nThreads,self.requestedNumberNodes))

    def createDistributedClient(self):
        self.cluster = LocalCluster(n_workers=self.nThreads,
                                # processes=True,
                                # threads_per_worker=1,
                                # memory_limit='2GB',
                                # ip='tcp://localhost:8787',
                                ) 
        self.client = Client(self.cluster)
        

# =============================================================================
# FUNCTIONS
# =============================================================================
def writeString2File(fileName, list2output, attribute="a"):
    with open(fileName, attribute) as fileHandle:
        fileHandle.write("{}\n".format(list2output))


def saveJSON(fileName, data):
    with open(fileName, "w") as f:
        json.dump(data, f, ensure_ascii=False, sort_keys=True, indent=4)


def loadJSON(fileName):
    if path.exists(fileName):
        with open(fileName) as json_file:
            data = json.load(json_file)
    else:
        data = {}
    return data


def isnotebook():
    """
    This function detects if you are running on an ipython console or in the shell.
    It is used to either kill plots or leave them open.

    Returns
    -------
    TYPE Boolean
        true if running in Jupyter or Ipython consoles within spyder.
        false otherwise (terminal)

    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True  # Jupyter notebook or qtconsole
        elif shell == "TerminalInteractiveShell":
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False  # Probably standard Python interpreter


def RT2fileName(param, referenceBarcode):
    """
    Finds the files in a list that contain the ReferenceBarcode in their name
    Also returs the ROI of each file in this list


    Parameters
    ----------
    param : class
        parameters class.
    referenceBarcode : string
        reference barcode name

    Returns
    -------
    fileNameReferenceList : list
        list of files with reference barcode in their name
    ROIList : list
        list of ROIs.

    """
    fileNameReferenceList = []
    ROIList = {}

    for file in param.fileList2Process:
        if referenceBarcode in file.split("_"):
            fileNameReferenceList.append(file)
            # ROIList[file] = os.path.basename(file).split("_")[positionROIinformation]
            fileParts=param.decodesFileParts(os.path.basename(file))
            ROIList[file] = fileParts['roi']
    return fileNameReferenceList, ROIList


def ROI2FiducialFileName(param, file, barcodeName):
    """
    Produces list of fiducial files that need to be loaded from a specific DAPI/barcode image


    Parameters
    ----------
    param : class
        parameters class.
    file : string
        filename.
    barcodeName : string
        barcode name to assign a fiducial to.

    Returns
    -------
    candidates : TYPE
        DESCRIPTION.

    """
    # gets rootFolder
    rootFolder = os.path.dirname(file)
    ROI = param.decodesFileParts(os.path.basename(file))['roi']
    
    channelFiducial = param.param["acquisition"]["fiducialBarcode_channel"]

    # looks for referenceFiducial file in folder
    listFiles = glob.glob(rootFolder + os.sep + "*.tif")

    candidates = [
        x
        for x in listFiles
        if (barcodeName+"_" in x)
        and (ROI == param.decodesFileParts(os.path.basename(x))['roi'])
        and (channelFiducial in os.path.basename(x))
    ]

    return candidates

def retrieveNumberUniqueBarcodesRootFolder(rootFolder, parameterFile, ext='tif'):
    """
    given a directory and a Parameter object, it returns the number of unique cycles/barcodes detected

    Parameters
    ----------
    rootFolder : string
    param : string
        parameterFile 
    ext : string, optional
        File extension. The default is 'tif'.

    Returns
    -------
    int
        number of unique cycles.

    """
    def unique(list1): 
        ''' function to get unique values'''
        # intilize a null list 
        unique_list = [] 
          
        # traverse for all elements 
        for x in list1: 
            # check if exists in unique_list or not 
            if x not in unique_list: 
                unique_list.append(x) 
    
        return unique_list      
            
    allFilesinRootFolder=glob.glob(rootFolder+os.sep+"*"+ext)
    
    param = Parameters(rootFolder, rootFolder+parameterFile)
       
    ROIs, RTs = [], []
    for x in allFilesinRootFolder:
        fileParts = param.decodesFileParts(x)
        ROIs.append(fileParts["roi"])
        RTs.append(fileParts["cycle"])
    
    numberUniqueCycles=len(unique(RTs))
        
    return numberUniqueCycles