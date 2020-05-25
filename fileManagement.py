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


# =============================================================================
# CLASSES
# =============================================================================


class log:
    def __init__(self, rootFolder="./", fileNameRoot="HiM_analysis"):
        now = datetime.now()
        dateTime = now.strftime("%d%m%Y_%H%M%S")
        self.fileName = rootFolder + os.sep + fileNameRoot + dateTime + ".log"
        self.fileNameMD = self.fileName.split(".")[0] + ".md"

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
        print(self.getFullString(text, status))
        self.save("\n" + text, status)

    # returns formatted line to be outputed
    def getFullString(self, text="", status="info"):
        now = datetime.now()
        return "{}|{}>{}".format(now.strftime("%d/%m/%Y %H:%M:%S"), status, text)

    def addSimpleText(self, title):
        print("{}".format(title))
        # with open(self.fileName, 'a') as file:
        #    file.write(title)
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
            if os.path.isdir(folder)
            and len(glob.glob(folder + os.sep + "*." + extension)) > 0
        ]
        # os.path.name(folder)[0]!='F']
        if len(hfolders) > 0:
            self.listFolders = hfolders
        else:
            self.listFolders = []

        # checks if there are files with the required extension in the root folder provided
        if (
            os.path.isdir(self.masterFolder)
            and len(glob.glob(self.masterFolder + os.sep + "*." + extension)) > 0
        ):
            # self.listFolders=self.masterFolder
            self.listFolders.append(self.masterFolder)

        print("Detected {} folders with images".format(len(self.listFolders)))

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
        # self.zProjectFolder=filesFolder+os.sep+param.param['zProject']['folder']
        # self.createSingleFolder(self.zProjectFolder)

        self.outputFolders["zProject"] = (
            filesFolder + os.sep + param.param["zProject"]["folder"]
        )
        self.outputFolders["alignImages"] = (
            filesFolder + os.sep + param.param["alignImages"]["folder"]
        )
        self.outputFolders["segmentedObjects"] = (
            filesFolder + os.sep + param.param["segmentedObjects"]["folder"]
        )
        self.outputFolders["buildsPWDmatrix"] = filesFolder + os.sep + "buildsPWDmatrix"
        self.outputFolders["projectsBarcodes"] = (
            filesFolder + os.sep + param.param["projectsBarcodes"]["folder"]
        )

        self.createSingleFolder(self.outputFolders["zProject"])
        self.createSingleFolder(self.outputFolders["alignImages"])
        self.createSingleFolder(self.outputFolders["segmentedObjects"])
        self.createSingleFolder(self.outputFolders["buildsPWDmatrix"])
        self.createSingleFolder(self.outputFolders["projectsBarcodes"])

        # self.outputFiles['zProject']=self.outputFolders['zProject']+os.sep+param.param['zProject']['outputFile']
        self.outputFiles["alignImages"] = (
            self.outputFolders["alignImages"]
            + os.sep
            + param.param["alignImages"]["outputFile"]
        )
        self.outputFiles["dictShifts"] = (
            self.masterFolder + os.sep + param.param["alignImages"]["outputFile"]
        )
        self.outputFiles["segmentedObjects"] = (
            self.outputFolders["segmentedObjects"]
            + os.sep
            + param.param["segmentedObjects"]["outputFile"]
        )
        self.outputFiles["buildsPWDmatrix"] = (
            self.outputFolders["buildsPWDmatrix"] + os.sep + "buildsPWDmatrix"
        )
        self.outputFiles["projectsBarcodes"] = (
            self.outputFolders["projectsBarcodes"]
            + os.sep
            + param.param["projectsBarcodes"]["outputFile"]
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
            "image": {
                "currentPlane": 1,
                # 'contrastMin': 0.1,
                # 'contrastMax': 0.9,
                "claheGridH": 8,
                "claheGridW": 8,
            },
            "acquisition": {
                "label": "DAPI",  # barcode, fiducial
                "positionROIinformation": 3,
            },
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
                "background_method": "inhomogeneous",  # flat or inhomogeneous
                "background_sigma": 3.0,  # used to remove inhom background
                "threshold_over_std": 1.0,  # threshold used to detect sources
                "fwhm": 3.0,  # source size in px
                "brightest": 1100,  # max number of objects segmented per FOV
                "intensity_min": 0,  # min int to keep object
                "intensity_max": 59,  # max int to keeep object
                "area_min": 50,  # min area to keeep object
                "area_max": 500,  # max area to keeep object
            },
        }
        self.initializeStandardParameters()
        self.paramFile = rootFolder + os.sep + label
        self.loadParametersFile(self.paramFile)
        self.param["rootFolder"] = rootFolder

    def get_param(self, param=False):
        if not param:
            return self.param
        else:
            return self.param[param]

    def initializeStandardParameters(self):
        with open(self.paramFile, "w") as f:
            # json.dump(json.dumps(self.param), f, ensure_ascii=False, indent=4)
            json.dump(self.param, f, ensure_ascii=False, sort_keys=True, indent=4)
        print(
            "Model parameters file saved to: {}".format(
                os.getcwd() + os.sep + self.paramFile
            )
        )

    def loadParametersFile(self, fileName):
        if path.exists(fileName):
            with open(fileName) as json_file:
                self.param = json.load(json_file)

            print("Parameters file read: {}".format(fileName))

    # method returns label specific filenames from filename list
    def files2Process(self, filesFolder):

        # finds if there is 2 or 3 channels for DAPI
        fileList2Process = [
            file
            for file in filesFolder
            if file.split("_")[-1].split(".")[0] == "ch02" and "DAPI" in file.split("_")
        ]

        if len(fileList2Process) > 0:
            channelDAPI_fiducial = "ch02"
            channelDAPI_RNA = "ch01"
        else:
            channelDAPI_fiducial = "ch01"
            channelDAPI_RNA = "ch04"

        # selects DAPI files
        if self.param["acquisition"]["label"] == "DAPI":
            self.fileList2Process = [
                file
                for file in filesFolder
                if file.split("_")[-1].split(".")[0] == "ch00"
                and "DAPI" in file.split("_")
            ]

        # selects DAPIch2 files
        if self.param["acquisition"]["label"] == "RNA":
            self.fileList2Process = [
                file
                for file in filesFolder
                if file.split("_")[-1].split(".")[0] == channelDAPI_RNA
                and "DAPI" in file.split("_")
            ]

        # selects barcode files
        elif self.param["acquisition"]["label"] == "barcode":
            self.fileList2Process = [
                file
                for file in filesFolder
                if len([i for i in file.split("_") if "RT" in i]) > 0
                and file.split("_")[-1].split(".")[0] == "ch01"
            ]

        # selects fiducial files
        elif self.param["acquisition"]["label"] == "fiducial":
            self.fileList2Process = [
                file
                for file in filesFolder
                if (
                    len([i for i in file.split("_") if "RT" in i]) > 0
                    and file.split("_")[-1].split(".")[0] == "ch00"
                )
                or (
                    "DAPI" in file.split("_")
                    and file.split("_")[-1].split(".")[0] == channelDAPI_fiducial
                )
            ]


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


def RT2fileName(fileList2Process, referenceBarcode, positionROIinformation=3):
    """
    Finds the files in a list that contain the ReferenceBarcode in their name
    Also returs the ROI of each file in this list

    Parameters
    ----------
    fileList2Process : TYPE
        DESCRIPTION.
    referenceBarcode : TYPE
        DESCRIPTION.
    positionROIinformation : TYPE, optional
        DESCRIPTION. The default is 3.

    Returns
    -------
    fileNameReferenceList : TYPE
        List of filenames with referenceBarcode in their RT field
    ROIList : TYPE
        Dictionary of ROIs for each file in list

    """
    fileNameReferenceList = []
    ROIList = {}

    for file in fileList2Process:
        if referenceBarcode in file.split("_"):
            fileNameReferenceList.append(file)
            ROIList[file] = os.path.basename(file).split("_")[positionROIinformation]

    return fileNameReferenceList, ROIList
