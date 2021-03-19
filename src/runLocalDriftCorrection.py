#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:48:55 2020

@author: marcnol
"""

""
# =============================================================================
# IMPORTS
# =============================================================================

from datetime import datetime

from fileProcessing.fileManagement import Parameters

# from imageProcessing.localDriftCorrection import localDriftCorrection
from fileProcessing.functionCaller import HiMfunctionCaller,HiM_parseArguments

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    runParameters=HiM_parseArguments()    
    runParameters["localAlignment"] =True
    
    HiM = HiMfunctionCaller(runParameters, sessionName="localDriftCorrection")
    HiM.initialize()
    session1, log1=HiM.session1, HiM.log1
    
    HiM.lauchDaskScheduler(threadsRequested = runParameters["threads"],maximumLoad=0.8)

    for ilabel in range(len(HiM.labels2Process)):
        HiM.log1.addSimpleText("**Analyzing label: {}**".format(HiM.labels2Process[ilabel]["label"]))

        # sets parameters
        param = Parameters(runParameters["rootFolder"], HiM.labels2Process[ilabel]["parameterFile"])
        param.param['parallel']=HiM.parallel

        # [local drift correction]
        HiM.localDriftCorrection(param, ilabel)

        print("\n")
        del param

    # exits
    HiM.session1.save(HiM.log1)
    HiM.log1.addSimpleText("\n===================={}====================\n".format("Normal termination"))

    if runParameters["parallel"]:
        HiM.client.close()
        HiM.cluster.close()   

    del HiM
    
    print("Elapsed time: {}".format(datetime.now() - begin_time))

    # runParameters=parseArguments()    

    # print("parameters> rootFolder: {}".format(runParameters["rootFolder"]))
    # sessionName = "localDriftCorrection"

    # labels2Process = [
    #     {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
    #     {"label": "barcode", "parameterFile": "infoList_barcode.json"},
    #     {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
    # ]

    # # session
    # session1 = session(runParameters["rootFolder"], sessionName)

    # # setup logs
    # log1 = log(rootFolder=runParameters["rootFolder"],parallel=runParameters["parallel"])
    # # labels2Process indeces: 0 fiducial, 1:
    # labelParameterFile = labels2Process[2]["parameterFile"]
    # param = Parameters(runParameters["rootFolder"], labelParameterFile)

    # dataFolder = folders(param.param["rootFolder"])

    # if runParameters["parallel"]:
    #     daskClusterInstance = daskCluster(20)
    #     print("Go to http://localhost:8787/status for information on progress...")
        
    #     cluster = LocalCluster(n_workers=daskClusterInstance.nThreads,
    #                         # processes=True,
    #                         # threads_per_worker=1,
    #                         # memory_limit='2GB',
    #                         # ip='tcp://localhost:8787',
    #                         ) 
    #     client = Client(cluster)
    
    # for currentFolder in dataFolder.listFolders:
    #     filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
    #     dataFolder.createsFolders(currentFolder, param)

    #     # generates lists of files to process
    #     param.files2Process(filesFolder)
    #     if runParameters["parallel"]:
    #         param.param['parallel']=True
    #     else:
    #         param.param['parallel']=False
            
    #     for fileName in param.fileList2Process:
    #         session1.add(fileName, sessionName)

    #     result = client.submit(localDriftCorrection,param, log1, session1)
    #     errorCode, dictShift, alignmentResultsTable  = client.gather(result)

    # if errorCode != 0:
    #     print("Error code reported: {}".format(errorCode))
    # else:
    #     print("normal termination")

    # if runParameters["parallel"]:
    #     client.close()
    #     cluster.close()   
        
    # print("Elapsed time: {}".format(datetime.now() - begin_time))
