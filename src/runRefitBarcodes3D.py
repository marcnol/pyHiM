#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 12:50:55 2020

@author: marcnol
"""

from datetime import datetime
from fileProcessing.fileManagement import Parameters
from fileProcessing.functionCaller import HiMfunctionCaller, HiM_parseArguments

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    runParameters=HiM_parseArguments()    
    runParameters["refit"]=True
    
    HiM = HiMfunctionCaller(runParameters, sessionName="refitBarcodes3D")
    HiM.initialize()
    session1, log1=HiM.session1, HiM.log1
    
    HiM.lauchDaskScheduler()

    for ilabel in range(len(HiM.labels2Process)):
        HiM.log1.addSimpleText("**Analyzing label: {}**".format(HiM.labels2Process[ilabel]["label"]))

        # sets parameters
        param = Parameters(runParameters["rootFolder"], HiM.labels2Process[ilabel]["parameterFile"])
        param.param['parallel']=HiM.parallel

        # [refits spots in 3D]
        HiM.refitBarcodes(param, ilabel)

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
        