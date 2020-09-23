#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:31:14 2020

@author: marcnol
"""

# =============================================================================
# IMPORTS
# =============================================================================

from datetime import datetime
from fileProcessing.fileManagement import Parameters
from fileProcessing.functionCaller import HiMfunctionCaller, HiM_parseArguments

# to remove in a future version
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    begin_time = datetime.now()

    runParameters=HiM_parseArguments()    

    HiM = HiMfunctionCaller(runParameters, sessionName="SegmentsMasks")
    HiM.initialize()
    session1, log1=HiM.session1, HiM.log1
    
    HiM.lauchDaskScheduler()

    for ilabel in range(len(HiM.labels2Process)):
        HiM.log1.addSimpleText("**Analyzing label: {}**".format(HiM.labels2Process[ilabel]["label"]))

        # sets parameters
        param = Parameters(runParameters["rootFolder"], HiM.labels2Process[ilabel]["parameterFile"])
        param.param['parallel']=HiM.parallel
        
        # [segments DAPI and spot masks]
        HiM.segmentMasks(param, ilabel)

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
