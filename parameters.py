import pandas as pd
import numpy as np
import json
import os.path
from os import path

class Parameters:
    def __init__(self):
        #self.parameters = param
        self.paramFile = "infoList.inf"
        self.param = {
                        'image': {
                                    'currentPlane': 1,
                                    # 'contrastMin': 0.1,
                                    # 'contrastMax': 0.9,
                                    'claheGridH': 8,
                                    'claheGridW': 8
                        },
                        'options': {
                                    'autoPlanes': True
                        },
                        'process': {
                                    'zmin': 1,
                                    'zmax': 59,
                                    'zwindows': 10,
                                    'windowSecurity': 2,
                                    'zProjectOption': 'sum'  # sum or MIP
                        }
        }


    def get_param(self, param=False):
        if not param:
            return self.param
        else:
            return self.param[param]

    def initializeStandardParameters(self):
        with open(self.paramFile, 'w') as f:
            json.dump(json.dumps(self.param), f, ensure_ascii=False)
            
    def loadParametersFile(self, fileName):
        if path.exists(fileName):
            with open(fileName) as json_file:
                data = json.load(json_file)
                
            self.param=json.loads(data)
            print("Parameters file read: {}".format(fileName))
        
        
    '''
    param = {
                        'image': {
                                    'currentPlane': 1,
                                    # 'contrastMin': 0.1,
                                    # 'contrastMax': 0.9,
                                    'claheGridH': 8,
                                    'claheGridW': 8
                        },
                        'options': {
                                    'autoPlanes': True
                        },
                        'process': {
                                    'zmin': 1,
                                    'zmax': 59,
                                    'zwindows': 10,
                                    'windowSecurity': 2,
                                    'zProjectOption': 'sum'  # sum or MIP
                        }
        }
    '''
    '''
    def read_params_file(self, rawDir):
        print(rawDir)
        pfile = pd.read_csv(rawDir+"/"+paramFile,
                            sep=' ',
                            header=None,
                            names=[0, 1, 2, 3])
        for i, row in pfile.iterrows():
            if np.isnan(row[2]):
                if row[1] == 'true':
                    self.parameters[row[0]] = True
                elif row[1] == 'false':
                    self.parameters[row[0]] = False
                else:
                    self.parameters[row[0]] = row[1]
            elif np.isnan(row[3]):
                self.parameters[row[0]] = [row[1]]+[row[2]]
            else:
                self.parameters[row[0]] = [row[1]]+[row[2]]+[row[3]]

        return self.parameters
    ''' 