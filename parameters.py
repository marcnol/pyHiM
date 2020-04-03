import pandas as pd
import numpy as np

paramFile = "infoList.inf"

class Parameters:
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

    def __init__(self):
        self.parameters = {}

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

    def get_param(self, param=False):
        if not param:
            return self.param
        else:
            return self.param[param]
