#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 16:14:13 2020

@author: aguimera
"""

import pickle
import pandas as pd
import seaborn as sns
import quantities as pq
import numpy as np

# load dataset from file
FileIn = 'DbData.pkl'
dfRaw = pickle.load(open(FileIn, 'rb'))

# %% Define the desidered values to calculate

# Array values on measured Vgs range
Vgs = np.linspace(-0.1, 0.5, 100) * pq.V
ArrayValues = {
              'Ids': {'Param': 'Ids',
                      'Vds': None,
                      'Vgs': Vgs,
                      'Ud0Norm': False,
                      'Units': 'uA',
                      },
              'GM': {'Param': 'GM',
                      'Vgs': Vgs,
                      'Vds': None,
                      'Ud0Norm': False,
                      'Units': 'mS',
                      },
              'GMV': {'Param': 'GMV',
                      'Vgs': Vgs,
                      'Vds': None,
                      'Ud0Norm': False,
                      'Units': 'mS/V'
                      },
              'Irms': {'Param': 'Irms',
                        'Vds': None,
                        'Vgs': Vgs,
                        'Ud0Norm': False,
                        'Units': 'uA',
                        },
              'Vrms': {'Param': 'Vrms',
                        'Vds': None,
                        'Vgs': Vgs,
                        'Ud0Norm': False,
                        'Units': 'uV'
                        },
              'NoA': {'Param': 'NoA',
                      'Vds': None,
                      'Vgs': Vgs,
                      'Ud0Norm': False,
                      # !!!!! First call to noise, force fitting
                      'FFmin': 5,
                      'FFmax': 7e3,
                      'Units': 'A**2',
                      },
              'NoB': {'Param': 'NoB',
                      'Vds': None,
                      'Vgs': Vgs,
                      'Ud0Norm': False,
                      },
              'NoC': {'Param': 'NoC',
                      'Vds': None,
                      'Vgs': Vgs,
                      'Ud0Norm': False,
                      'Units': 'A**2',
                      },
              }

VgsNorm = np.linspace(-0.4, 0.3, 100) * pq.V
ArrayValuesNorm = {
              'IdsNorm': {'Param': 'Ids',
                          'Vds': None,
                          'Vgs': VgsNorm,
                          'Ud0Norm': True,
                          'Units': 'uA',
                          },
              'GMNorm': {'Param': 'GM',
                          'Vgs': VgsNorm,
                          'Vds': None,
                          'Ud0Norm': True,
                          'Units': 'mS',
                          },
              'GMVNorm': {'Param': 'GMV',
                          'Vgs': VgsNorm,
                          'Vds': None,
                          'Ud0Norm': True,
                          'Units': 'mS/V'
                          },
              'IrmsNorm': {'Param': 'Irms',
                            'Vds': None,
                            'Vgs': VgsNorm,
                            'Ud0Norm': True,
                            'Units': 'uA',
                            },
              'VrmsNorm': {'Param': 'Vrms',
                            'Vds': None,
                            'Vgs': VgsNorm,
                            'Ud0Norm': True,
                            'Units': 'uV'
                            },
              'NoANorm': {'Param': 'NoA',
                          'Vds': None,
                          'Vgs': VgsNorm,
                          'Ud0Norm': True,
                          # !!!!! First call to noise, force fitting
                          # 'FFmin': 5,
                          # 'FFmax': 7e3,
                          'Units': 'A**2',
                          },
              'NoBNorm': {'Param': 'NoB',
                          'Vds': None,
                          'Vgs': VgsNorm,
                          'Ud0Norm': True,
                          },
              'NoCNorm': {'Param': 'NoC',
                          'Vds': None,
                          'Vgs': VgsNorm,
                          'Ud0Norm': True,
                          'Units': 'A**2',
                          },
              }

# Scalar values
ScalarValues = {'CNP': {'Param': 'Ud0',
                        'Units': 'mV'
                        },
                'Ids01': {'Param': 'Ids',
                          'Vgs': -0.1*pq.V,
                          'Vds': None,
                          'Ud0Norm': True,
                          'Units': 'uA'
                          },
                'GM01': {'Param': 'GMV',
                          'Vgs': -0.1*pq.V,
                          'Vds': None,
                          'Ud0Norm': True,
                          'Units': 'mS/V',
                          },
                'Vrms01': {'Param': 'Vrms',
                            'Vgs': -0.1*pq.V,
                            'Vds': None,
                            'Ud0Norm': True,
                            'Units': 'uV',
                            },
                'IdsCNP': {'Param': 'Ids',
                            'Vgs': 0*pq.V,
                            'Vds': None,
                            'Ud0Norm': True,
                            'Units': 'uA'
                            },
                'NoA01': {'Param': 'NoA',
                          'Vgs': -0.1*pq.V,
                          'Vds': None,
                          'Ud0Norm': True,
                          'Units': 'A**2',
                            },
                'NoB01': {'Param': 'NoB',
                          'Vgs': -0.1*pq.V,
                          'Vds': None,
                          'Ud0Norm': True,
                          # 'Units': 'A**2',
                            },
                'NoC01': {'Param': 'NoC',
                          'Vgs': -0.1*pq.V,
                          'Vds': None,
                          'Ud0Norm': True,
                          'Units': 'A**2',
                            },

                }


# %% add values to dataframe

FileOut = 'DbData2.pkl'

Parameters = ArrayValues.copy()
Parameters.update(ScalarValues)
Parameters.update(ArrayValuesNorm)

df1 = dfRaw.copy()

PdSeries = []
for index, row in df1.iterrows():
    char = row['CharCl']
    vals = {}
    for parn, park in Parameters.items():
        try:
            val = char.Get(**park)
        except:
            print('Error', parn, char.Name)
            continue
        if val is None:
            continue
        if val.size > 1:
            vals[parn] = val
        else:
            vals[parn] = val.flatten()
    PdSeries.append(row.append(pd.Series(vals)))

dfDat = pd.concat(PdSeries, axis=1).transpose()

ScalarTypes = {}
for p in ScalarValues:
    ScalarTypes[p] = np.float
for p in ('Width', 'Length', 'Pass', 'Area', 'Ph', 'IonStrength'):
    ScalarTypes[p] = np.float
dfDat = dfDat.astype(ScalarTypes)

# adding extra variables for easy performance
ExtraVars = (ScalarValues, ArrayValues, ArrayValuesNorm, Vgs, VgsNorm)
pickle.dump((dfDat, ExtraVars), open(FileOut, 'wb'))

