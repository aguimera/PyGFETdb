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
FileIn = 'DbRaw.pkl'
dfRaw = pickle.load(open(FileIn, 'rb'))

# %% Define the desidered values to calculate
# Vds has no to be difined, one row is created for each Vds value found

# Array values on measured Vgs range
Vgs = np.linspace(-0.1, 0.5, 100) * pq.V
ArrayValues = {
              'Ids': {'Param': 'Ids',
                      'Vgs': Vgs,
                      'Ud0Norm': False,
                      'Units': 'uA',
                      },
              'GM': {'Param': 'GM',
                     'Vgs': Vgs,
                     'Ud0Norm': False,
                     'Units': 'mS',
                     },
              'GMV': {'Param': 'GMV',
                      'Vgs': Vgs,
                      'Ud0Norm': False,
                      'Units': 'mS/V'
                      },
              }

VgsNorm = np.linspace(-0.4, 0.3, 100) * pq.V
ArrayValuesNorm = {
              'IdsNorm': {'Param': 'Ids',
                          'Vgs': VgsNorm,
                          'Ud0Norm': True,
                          'Units': 'uA',
                          },
              'GMNorm': {'Param': 'GM',
                         'Vgs': VgsNorm,
                         'Ud0Norm': True,
                         'Units': 'mS',
                         },
              'GMVNorm': {'Param': 'GMV',
                          'Vgs': VgsNorm,
                          'Ud0Norm': True,
                          'Units': 'mS/V'
                          },
              }

# Scalar values
ScalarValues = {'CNP': {'Param': 'Ud0',
                        'Units': 'mV'
                        },
                'Ids01': {'Param': 'Ids',
                          'Vgs': -0.1*pq.V,
                          'Ud0Norm': True,
                          'Units': 'uA'
                          },
                'GM01': {'Param': 'GMV',
                         'Vgs': -0.1*pq.V,
                         'Ud0Norm': True,
                         'Units': 'mS/V',
                         },
                }


# %% add values to dataframe

FileOut = 'DbData.pkl'

Parameters = ArrayValues.copy()
Parameters.update(ScalarValues)
Parameters.update(ArrayValuesNorm)

df1 = dfRaw.copy()

PdSeries = []
for index, row in df1.iterrows():
    char = row['CharCl']
    Vds = char.GetVds()
    for vds in Vds:
        vals = {}
        vals['Vds'] = vds.flatten()
        for parn, park in Parameters.items():
            park['Vds'] = vds
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

# %% Formating and saving
# Data types definition for each column
# default type is 'category'

DataTypes = {}
for col in dfRaw.keys():    
    if col in ('Width', 'Length', 'Pass', 'Area', 'Ph', 'IonStrength', 'AnalyteCon'):
        DataTypes[col] = np.float
    else:
        DataTypes[col] = 'category'

DataTypes['CharCl'] = object
DataTypes['Vds'] = np.float
DataTypes['IsOk'] = bool

for p in ArrayValues:
    DataTypes[p] = object
for p in ArrayValuesNorm:
    DataTypes[p] = object

for p in ScalarValues:
    DataTypes[p] = np.float

dfDat = dfDat.astype(DataTypes)

# Save dtaframe to file
# adding extra variables for easy performance
ExtraVars = (ScalarValues, ArrayValues, ArrayValuesNorm, Vgs, VgsNorm)
pickle.dump((dfDat, ExtraVars), open(FileOut, 'wb'))

# Generate output summary log
dfDat.describe(include=['category', ]).to_excel('Log.xls')
Desc = dfDat.describe(exclude=[object, ])
Desc.to_excel('Log1.xls')





