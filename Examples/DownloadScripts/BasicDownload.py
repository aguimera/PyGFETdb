#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 14:41:57 2022

@author: aguimera
"""
from PyGFETdb.DBInterface import GetFromDB, Data2Pandas, CalcElectricalParams
from PyGFETdb.DBInterface import ClassQueries, pdAttr
import pandas as pd

DevicesList = ('B12744W3-Xip6NS',
               )

CharTable = 'ACcharacts'
Conditions = {'Devices.name = ': DevicesList,
              }

GroupBase = {'Conditions': Conditions,
             'Table': CharTable,
             'Last': False,
             'GetGate': True,
             }


# %% Get data from DataBase

Data, Trts = GetFromDB(**GroupBase)
dfRaw = Data2Pandas(Data)

dfRaw.to_pickle('RawDat.pkl')

# %% Calculate basic electrical parameters
"""
Basic definitions are implemented in DBInterface module

ClassQueries -> default electrical parameters to calculate

ClsQueries --> dictiory where key is the column Name, 
                              value is the kwargs for Get function

pdAttr --> metadata added as attrr in pandas dataset

Example
Vgs Vector definitions
Vgs = np.linspace(-0.1, 0.5, 100) * pq.V
VgsNorm = np.linspace(-0.4, 0.3, 100) * pq.V

Queries definition
ArrayQueries = {'Ids': {'Param': 'Ids',   ### Vgs vector
                        'Vds': None,
                        'Vgs': Vgs,
                        'Ud0Norm': False,
                        'Units': 'uA',
                        },
                'IdsNorm': {'Param': 'Ids', ### CNP norm vector
                            'Vds': None,
                            'Vgs': VgsNorm,
                            'Ud0Norm': True,
                            'Units': 'uA',
                            },
                'Ids01': {'Param': 'Ids', ### Scalar value
                          'Vgs': -0.1*pq.V,
                          'Vds': None,
                          'Ud0Norm': True,
                          'Units': 'uA'
                          },

Atributtes definition, recomended for further analysis consistency
pdAttr = {'Vgs': Vgs,
          'VgsNorm': VgsNorm,
          'ScalarCols': list of scalar columns, 
          'ArrayCols': list of vector columns,
          }

"""

dfRaw = pd.read_pickle('RawDat.pkl')
dfDat = CalcElectricalParams(dbRaw=dfRaw,
                             ClsQueries=ClassQueries,
                             dfAttr=pdAttr)
                             
dfDat.to_pickle('Data.pkl')


    
    
