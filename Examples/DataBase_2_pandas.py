#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:39:14 2020

@author: aguimera
"""

# import pickle
import pandas as pd
import PyGFETdb.DBSearch as DBs
import pickle

# %% Define search


# WafersList = ('B12752W2',
#                )

# CharTable = 'ACcharacts'
# Conditions = {'Wafers.name = ': WafersList,
#               'CharTable.Comments !=': ('Modulated', ),
#               'CharTable.MeasDate < ': ('2019-12-01 00:00:00', ),
#               'Trts.Name NOT LIKE ': ('%Col1', ),
#              }

# GroupBase = {'Conditions': Conditions,
#              'Table': CharTable,
#              'Last': True,
#              'GetGate': True,
#              }

WafersList = ('B12708W2',
              'B12752W2',
               )

CharTable = 'ACcharacts'
Conditions = {'Wafers.name = ': WafersList,
              'CharTable.Comments !=': ('Modulated', ),
              'CharTable.Comments  !=': ('py3', ),
              'CharTable.MeasDate < ': ('2019-12-01 00:00:00', ),
              'Trts.Name NOT LIKE ': ('%Col1', ),
              }

GroupBase = {'Conditions': Conditions,
             'Table': CharTable,
             'Last': True,
             'GetGate': True,
             }


# %% Get data from DataBase

Data, Trts = DBs.GetFromDB(**GroupBase)

# %% Create a pandas dataset
FileOut = 'DbData.pkl'

# Generate one pandas series for each characterization
pdSeries = []
for tn, dats in Data.items():
    for dd in dats:
        pdser = {}
        pdser['CharCl'] = dd
        pdser['IsOk'] = dd.IsOK
        pdser['Wafer'] = dd.Wafers['Name']
        pdser['Device'] = dd.Devices['Name']
        pdser['TrtName'] = tn
        pdser['Date'] = dd.GetDateTime()
        for k, v in dd.TrtTypes.items():
            pdser[k] = v
        for k, v in dd.Info.items():
            pdser[k] = v
        pdSeries.append(pdser)
            
# Create pandas dataframe
pdSeries = [pd.Series(ser) for ser in pdSeries]
dfRaw = pd.concat(pdSeries, axis=1).transpose()

# Save into pickle file
pickle.dump(dfRaw, open(FileOut, 'wb'))


