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
import numpy as np

# %% Define search

WafersList = ('B13044W1',
               )

CharTable = 'DCcharacts'
Conditions = {'Wafers.name = ': WafersList,
              }

GroupBase = {'Conditions': Conditions,
             'Table': CharTable,
             'Last': False,
             'GetGate': True,
             }


# %% Get data from DataBase

Data, Trts = DBs.GetFromDB(**GroupBase)

# %% Create a pandas dataset
FileOut = 'DbRaw.pkl'

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
            if k == 'Name':
                pdser['TrtName'] = v
            else:
                pdser[k] = v
        for k, v in dd.Info.items():
            if k ==  'Gate_id':
                continue
            pdser[k] = v
        pdSeries.append(pdser)
            
# Create pandas dataframe
pdSeries = [pd.Series(ser) for ser in pdSeries]
dfRaw = pd.concat(pdSeries, axis=1).transpose()

# %% Formating and saving
# Data types definition for each column
# default type is 'category'
# Genarates a log xls file with dataset description

DataTypes = {}
for col in dfRaw.keys():    
    if col in ('Width', 'Length', 'Pass', 'Area', 'Ph', 'IonStrength', 'AnalyteCon'):
        DataTypes[col] = np.float
    else:
        DataTypes[col] = 'category'

DataTypes['CharCl'] = object
DataTypes['IsOk'] = bool

dfRaw = dfRaw.astype(DataTypes)

# Save into pickle file
pickle.dump(dfRaw, open(FileOut, 'wb'))

# Generate output summary log
dfRaw.describe(include=['category', ]).to_excel('Log_db.xls')
Desc = dfRaw.describe(exclude=[object, ])
Desc.to_excel('Log1_db.xls')