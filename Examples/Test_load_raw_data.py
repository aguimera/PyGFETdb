# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 23:15:19 2016

@author: aguimera
"""

import matplotlib.pyplot as plt
import gc
import PyGFETdb.DBCore as Pydb
import PyGFETdb.DataStructures as PyFETData
import PyGFETdb.AnalyzeData as PyFETAnalyze
from PyGFETdb.DataClass import DataCharAC, DataCharDC
import glob
import pickle
import pandas as pd
import numpy as np
from PyGFETdb.DBInterface import Data2Pandas, ClassQueries, pdAttr
import quantities as pq


filen='B12744W3-Xip6NS-testN.pkl'


DataIn = pickle.load(open(filen, 'rb'), encoding='latin')
DevDCVals = DataIn['DevDC']
DevACVals = DataIn['DevAC']

DictClDC = {}
DictClAC = {}

for chn, dat in DevDCVals.items():
    if chn == 'Gate':
        continue
    DCclass = DataCharDC(dat)
    DictClDC[chn] = (DCclass, )
    acDat = DCclass.__dict__.copy()
    acDat.update(DevACVals[chn])
    acDat.pop('Ids')
    acDat['Vgs'] = acDat.pop('VgsAC')
    acDat['Vds'] = acDat.pop('VdsAC')
    DictClAC[chn] = (DataCharAC(acDat), )

# Char = DictClAC['Ch01'][0]
# Char.Get(**ClassQueries['Vrms01'])

# Char.Get(**{'Param': 'Vrms',
#             'Vgs': -0.1*pq.V,
#             'Ud0Norm': True,
#             # 'Units': 'uV',
#             })

dbRaw = Data2Pandas(DictClAC)
    
PdSeries = []
for index, row in dbRaw.iterrows():
    char = row['CharCl']
    Vds = char.GetVds()
    for vds in Vds:
        vals = {}
        vals['Vds'] = vds.flatten()
        for parn, park in ClassQueries.items():
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

DataTypes = dbRaw.dtypes
for col in pdAttr['ScalarCols']:
    DataTypes[col] = np.float
for col in pdAttr['ArrayCols']:
    DataTypes[col] = object
dfDat = dfDat.astype(DataTypes)
dfDat.attrs.update(pdAttr)

#%% Plots
import math

dfDat.query("IsOk == True", inplace=True)

Pars = ('Ids','GM', 'NoA', 'NoB', 'NoC', 'Irms', 'Vrms')
xvar = dfDat.attrs['Vgs']
ArVals = ClassQueries

nRows = 2
nCols = math.ceil(len(Pars)/2)
fig, Axs = plt.subplots(nrows=nRows, ncols=nCols, sharex=True)
Axs = Axs.flatten()

for ip, par in enumerate(Pars):
    Vals = np.array([])
    for index, row in dfDat.iterrows():
        v = row[par]
        Vals = np.vstack((Vals, v)) if Vals.size else v
    Vals = Vals.magnitude.transpose()
    Axs[ip].plot(xvar, Vals, 'k', alpha=0.2)
    mean = np.nanmean(Vals, axis=1)
    std = np.nanstd(Vals, axis=1)
    Axs[ip].plot(xvar, mean, 'r', lw=1.5)
    Axs[ip].fill_between(xvar, mean+std, mean-std, color='r', alpha=0.2)
    if 'Units' in ArVals[par]:
        Axs[ip].set_ylabel(par + '[' + ArVals[par]['Units'] + ']')
    else:
        Axs[ip].set_ylabel(par)
#%%
from matplotlib.backends.backend_pdf import PdfPages

PDF = PdfPages('multipage_pdf.pdf')

for index, row in dfDat.iterrows():
    fig, (axM, axP) = plt.subplots(2,1)
    char = row['CharCl']    
    axM.loglog(char.GetFgm(), char.GetGmMag()/2)
    axM.plot(char.GetFgm()[0], np.abs(char.GetGM()), 'k*')
    axP.semilogx(char.GetFgm(), char.GetGmPh())
    PDF.savefig(fig)
    plt.close('all')
    
PDF.close()



