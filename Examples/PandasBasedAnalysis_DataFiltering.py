#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 18:10:56 2020

@author: aguimera
"""

import pickle
import pandas as pd
import seaborn as sns
import quantities as pq
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.artist import ArtistInspector
import matplotlib.colors as mcolors
import itertools



def UpdateTreeDictProp(obj, prop):
    ains = ArtistInspector(obj)
    validp = ains.get_setters()
    for p in prop.keys():
        if p in validp:
            obj.set(**{p: prop[p]})
        else:
            obj2 = getattr(obj, 'get_' + p)()
            UpdateTreeDictProp(obj2, prop[p])
            

# load dataset from file
FileIn = 'DbData2.pkl'
(dfDat, ExtraVars) = pickle.load(open(FileIn, 'rb'))
(ScalarValues, ArrayValues, ArrayValuesNorm, Vgs, VgsNorm) = ExtraVars

# %% Data filtering and plotting multiple vars
# Boxplots
# Copy data
dSel = dfDat.copy()

# Define queries
queries = (
            "CNP > 200 ",
            "Ids01 < 30 ",
            "NoC01 < 1e-21 ",
           )
for query in queries:
    dSel.query(query, inplace=True)

# generate figure
LogAxis = ('Irms01', 'Vrms01', 'NoA01', 'NoC01')
nRows = 2
nCols = math.ceil(len(ScalarValues.keys())/nRows)
fig, Axs = plt.subplots(nRows, nCols)
Axs = Axs.flatten()

for ic, Par in enumerate(ScalarValues.keys()):
    ax = Axs[ic]
    bg = sns.boxplot(x='Wafer',
                     y=Par,
                     hue='Width',
                     ax=ax,
                     data=dSel,
                     )
    bg.set_xticklabels(bg.get_xticklabels(), rotation=20)
    bg.legend(title='Width', loc='upper left', fontsize='xx-small')
    if Par in LogAxis:
        bg.set_yscale('log')


#%% values correlations

# Copy data
dSel = dfDat.copy()

# Define queries
queries = (
            "CNP > 200 ",
            "Ids01 < 30 ",            
           )
for query in queries:
    dSel.query(query, inplace=True)

# Calculate logvalues
dSel['LogVrms'] = np.log10(dSel['Vrms01'])
dSel['LogNoA'] = np.log10(dSel['NoA01'])
dSel['LogNoC'] = np.log10(dSel['NoC01'])

sns.pairplot(dSel,
             vars=('LogVrms', 'LogNoA', 'LogNoC', 'Ids01', 'GM01'),
             hue="Wafer")

#%% Data filtering and plotting multiple array variables
# Copy data
dSel = dfDat.copy()

# Define queries
queries = (
            "CNP > 200 ",
            "Ids01 < 30 ",
            "NoC01 < 1e-21 ",
            "Wafer == 'B12708W2' ",
            "IsOk == True"
           )
for query in queries:
    dSel.query(query, inplace=True)

xlab = 'Vgs [V]'
Vals2Plot = ArrayValues
xVals = Vgs
LogAxis = ('Irms', 'Vrms', 'NoA', 'NoC')

AxisProp = {}
for parn, park in Vals2Plot.items():   
    AxProp = {
              'autoscaley_on': True,
              'xlabel': xlab,
              }
    if 'Units' in park:
        AxLabel = '{} [{}]'.format(park['Param'], park['Units'])
    else:
        AxLabel = park['Param']
    AxProp['ylabel'] = AxLabel
    if park['Param'] in LogAxis:
        AxProp['yscale'] = 'log'
    AxisProp[parn] = AxProp

nRows = 2
nCols = math.ceil(len(Vals2Plot.keys())/nRows)
fig, Axs = plt.subplots(nRows, nCols)
Axs = Axs.flatten()
LoopColors = itertools.cycle(mcolors.TABLEAU_COLORS)

dfg = dSel[dfDat.IsOk == True].groupby('Device')
for ic, ParN in enumerate(Vals2Plot.keys()):
    Ax = Axs[ic]
    UpdateTreeDictProp(Ax, AxisProp[ParN])
    for ic, g in enumerate(dfg.groups):
        col = next(LoopColors)
        gg = dfg.get_group(g)
        Val = np.array([])
        for v in gg[ParN]:
            try:
                Val = np.vstack((Val, v)) if Val.size else v
            except:
                print(v, gg.TrtName)
        Val = Val.magnitude.transpose()
        Mean = np.nanmean(Val, axis=1)
        Std = np.nanstd(Val, axis=1)
        Ax.plot(xVals, Val, color=col, alpha=0.01)
        Ax.plot(xVals, Mean, color=col, lw=2, label=g)
        Ax.fill_between(xVals, Mean-Std, Mean+Std, color=col, alpha=0.2)
    # Ax.legend(fontsize='small')
     
#%% Look for thermal noise
# Copy data

dSel = dfDat.copy()
queries = (
            "NoC01 > 1e-22 ",
            "Device != 'B12708W2-M3'"
           )
for query in queries:
    dSel.query(query, inplace=True)

dgroups = dSel.groupby('TrtName')

fig, Ax = plt.subplots()
# LoopColors = itertools.cycle(mcolors.CSS4_COLORS)
LoopColors = itertools.cycle(mcolors.TABLEAU_COLORS)


for ic, gn in enumerate(dgroups.groups):
    gg = dgroups.get_group(gn)    
    
    CharCl = gg.CharCl.values[0]
    fpsd = CharCl.GetFpsd()
    PSD = CharCl.GetPSD()

    col=next(LoopColors)
    Ax.loglog(fpsd, PSD, color=col, alpha=0.05)
    Ax.loglog(fpsd, np.mean(PSD, axis=1), lw=1.5, color=col, label=gn)
    
Ax.legend(fontsize='x-small')





