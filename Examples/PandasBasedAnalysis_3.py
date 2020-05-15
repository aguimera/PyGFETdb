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


xlab = 'Vgs [V]'
Vals2Plot = ArrayValues
xVals = Vgs
LogAxis = ('Irms', 'Vrms', 'NoA', 'NoC')


# ddg = df.groupby('DevName')
# # dplt = ddg.get_group('B10803W21-Xip1')
# # dplt = ddg.get_group('B10803W17-Xip4')
# dplt = ddg.get_group('B10803W16-Xip5')



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

dfg = dfDat[dfDat.IsOk == True].groupby('Device')

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
    Ax.legend(fontsize='small')
     


