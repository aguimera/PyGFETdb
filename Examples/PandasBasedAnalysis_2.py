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
import matplotlib as mpl
from matplotlib.artist import ArtistInspector
import matplotlib.colors as mcolors
import itertools
from PyGFETdb.DataClass import FnoiseTh 
  

# load dataset from file
FileIn = 'DbData2.pkl'
(dfDat, ExtraVars) = pickle.load(open(FileIn, 'rb'))
(ScalarValues, ArrayValues, ArrayValuesNorm, Vgs, VgsNorm) = ExtraVars

# %% plotting array values
df2Plot = dfDat.copy()
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=mpl.colors.CSS4_COLORS.values())

# data selection
df2Plot.query("IsOk == True", inplace=True)

# Plot Ids 
plt.figure()
Vals = np.array([])
for index, row in df2Plot.iterrows():
    v = row['Ids']
    print(Vals.shape, v.shape)
    Vals = np.vstack((Vals, v)) if Vals.size else v
Vals = Vals.transpose() 

plt.plot(Vgs, Vals)
plt.xlabel('Vgs [V]')
plt.ylabel('Ids [uA]')

# Plot ids normalized
plt.figure()
Vals = np.array([])
for index, row in df2Plot.iterrows():
    v = row['IdsNorm']
    Vals = np.vstack((Vals, v)) if Vals.size else v
Vals = Vals.transpose()
   
plt.plot(VgsNorm, Vals)
plt.xlabel('Vgs - CNP [V]')
plt.ylabel('Ids [uA]')
    

# %% Grouping and calculate means

df2Plot = dfDat.copy()

df2Plot.query("IsOk == True", inplace=True)

print(df2Plot.Device.value_counts())

# dgroups = df2Plot.groupby('Device')
dgroups = df2Plot.groupby('Wafer')

# Plot Ids 
for gn in dgroups.groups:
    gg = dgroups.get_group(gn)

    plt.figure()
    Vals = np.array([])
    for index, row in gg.iterrows():
        v = row['IdsNorm']    
        Vals = np.vstack((Vals, v)) if Vals.size else v
    Vals = Vals.magnitude.transpose() 
    
    plt.plot(VgsNorm, Vals, alpha=0.1)
    mean = np.nanmean(Vals, axis=1)
    std = np.nanstd(Vals, axis=1)
    plt.plot(VgsNorm, mean, lw=1.5)
    plt.fill_between(VgsNorm, mean+std, mean-std, alpha=0.4)
    plt.title(gn)
    

# %% debuggin noise fitting

df2Plot = dfDat.copy()

df2Plot.query("IsOk == True", inplace=True)
df2Plot.query("Device == 'B12708W2-S9' ", inplace=True)

Pars = ('NoANorm', 'NoBNorm', 'NoCNorm', 'IrmsNorm', 'VrmsNorm')

for par in Pars:
    plt.figure()
    Vals = np.array([])
    for index, row in df2Plot.iterrows():
        v = row[par]
        print(Vals.shape, v.shape)
        Vals = np.vstack((Vals, v)) if Vals.size else v
    Vals = Vals.magnitude.transpose() 
    
    plt.plot(VgsNorm, Vals, alpha=0.1)
    mean = np.nanmean(Vals, axis=1)
    std = np.nanstd(Vals, axis=1)
    plt.semilogy(VgsNorm, mean, lw=1.5)
    plt.fill_between(VgsNorm, mean+std, mean-std, alpha=0.4)


# %%
df2Plot = dfDat.copy()
df2Plot.query("IsOk == True", inplace=True)
df2Plot.query("Device == 'B12708W2-S9' ", inplace=True)


dgroups = df2Plot.groupby('TrtName')

for gn in dgroups.groups:
    gg = dgroups.get_group(gn)

    plt.figure()    
    CharCl = gg.CharCl.values[0]
    fpsd = CharCl.GetFpsd()
    PSD = CharCl.GetPSD()
    NoA = CharCl.GetNoA().flatten()
    NoB = CharCl.GetNoB().flatten()
    NoC = CharCl.GetNoC().flatten()   
    for p, a, b, c in zip(PSD.transpose(), NoA, NoB, NoC):
        plt.loglog(fpsd, p)
        fp = FnoiseTh(fpsd, a, b, c)
        plt.loglog(fpsd, fp, '--k')
        