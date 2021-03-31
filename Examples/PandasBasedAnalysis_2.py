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
# Create a copy of the data frame to aviod destroy values in the selections
df2Plot = dfDat.copy()
# Increase the default color lines in matplotlib
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=mpl.colors.CSS4_COLORS.values())

# data selection
df2Plot.query("IsOk == True", inplace=True)

# Plot Ids 
plt.figure()
Vals = np.array([])
for index, row in df2Plot.iterrows():
    v = row['Ids']
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

# Print unique values for a column
print(df2Plot.Device.value_counts())
# Create a colors cycle
Colors = itertools.cycle(mpl.colors.TABLEAU_COLORS)

# Grouping rows for unique values of some column
# dgroups = df2Plot.groupby('Device')
dgroups = df2Plot.groupby('Wafer')

plt.figure()
# Plot Ids for groups
for gn in dgroups.groups:
    gg = dgroups.get_group(gn)

    Vals = np.array([])
    for index, row in gg.iterrows():
        v = row['IdsNorm']
        Vals = np.vstack((Vals, v)) if Vals.size else v
    Vals = Vals.magnitude.transpose()

    Col = next(Colors)
    plt.plot(VgsNorm, Vals, color=Col, alpha=0.2)
    mean = np.nanmean(Vals, axis=1)
    std = np.nanstd(Vals, axis=1)
    plt.plot(VgsNorm, mean, color=Col, lw=1.5, label=gn)
    plt.fill_between(VgsNorm, mean+std, mean-std, color=Col, alpha=0.4)
plt.legend()
plt.xlabel('Vgs - CNP [V]')
plt.ylabel('Ids [uA]')

# %% Noise Parameters

df2Plot = dfDat.copy()


df2Plot.query("IsOk == True", inplace=True)
df2Plot.query("Device == 'B12708W2-S9' ", inplace=True)

# Pars = ('NoANorm', 'NoBNorm', 'NoCNorm', 'IrmsNorm', 'VrmsNorm')
# xvar = VgsNorm
# ArVals = ArrayValuesNorm
Pars = ('NoA', 'NoB', 'NoC', 'Irms', 'Vrms')
xvar = Vgs
ArVals = ArrayValues

nRows = 2
nCols = math.ceil(len(Pars)/2)
fig, Axs = plt.subplots(nrows=nRows, ncols=nCols, sharex=True)
Axs = Axs.flatten()

for ip, par in enumerate(Pars):
    Vals = np.array([])
    for index, row in df2Plot.iterrows():
        v = row[par]
        Vals = np.vstack((Vals, v)) if Vals.size else v
    Vals = Vals.magnitude.transpose()
    Axs[ip].plot(xvar, Vals, 'k', alpha=0.2)
    mean = np.nanmean(Vals, axis=1)
    std = np.nanstd(Vals, axis=1)
    Axs[ip].semilogy(xvar, mean, 'r', lw=1.5)
    Axs[ip].fill_between(xvar, mean+std, mean-std, color='r', alpha=0.2)
    if 'Units' in ArVals[par]:
        Axs[ip].set_ylabel(par + '[' + ArVals[par]['Units'] + ']')
    else:
        Axs[ip].set_ylabel(par)

# %% debug 1/f fitting

df2Plot = dfDat.copy()
df2Plot.query("IsOk == True", inplace=True)
df2Plot.query("Device == 'B12708W2-S9' ", inplace=True)

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=mpl.colors.CSS4_COLORS.values())
mpl.rcParams['axes.grid'] = True

dgroups = df2Plot.groupby('TrtName')

nRows = 4
nCols = math.ceil(len(dgroups)/nRows)
fig, Axs = plt.subplots(nrows=nRows, ncols=nCols, sharex=True, sharey=True)
Axs = Axs.flatten()

for ic, gn in enumerate(dgroups.groups):
    gg = dgroups.get_group(gn)

    CharCl = gg.CharCl.values[0]
    fpsd = CharCl.GetFpsd()
    PSD = CharCl.GetPSD()
    FitPSD = CharCl.GetFitPSD()
    for p, fp in zip(PSD.transpose(), FitPSD.transpose()):
        Axs[ic].loglog(fpsd, p)        
        Axs[ic].loglog(fpsd, fp, '--k')
    Axs[ic].set_title(gn, fontsize='small')

# %% debug 1/f fitting

df2Plot = dfDat.copy()
df2Plot.query("IsOk == True", inplace=True)
df2Plot.query("TrtName=='B12708W2-S9-Ch06' or TrtName=='B12708W2-S9-Ch05'", inplace=True)

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=mpl.colors.CSS4_COLORS.values())
mpl.rcParams['axes.grid'] = True

dgroups = df2Plot.groupby('TrtName')

nRows = 1
nCols = math.ceil(len(dgroups)/nRows)
fig, Axs = plt.subplots(nrows=nRows, ncols=nCols, sharex=True, sharey=True)
Axs = Axs.flatten()

for ic, gn in enumerate(dgroups.groups):
    gg = dgroups.get_group(gn)

    CharCl = gg.CharCl.values[0]
    fpsd = CharCl.GetFpsd()
    PSD = CharCl.GetPSD()
    FitPSD = CharCl.GetFitPSD()
    for p, fp in zip(PSD.transpose(), FitPSD.transpose()):
        Axs[ic].loglog(fpsd, p)        
        Axs[ic].loglog(fpsd, fp, '--k')
    Axs[ic].set_title(gn, fontsize='small')

