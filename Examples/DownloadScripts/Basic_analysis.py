#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 14:41:57 2022

@author: aguimera
"""
import pandas as pd
import seaborn as sns
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import itertools
import numpy as np
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages


# %% load data

dfDat = pd.read_pickle('Data.pkl')
dfAttrs = dfDat.attrs

print('Metadata found', dfAttrs.keys())

# %% Quick data filtering

dSel = dfDat.copy()

queries = ("RdsCNP > 2 ",
           "IsOk == True",
           # "CNP < 350 ",
           # "CNP > 100 ",
           # "Irms01 < 6",
            )

for query in queries:
    dSel.query(query, inplace=True)

# %% Boxplots representations
"""
Boxplot representation of several scalar values

Objective manage autmatically the labels and units

can work with boxplot, swarmplot .... and similars

"""

# columns to plot
PlotPars = {
            'CNP': {'Column': 'CNP', 'Ylog': False},
            'Ids': {'Column': 'Ids01', 'Ylog': False},
            'GM': {'Column': 'GM01', 'Ylog': False},
            'Rds': {'Column': 'RdsCNP', 'Ylog': False},
            'Irms': {'Column': 'Irms01', 'Ylog': True},
            'Vrms': {'Column': 'Vrms01', 'Ylog': True},
            }

# figure gneration
nRows = 2
nCols = math.ceil(len(PlotPars)/nRows)
fig, Axs = plt.subplots(nrows=nRows, ncols=nCols)
Axs = Axs.flatten()

# plotting
AxsDict = {}
for ip, (parn, dic) in enumerate(PlotPars.items()):
    ax = Axs[ip]
    Column = dic['Column']

    bg = sns.boxplot(x='Device',
                     y=Column,
                     hue='Date',
                     data=dSel,
                     ax=Axs[ip])

    # Format axis
    AxsDict[parn] = bg
    # bg.set_title(parn)
    slabel = '{}[{}]'.format(parn, dfAttrs['ColUnits'][Column])
    bg.set_ylabel(slabel)
    if dic['Ylog']:
        ax.set_yscale('log')
    bg.set_xticklabels(bg.get_xticklabels(), rotation=20, fontsize='small')
    bg.legend(loc='upper left', fontsize='xx-small')

    
# %% Vectors plotting
"""
Vector representation of several vector values

Objective manage autmatically the labels and units

"""

# Colors iteration
Colors = itertools.cycle(mpl.colors.TABLEAU_COLORS)

# Grouping rows for unique values of some column
# dgroups = df2Plot.groupby('Device')
dgroups = dSel.groupby('Date')

# Columns to plot
PlotPars = {
            'Ids': {'Column': 'Ids', 'Ylog': False},
            'GM': {'Column': 'GMNorm', 'Ylog': False},
            'Vrms': {'Column': 'VrmsNorm', 'Ylog': True},
            'Irms': {'Column': 'IrmsNorm', 'Ylog': True},
            'NoA': {'Column': 'NoANorm', 'Ylog': True},
            'NoB': {'Column': 'NoBNorm', 'Ylog': False},
            'NoC': {'Column': 'NoCNorm', 'Ylog': True},
            }

# figure generation
nRows = 2
nCols = math.ceil(len(PlotPars)/nRows)
fig, Axs = plt.subplots(nrows=nRows, ncols=nCols)
Axs = Axs.flatten()

# Format axis
AxsDict = {}
for ip, (parn, dic) in enumerate(PlotPars.items()):
    ax = Axs[ip]
    AxsDict[parn] = ax
    par = dic['Column']

    slabel = '{}[{}]'.format(parn, dfAttrs['ColUnits'][Column])
    ax.set_ylabel(slabel)

    if dic['Ylog']:
        ax.set_yscale('log')

    # select Vgs vector
    if par.endswith('Norm'):
        dic['xVar'] = dfAttrs['VgsNorm']
        ax.set_xlabel('Vgs - CNP [V]')
    else:
        ax.set_xlabel('Vgs [V]')
        dic['xVar'] = dfAttrs['Vgs']

# Plotting
for gn in dgroups.groups:
    gg = dgroups.get_group(gn)
    Col = next(Colors)

    for parn, dic in PlotPars.items():
        Vals = np.array([])
        for index, row in gg.iterrows():
            v = row[dic['Column']]
            if type(v) == float:
                continue
            try:
                Vals = np.vstack((Vals, v)) if Vals.size else v
            except:
                pass

        xVar = dic['xVar']
        ax = AxsDict[parn]

        Vals = Vals.magnitude.transpose()
        mean = np.nanmedian(Vals, axis=1)
        std = stats.tstd(Vals, axis=1) # changed to mad for robustness

        ax.plot(xVar, Vals, color=Col, alpha=0.1)
        ax.plot(xVar, mean, color=Col, lw=1.5, label=gn)
        if not dic['Ylog']:
            ax.fill_between(xVar, mean+std, mean-std, color=Col, alpha=0.2)

# %% Arrays plotting, advanced
"""
Plotting of raw PSD or Bode
As it will generate a lot of figure is better to generate a PDF output
"""

PDF = PdfPages('BodeGraph.pdf')

for index, row in dfDat.iterrows():
    fig, (axM, axP) = plt.subplots(2,1)
    char = row['CharCl']    
    axM.loglog(char.GetFgm(), char.GetGmMag()/2)
    axM.plot(char.GetFgm()[0], np.abs(char.GetGM()), 'k*')
    axP.semilogx(char.GetFgm(), char.GetGmPh())
    PDF.savefig(fig)
    plt.close('all')
    
PDF.close()




