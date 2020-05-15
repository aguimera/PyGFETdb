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
from PyGFETdb.DataClass import FnoiseTh 
  

# load dataset from file
FileIn = 'DbData2.pkl'
(dfDat, ExtraVars) = pickle.load(open(FileIn, 'rb'))
(ScalarValues, ArrayValues, ArrayValuesNorm, Vgs, VgsNorm) = ExtraVars

# %% plotting scalar values

# Boxplots
plt.figure()
sns.boxplot(x='Vrms01',
            y='Device',
            hue='Width',
            data=dfDat)

# Values distribution
fig, ax = plt.subplots()
sns.kdeplot(data=np.log10(dfDat.query("Width == 50e-6")['Vrms01']),
            shade=True,                                
            label='W=50um',
            ax=ax)
sns.kdeplot(data=np.log10(dfDat.query("Width == 100e-6")['Vrms01']),
            shade=True,                                
            label='W=100um',
            ax=ax)
ax.set_xlabel('log(Vrms)')

fig, ax = plt.subplots()
sns.kdeplot(data=np.log10(dfDat.query("Wafer == 'B12708W2'")['Vrms01']),
            shade=True,                                
            label='B12752W2',
            ax=ax)
sns.kdeplot(data=np.log10(dfDat.query("Wafer == 'B12752W2'")['Vrms01']),
            shade=True,                                
            label='B12752W2',
            ax=ax)
ax.set_xlabel('log(Vrms)')
