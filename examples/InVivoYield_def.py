#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: dragc

"""

import matplotlib.pyplot as plt
import quantities as pq

import PyGFETdb.DBSearch as DbSe
import PyGFETdb.GlobalFunctions as g
from PyGFETdb import qty, multithrds, Multiprocessing as mp

qtyDebug = False  # True  # Set to false to print less debug info of the rescaling examples

# qty.setActive(False)  # Uncomment to deactivate Quantity Support

plt.close('all')


Wafers = (
    'B12708W2',  # (in vivo Rob, slices Mavi) Very good
    'B12142W46',  # (in vivo Rob) # High doping
    'B12142W15',  # (in vivo Rob)
    'B11870W8',  # (IDIBAPS implants)
    'B11601W4',
)

DataSelectionConfig = [
    {'Param': 'Ud0',  # Parameter to evaluate
     'Range': (0.2, 0.5),  # Range of allowed values, (Min, Max)
     'Name': 'UD0y'},

    #                            {'Param': 'Rds', # Parameter to evaluate
    #                              'Range': (1e3, 11e3), # Range of allowed values, (Min, Max)
    #                              'Function': np.max, # Funtion to apply into given results
    #                              'InSide': True,  # Inside or outside the range, Default True
    #                              'ParArgs': {'Vgs': None, # Bias point to evaluate
    #                                          'Vds': None,
    #                                          'Ud0Norm': False},
    #                              'Name': 'RDSy'}, # OPTIONAL name to the selection group

    {'Param': 'GMV',  # Parameter to evaluate
     'Range': (1e-4, 1e-2),  # Range of allowed values, (Min, Max)
     'ParArgs': {'Vgs': -0.1,  # Bias point to evaluate
                 'Vds': None,
                 'Ud0Norm': True,
                 },
     'Name': 'GMVy'},

    {'Param': 'Vrms',  # Parameter to evaluate
     'ParArgs': {
         'Vgs': -0.1,  # Bias point to evaluate
         'Vds': None,
         'Ud0Norm': True,
         'NFmin': 10,
         'NFmax': 1000,
     },
     'Range': (5e-6, 0.6e-4),  # Range of allowed values, (Min, Max)
     #                              'Function': np.min, # Funtion to apply into given results

     }
]

Conditions = {'Wafers.Name = ': Wafers,
              'Devices.Name != ': ('B12708W2-M6',),
              'Devices.Name  != ': ('B12142W15-DUT',)}
CharTable = 'ACcharacts'

GrBase = {'Conditions': Conditions,
          'Table': CharTable,
          'Last': True,
          'DataSelectionConfig': DataSelectionConfig
          }

Devices = DbSe.FindCommonValues(Parameter='Devices.Name',
                                Conditions=Conditions)

GrWs = DbSe.GenGroups(GrBase, 'Wafers.Name', LongName=False)
GrDevs = DbSe.GenGroups(GrBase, 'Devices.Name', LongName=False)

arguments = {
    'arg0': {
        'Param': 'Vrms',
        'Vgs': -0.1,
        'Ud0Norm': True,
        'yscale': 'log',
        'NFmin': 10,
        'NFmax': 1000,
        'Units': 'uV'
    },
    'arg1': {
        'Param': 'Rds',
        'Vgs': 0,
        'Ud0Norm': True,
        'yscale': 'log',
        'Units': 'uV/A'
    },
    'arg2': {
        'Param': 'GMV',
        'Vgs': -0.1,
        'Ud0Norm': True,
        'yscale': 'log',
        'Units': 'uS/V'
    },
    'arg3': {
        'Param': 'Ud0',
        'Units': 'uV'
    },
    'arg4': {
        'Param': 'Ids',
        'Vgs=': -0.1,
        'Ud0Norm': True,
        'yscale': 'log',
        'Units': 'uA'
    }
}

if multithrds:
    search = mp.SearchDB_MP
    getparams = mp.GetParams_MP
else:
    search = mp.SearchDB
    getparams = mp.GetParams


ResultsDB = search(GrWs)
argParams = {'ResultsDB': dict(ResultsDB), 'GrWfs': GrWs, 'arguments': arguments, 'args': arguments}
ResultsParams = getparams(**argParams)
Vals = g.PlotGroup(ResultsParams, GrWs, arguments)

fig, ax = plt.subplots()

Colors = ('r', 'g', 'b', 'm', 'y', 'k')

args = {'0': {
    'Param': 'Vrms',
    'Vgs': -0.1,
    'Ud0Norm': True,
    'Units': 'uV',
}}
GrWs = DbSe.GenGroups(GrBase, 'Wafers.Name', LongName=False)
Results = {}
for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
    GrDs = DbSe.GenGroups(Grwc, 'Devices.Name', LongName=False)
    Col = Colors[iWf]
    Results[Grwn] = {}
    ResultsDB = search(GrDs)
    ResultsParams = getparams(ResultsDB, GrDs, args)
    for iDev, (Grn, Grc) in enumerate(sorted(ResultsParams['0'].items())):  # Param 0
        quantities = Grc  # Param 0
        if qty.isActive():
            Vals = pq.Quantity(quantities)  # np.array also works fine
        else:
            Vals = quantities * 1e6

        Results[Grwn][Grn] = Vals
    g._BoxplotValsGroup(ax, Col, iWf, Vals.transpose())

ax.set_ylabel('[uVrms]', fontsize='large')
ax.set_xlabel('Probes', fontsize='large')
ax.set_title('Noise (10Hz-1kHz)', fontsize='large')

# %%
fig, ax2 = plt.subplots()

Colors = ('r', 'g', 'b', 'm', 'y', 'k')
xLab = []
xPos = []
for iWf, (wn, dd) in enumerate(Results.items()):
    work = []
    #    (work.shape[0]/16)*100
    Col = Colors[iWf]
    for dn, d in dd.items():
        if len(d):
            work.append((d.shape[1] / 16) * 100)

    xLab.append(wn)
    xPos.append(iWf)
    g._BoxplotValsGroup(ax2, Col, iWf, work)

plt.xticks(xPos, xLab, rotation=45, fontsize='small')
ax2.set_ylabel('Yield [%]', fontsize='large')
ax2.set_xlabel('Wafers', fontsize='large')
ax2.set_title('Working gSGFETs (16 x Probe)', fontsize='large')
#
#
plt.show()
# os.system("read")
#
#
#
