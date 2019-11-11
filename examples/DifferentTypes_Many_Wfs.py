#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: dragc

"""
import os

import matplotlib.pyplot as plt
import quantities as pq
from matplotlib.patches import Patch

import PyGFETdb.DBSearch as DbSe
import PyGFETdb.GlobalFunctions as g
from PyGFETdb import qty, multithrds, Multiprocessing as mp

qtyDebug = False  # True  # Set to false to print less debug info of the rescaling examples

# qty.setActive(False)  # Uncomment to deactivate Quantity Support

plt.close('all')

Colors = ('r', 'g', 'b', 'm', 'y', 'k')

Wafers1 = (
    'B12708W2',  # (in vivo Rob, slices Mavi) Very good
    'B12142W46',  # (in vivo Rob) # High doping
    'B12142W15',  # (in vivo Rob)
    'B11870W8',  # (IDIBAPS implants)
    'B11601W4',
)
Wafers2 = (
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

     }
]

Conditions1 = {'Wafers.Name = ': Wafers1,
               'Devices.Name != ': ('B12708W2-M6',),
               'Devices.Name  != ': ('B12142W15-DUT',),
               }

Conditions2 = {'Wafers.Name = ': Wafers2,
               'Devices.Name != ': ('B12708W2-M6',),
               'Devices.Name  != ': ('B12142W15-DUT',),
               }

CharTable = 'ACcharacts'

GrBase1 = {'Conditions': Conditions1,
           'Table': CharTable,
           'Last': True,
           'DataSelectionConfig': DataSelectionConfig
           }

GrBase2 = {'Conditions': Conditions2,
           'Table': CharTable,
           'Last': True,
           'DataSelectionConfig': DataSelectionConfig
           }

Devices = DbSe.FindCommonValues(Parameter='Devices.Name',
                                Conditions=Conditions1)

GrWs = DbSe.GenGroups(GrBase1, 'Wafers.Name', LongName=False)
GrDevs = DbSe.GenGroups(GrBase1, 'Devices.Name', LongName=False)
GrTypes = DbSe.GenGroups(GrBase1, 'TrtTypes.Name', LongName=False)

legendtitle = 'Wafers'
xlabel = 'Types'
handles = list((Patch(color=Colors[i], label=Wafers2[i]) for i in range(0, len(Wafers2))))

arguments = {
    'arg0': {
        'Param': 'Vrms',
        'Vgs': -0.1,
        'Ud0Norm': True,
        'yscale': 'log',
        'NFmin': 10,
        'NFmax': 1000,
        'Units': 'uV',
        'title': 'Vrms',
        'legendTitle': legendtitle,
        'xlabel': xlabel,
        'handles': handles,
    },
    'arg1': {
        'Param': 'Rds',
        'Vgs': 0,
        'Ud0Norm': True,
        'yscale': 'log',
        'Units': 'uV/A',
        'title': 'Rds',
        'legendTitle': legendtitle,
        'xlabel': xlabel,
        'handles': handles,
    },
    'arg2': {
        'Param': 'GMV',
        'Vgs': -0.1,
        'Ud0Norm': True,
        'yscale': 'log',
        'Units': 'uS/V',
        'title': 'GMV',
        'legendTitle': legendtitle,
        'xlabel': xlabel,
        'handles': handles,
    },
    'arg3': {
        'Param': 'Ud0',
        'Units': 'uV',
        'title': 'Ud0',
        'legendTitle': legendtitle,
        'xlabel': xlabel,
        'handles': handles,
    },
    'arg4': {
        'Param': 'Ids',
        'Vgs=': -0.1,
        'Ud0Norm': True,
        'yscale': 'log',
        'Units': 'uA',
        'title': 'Ids',
        'legendTitle': legendtitle,
        'xlabel': xlabel,
        'handles': handles,
    }
}


if multithrds:
    search = mp.SearchDB_MP
    getparams = mp.GetParams_MP
else:
    search = mp.SearchDB
    getparams = mp.GetParams

ResultsDB = search(GrTypes)
argParams = {'ResultsDB': dict(ResultsDB), 'GrWfs': GrTypes, 'arguments': arguments, 'args': arguments}
ResultsParams = getparams(**argParams)
Vals = g.PlotGroup(ResultsParams, GrTypes, arguments)


args = {'arg1': {
    'Param': 'Vrms',
    'Vgs': -0.1,
    'Ud0Norm': True,
    'Units': 'uV',
}}

# DATABASE SEARCH:
GrWs = DbSe.GenGroups(GrBase2, 'Wafers.Name', LongName=False)
for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
    GrTypes = DbSe.GenGroups(Grwc, 'TrtTypes.Name', LongName=False)
    ResultsDB = search(GrTypes)
    ResultsParams[Grwn] = getparams(ResultsDB, GrTypes, args)

# DATA CLASSIFICATION
ResultsPer_Wf_Type = {}
types = []
devs = {}
resPerType = {}
for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
    ResultsPer_Wf_Type[Grwn] = {}
    for iGr, (Grn, Grc) in enumerate(sorted(ResultsParams.get(Grwn).items())):  # Param 0
        for iType, (TGrn, TGrc) in enumerate(Grc.items()):
            quantities = TGrc  # Param 0
            if qty.isActive():
                Vals = pq.Quantity(quantities)  # np.array also works fine
            else:
                Vals = quantities * 1e6
            if Vals is not None:
                ResultsPer_Wf_Type[Grwn][TGrn] = Vals
                work = resPerType.get(TGrn)
                if work is None:
                    work = [Vals]
                    devs[TGrn] = Vals.size
                else:
                    work.append(Vals)
                    devs[TGrn] += Vals.size

                resPerType[TGrn] = work

                if TGrn not in types:
                    types.append(TGrn)

fig, ax = plt.subplots()
# PLOT 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xLab = []
xPos = []
pos = 0
for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
    for nt, typename in enumerate(types):
        Col = Colors[iWf]
        vals = ResultsPer_Wf_Type[Grwn].get(typename)
        if vals is not None:
            xPos.append(pos)
            xLab.append(typename)
            g._BoxplotValsGroup(ax, Col, pos, vals.transpose())
            pos += 1
plt.xticks(xPos, xLab, rotation=45, fontsize='small')
ax.set_ylabel('[uVrms]', fontsize='large')
ax.set_xlabel(xlabel, fontsize='large')
ax.set_title('Noise (10Hz-1kHz)', fontsize='large')
g.Legend(ax, legendtitle, handles)

# PLOT 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig, ax2 = plt.subplots()
xLab = []
xPos = []
for nt, typename in enumerate(types):
    Col = Colors[nt]
    xLab.append(typename)
    xPos.append(nt)
    plot = []
    for item in resPerType[typename]:
        plot.append((item.shape[1] / devs[typename]) * 100)
    g._BoxplotValsGroup(ax2, Col, nt, plot)
plt.xticks(xPos, xLab, rotation=45, fontsize='small')
ax2.set_ylabel('Yield [%]', fontsize='large')

ax2.set_xlabel(xlabel, fontsize='large')
ax2.set_title('Working gSGFETs', fontsize='large')
g.Legend(ax2, legendtitle, handles)

# """
#
plt.show()
os.system("read")
#
#
#
