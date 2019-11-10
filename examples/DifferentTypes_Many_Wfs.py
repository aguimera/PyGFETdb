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

Wafers1 = (
    'B12708W2',  # (in vivo Rob, slices Mavi) Very good
    # 'B12142W46',  # (in vivo Rob) # High doping
    # 'B12142W15',  # (in vivo Rob)
    # 'B11870W8',  # (IDIBAPS implants)
    # 'B11601W4',
)
Wafers2 = (
    'B12708W2',  # (in vivo Rob, slices Mavi) Very good
    'B12142W46',  # (in vivo Rob) # High doping
    'B12142W15',  # (in vivo Rob)
    'B11870W8',  # (IDIBAPS implants)
    'B11601W4',
)
""" 
DataSelectionConfig = [
    {'Param': 'Ud0',  # Parameter to evaluate
     'Range': (200, 500),  # Range of allowed values, (Min, Max)
     'Name': 'UD0y',
     'ParArgs': {'Units': 'mV'}
     },

    #                            {'Param': 'Rds', # Parameter to evaluate
    #                              'Range': (1e3, 11e3), # Range of allowed values, (Min, Max)
    #                              'Function': np.max, # Funtion to apply into given results
    #                              'InSide': True,  # Inside or outside the range, Default True
    #                              'ParArgs': {'Vgs': None, # Bias point to evaluate
    #                                          'Vds': None,
    #                                          'Ud0Norm': False},
    #                              'Name': 'RDSy'}, # OPTIONAL name to the selection group

    {'Param': 'GMV',  # Parameter to evaluate
     'Range': (100, 10000),  # Range of allowed values, (Min, Max)
     'ParArgs': {'Vgs': -0.1,  # Bias point to evaluate
                 'Vds': None,
                 'Ud0Norm': True,
                 'Units': "uS/V"
                 },
     'Name': 'GMVy'},

    {'Param': 'Vrms',  # Parameter to evaluate
     'ParArgs': {
         'Vgs': -0.1,  # Bias point to evaluate
         'Vds': None,
         'Ud0Norm': True,
         'NFmin': 10,
         'NFmax': 1000,
         'Units': 'V'
     },
     'Range': (5e-6, 0.6e-4), '# Range of allowed values, (Min, Max)
     #                              'Function': np.min, # Funtion to apply into given results
     },
]
"""

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

arguments = {
    'arg0': {
        'Param': 'Vrms',
        'Vgs': -0.1,
        'Ud0Norm': True,
        'yscale': 'log',
        'NFmin': 10,
        'NFmax': 1000,
        'Units': 'uV',
        'title': 'Vrms of Waffer {}'.format(Wafers1[0]),
    },
    'arg1': {
        'Param': 'Rds',
        'Vgs': 0,
        'Ud0Norm': True,
        'yscale': 'log',
        'Units': 'uV/A',
        'title': 'Rds of Waffer {}'.format(Wafers1[0]),
    },
    'arg2': {
        'Param': 'GMV',
        'Vgs': -0.1,
        'Ud0Norm': True,
        'yscale': 'log',
        'Units': 'uS/V',
        'title': 'GMV of Waffer {}'.format(Wafers1[0]),
    },
    'arg3': {
        'Param': 'Ud0',
        'Units': 'uV',
        'title': 'Ud0 of Waffer {}'.format(Wafers1[0]),

    },
    'arg4': {
        'Param': 'Ids',
        'Vgs=': -0.1,
        'Ud0Norm': True,
        'yscale': 'log',
        'Units': 'uA',
        'title': 'Ids of Waffer {}'.format(Wafers1[0]),
    }
}

if multithrds:
    search = mp.SearchDB_MP
    getparams = mp.GetParams_MP
else:
    search = mp.SearchDB
    getparams = mp.GetParams
titlelist = []

ResultsDB = search(GrTypes)
argParams = {'ResultsDB': dict(ResultsDB), 'GrWfs': GrTypes, 'arguments': arguments, 'args': arguments}
ResultsParams = getparams(**argParams)
Vals = g.PlotGroup(ResultsParams, GrTypes, arguments)

fig, ax = plt.subplots()

Colors = ('r', 'g', 'b', 'm', 'y', 'k')

args = {'arg1': {
    'Param': 'Vrms',
    'Vgs': -0.1,
    'Ud0Norm': True,
    'Units': 'uV',
}}
GrWs = DbSe.GenGroups(GrBase2, 'Wafers.Name', LongName=False)
Results = {}
ResultsParams = {}
pos = 0
xLab = []
xPos = []
types = []
res = {}
work = None
devs = {}
for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
    GrTypes = DbSe.GenGroups(Grwc, 'TrtTypes.Name', LongName=False)
    Results[Grwn] = {}
    ResultsDB = search(GrTypes)
    ResultsParams[iWf] = getparams(ResultsDB, GrTypes, args)
    for iType, (Grn, Grc) in enumerate(sorted(ResultsParams[iWf]['arg1'].items())):  # Param 0
        if len(Grc):
            quantities = Grc  # Param 0
            if qty.isActive():
                Vals = pq.Quantity(quantities)  # np.array also works fine
            else:
                Vals = quantities * 1e6
            Results[Grwn][Grn] = Vals
            work = res.get(Grn)
            if work is None:
                work = [Vals]
                devs[Grn] = Vals.size
            else:
                work.append(Vals)
                devs[Grn] += Vals.size
            res[Grn] = work
            xPos.append(pos)
            xLab.append(Grn)
            types.append(Grn)
            g._BoxplotValsGroup(ax, Colors[iWf], pos, Vals.transpose())
            pos += 1

plt.xticks(xPos, xLab, rotation=45, fontsize='small')
ax.set_ylabel('[uVrms]', fontsize='large')
ax.set_xlabel('Types', fontsize='large')
ax.set_title('Noise (10Hz-1kHz) ( {} x Wafers)'.format(len(Wafers2)), fontsize='large')

# %%

fig, ax2 = plt.subplots()
xLab = []
xPos = []
pos = 0
for iGrn, (Grn, Grc) in enumerate((sorted(res.items()))):
    work = []
    for vals in Grc:
        if len(Grc):
            work.append((vals.shape[1] / devs[Grn]) * 100)
    xPos.append(pos)
    xLab.append(Grn)
    g._BoxplotValsGroup(ax2, Colors[iGrn], pos, work)
    pos += 1
plt.xticks(xPos, xLab, rotation=45, fontsize='small')
ax2.set_ylabel('Yield [%]', fontsize='large')
ax2.set_xlabel('Types', fontsize='large')
ax2.set_title('Working gSGFETs ( {} x Wafers)'.format(len(Wafers2)), fontsize='large')
# """
#
plt.show()
# os.system("read")
#
#
#
