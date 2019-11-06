#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: dragc

"""
import os

import matplotlib.pyplot as plt
import numpy as np
import quantities as pq

import PyGFETdb.DBAnalyze as DbAn
import PyGFETdb.DBSearch as DbSe
import PyGFETdb.GlobalFunctions as g
from PyGFETdb import qty, Thread, multithrds

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
     'Range': (5e-6, 0.6e-4),  # Range of allowed values, (Min, Max)
     #                              'Function': np.min, # Funtion to apply into given results
     },
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

args1 = {
    'Param': 'Vrms',
    'Vgs': -0.1,
    'Ud0Norm': True,
    'yscale': 'log',
    'NFmin': 10,
    'NFmax': 1000,
    'Units': 'uV'
}
args2 = {
    'Param': 'Rds',
    'Vgs': 0,
    'Ud0Norm': True,
    'yscale': 'log',
    'Units': 'uV/A'
}
args3 = {
    'Param': 'GMV',
    'Vgs': -0.1,
    'Ud0Norm': True,
    'yscale': 'log',
    'Units': 'uS/V'
}
args4 = {
    'Param': 'Ud0',
    'Units': 'uV'
}
args5 = {
    'Param': 'Ids',
    'Vgs=': -0.1,
    'Ud0Norm': True,
    'yscale': 'log',
    'Units': 'uA'
}

# GetParamsThread
args = [args1, args2, args3, args4, args5]
if multithrds:
    # GetFromDB
    ResultsDB = g.getFromDB(GrDevs)
    pool = Thread.PyFETdb(DbAn)
    g._GetParamsThread(pool, args, ResultsDB)
    ResultsParams = pool.getResults()
    del pool
    Vals = g.Plot(GrWs, ResultsParams, args)
else:
    Vals = {}
    for iarg, arg in enumerate(args):
        xLab = []
        xPos = []
        fig, Ax = plt.subplots()
        qtys = None
        for iGr, (Grn, Grc) in enumerate(sorted(GrWs.items())):
            group = {}
            Data, Trts = DbSe.GetFromDB(**Grc)
            ParamData = DbAn.GetParam(Data, **arg)
            if ParamData is not None:
                g.updateDictOfLists(group, Grn, ParamData)
                Vals.update({iarg: group})
                if qty.isActive():
                    ParamData = qty.flatten(ParamData)
                    qtys = np.array(ParamData)
                ParamData = np.array(ParamData)
                g.PlotValsGroup(Ax, xLab, xPos, iGr, Grn, ParamData, **arg)
        g.closePlotValsGroup(Ax, xLab, xPos, qtys, **arg)

fig, ax = plt.subplots()

Colors = ('r', 'g', 'b', 'm', 'y', 'k')

GrWs = DbSe.GenGroups(GrBase, 'Wafers.Name', LongName=False)
iDev = 0
Results = {}
for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
    GrDs = DbSe.GenGroups(Grwc, 'Devices.Name', LongName=False)
    Col = Colors[iWf]

    Results[Grwn] = {}
    for Grn, Grc in GrDs.items():
        iDev += 1
        Data, t = DbSe.GetFromDB(**Grc)

        quantities = DbAn.GetParam(Data,
                                   Param='Vrms',
                                   Vgs=-0.1,
                                   Ud0Norm=True,
                                   Units='uV',
                                   )

        # startExamples ###################################

        if qty.isActive() and qtyDebug:
            """
                Some examples with Quantities conversion and rescaling
            """

            # test 1
            Vals = qty.rescaleFromKey(quantities, "mV/S")
            print('Vals=', Vals)

            # test 2
            Vals = qty.toQuantity(quantities)

            # We only rescale Quantities with Units
            # Vals.rescale("mV/S") raises exceptions rescaling
            # Quantitites without units, so on those cases
            # we rescale manually
            if Vals.units != pq.dimensionless:
                Vals = Vals.rescale("mV/S")
            else:
                Vals *= 1e3
            print('Vals=', Vals)

        # endExamples ####################################

        """ 
        Finally, we lose the units, 
        but maintained compatibility
        with the rest of the script
        """

        # Final value of Vals
        if qty.isActive():
            Vals = pq.Quantity(quantities)  # np.array also works fine
        else:
            Vals = quantities * 1e6

        Results[Grwn][Grn] = Vals

        bplt = ax.boxplot(Vals.transpose(),
                          positions=(iDev,),
                          patch_artist=True,  # fill with color
                          widths=0.75,
                          sym='+',
                          labels=('',),
                          #                      notch=True,
                          )

        for element in ('boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps'):
            plt.setp(bplt[element], color=Col)

        for fl in bplt['fliers']:
            fl.set_markeredgecolor(Col)

        for patch in bplt['boxes']:
            patch.set(facecolor=Col)
            patch.set(alpha=0.5)

        # ax.set_yscale('log')
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
    Col = Colors[iWf]
    for dn, d in dd.items():
        if len(d):
            work.append((d.shape[1] / 16) * 100)

    xLab.append(wn)
    xPos.append(iWf)
    bplt = ax2.boxplot(work,
                       positions=(iWf,),
                       patch_artist=True,  # fill with color
                       widths=0.75,
                       sym='+')

    for element in ('boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps'):
        plt.setp(bplt[element], color=Col)

    for fl in bplt['fliers']:
        fl.set_markeredgecolor(Col)

    for patch in bplt['boxes']:
        patch.set(facecolor=Col)
        patch.set(alpha=0.5)

plt.xticks(xPos, xLab, rotation=45, fontsize='small')
ax2.set_ylabel('Yield [%]', fontsize='large')
ax2.set_xlabel('Wafers', fontsize='large')
ax2.set_title('Working gSGFETs (16 x Probe)', fontsize='large')
#
#
plt.show()
os.system("read")
#
#
#
