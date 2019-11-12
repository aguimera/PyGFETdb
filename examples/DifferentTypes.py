#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: dragc

"""

import matplotlib.pyplot as plt
import numpy as np
import quantities as pq
from matplotlib.patches import Patch

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
    # 'B12142W46',  # (in vivo Rob) # High doping
    # 'B12142W15',  # (in vivo Rob)
    #'B11870W8',  # (IDIBAPS implants)
    #'B11601W4',
)

Wafers3 = (
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

Conditions3 = {'Wafers.Name = ': Wafers3,
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

GrBase3 = {'Conditions': Conditions3,
           'Table': CharTable,
           'Last': True,
           'DataSelectionConfig': DataSelectionConfig
           }

Devices = DbSe.FindCommonValues(Parameter='Devices.Name',
                                Conditions=Conditions1)

GrWs = DbSe.GenGroups(GrBase1, 'Wafers.Name', LongName=False)
GrDevs = DbSe.GenGroups(GrBase1, 'Devices.Name', LongName=False)
GrTypes = DbSe.GenGroups(GrBase1, 'TrtTypes.Name', LongName=False)

# PLOT GLOBALS ####################################################################
Colors = ('r', 'g', 'b', 'm', 'y', 'k')
legendtitle = 'Wafers'
xlabel = 'Types'
handles1 = list((Patch(color=Colors[i], label=sorted(Wafers1)[i]) for i in range(0, len(Wafers1))))
handles2 = list((Patch(color=Colors[i], label=sorted(Wafers2)[i]) for i in range(0, len(Wafers2))))
handles3 = list((Patch(color=Colors[i], label=sorted(Wafers3)[i]) for i in range(0, len(Wafers3))))

# MULTIPROCESSING INITIALIZATION #################################################
if multithrds:
    search = mp.SearchDB_MP
    getparams = mp.GetParams_MP
else:
    search = mp.SearchDB
    getparams = mp.GetParams

# """ # Add '#' at the begging of this line to active
#############################
# PLOTS PER WAFER AND TYPE
############################
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
        'handles': handles1,
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
        'handles': handles1,
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
        'handles': handles1,
    },
    'arg3': {
        'Param': 'Ud0',
        'Units': 'uV',
        'title': 'Ud0',
        'legendTitle': legendtitle,
        'xlabel': xlabel,
        'handles': handles1,
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
        'handles': handles1,
    }
}

ResultsDB = search(GrTypes)
argParams = {'ResultsDB': dict(ResultsDB), 'GrWfs': GrTypes, 'arguments': arguments, 'args': arguments}
ResultsParams = getparams(**argParams)
Vals = g.PlotGroup(ResultsParams, GrTypes, arguments)


# END PLOTS PER WAFFER ################################################################
# """


# """ # Add '#' at the begging of this line to active
############################
# PLOTS PER TYPE AND WAFFER
###########################
def DBSearchPerType(GrBase, args):
    GrWs = DbSe.GenGroups(GrBase, 'Wafers.Name', LongName=False)
    ResultsParams = {}
    for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
        GrTypes = DbSe.GenGroups(Grwc, 'TrtTypes.Name', LongName=False)
        ResultsDB = search(GrTypes)
        ResultsParams[Grwn] = getparams(ResultsDB, GrTypes, args)
    return GrWs, ResultsParams


def DataClassification(GrWs, ResultsParams):
    clssfResults = {}
    ResultsPer_Wf_Type = {}
    types = []
    resPerType = {}
    for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
        for iarg, (narg, arg) in enumerate(sorted(ResultsParams.get(Grwn).items())):  # Params
            for iType, (TGrn, TGrc) in enumerate(arg.items()):
                if TGrn not in types:
                    types.append(TGrn)
                quantities = TGrc  # Param 0
                if qty.isActive():
                    Vals = pq.Quantity(quantities)  # np.array also works fine
                else:
                    Vals = quantities * 1e6
                if Vals is not None:
                    g.updateDictOfLists(resPerType, TGrn, Vals)
                    g.updateDictOfDicts(ResultsPer_Wf_Type, Grwn, TGrn, Vals)
            clssfResults[narg] = {
                'ResultsPer_Wf_Type': ResultsPer_Wf_Type,
                'types': types,
                'resPerType': resPerType,
            }
    return clssfResults


def PlotPerTypeNoise(types=None, ResultsPer_Wf_Type=None, handles=None, **kwargs):
    # PLOT 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig, ax = plt.subplots()
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
    g._closeBoxplotValsGroup(ax, xPos, xLab, xlabel, "[uVrms]", "Noise (10Hz-1kHz)", **kwargs)
    g.Legend(ax, legendtitle, handles)


def PlotPerTypeYield(ResultsPer_Wf_Type=None, handles=None, **kwargs):
    # PLOTS 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp = []
    for iWf, (wn, dd) in enumerate(ResultsPer_Wf_Type.items()):
        for iType, (Grn, Grc) in enumerate(sorted(dd.items())):  # Param 0
            temp.append(Grc.shape[1])
    total = np.sum(temp)

    fig, ax2 = plt.subplots()
    xLab = []
    xPos = []
    pos = 0

    for iWf, (wn, dd) in enumerate(ResultsPer_Wf_Type.items()):
        temp = []
        for iType, (Grn, Grc) in enumerate(sorted(dd.items())):  # Param 0
            temp.append(Grc.shape[1])
        totalWf = np.sum(temp)
        work = []
        for iType, (Grn, Grc) in enumerate(sorted(dd.items())):  # Param 0
            work.append((Grc.shape[1] / totalWf)*100)
            xLab.append(Grn)
            xPos.append(pos)
            Col = Colors[iWf]
            g._BoxplotValsGroup(ax2, Col, pos, work)
            pos += 1
    g._closeBoxplotValsGroup(ax2, xPos, xLab, xlabel, "Yield [%]", "Working SGFETs x Waffer", **kwargs)
    g.Legend(ax2, legendtitle, handles)


def PlotPerTypeYieldMultiWaffer(types=None, resPerType=None,
                                **kwargs):
    # PLOT 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig, ax2 = plt.subplots()
    xLab = []
    xPos = []
    for nt, typename in enumerate(types):
        Col = Colors[nt]
        xLab.append(typename)
        xPos.append(nt)
        work = []
        temp = []
        for item in resPerType[typename]:
            temp.append(item.shape[1])
        total = np.sum(temp)
        for item in resPerType[typename]:
            work.append((item.shape[1] / total) * 100)
        g._BoxplotValsGroup(ax2, Col, nt, work, **kwargs)
    g._closeBoxplotValsGroup(ax2, xPos, xLab, xlabel, "Yield [%]", "Working SGFETs x Type", **kwargs)


# ARGUMENT SPECIFICATION #############################################################################
args = {'arg1': {
    'Param': 'Vrms',
    'Vgs': -0.1,
    'Ud0Norm': True,
    'Units': 'uV',
    'legendTitle': legendtitle,
    'handles': handles2,
}}

# DATABASE SEARCH ####################################################################################
GrWs, ResultsParams = DBSearchPerType(GrBase2, args)
# DATA CLASSIFICATION ################################################################################
DataClasf = DataClassification(GrWs, ResultsParams)

#######
kwargs = args['arg1']
data = DataClasf['arg1']
# PLOT 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotPerTypeNoise(**data, **kwargs)
# PLOT 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotPerTypeYield(**data, **kwargs)

################################################################
# PLOT 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ARGUMENT SPECIFICATION #############################################################################
args = {'arg1': {
    'Param': 'Vrms',
    'Vgs': -0.1,
    'Ud0Norm': True,
    'Units': 'uV',
    'legendTitle': legendtitle,
    'handles': handles3,
}}

# DATABASE SEARCH ####################################################################################
GrWs, ResultsParams = DBSearchPerType(GrBase3, args)
# DATA CLASSIFICATION ################################################################################
DataClasf = DataClassification(GrWs, ResultsParams)

#######
kwargs = args['arg1']
data = DataClasf['arg1']
PlotPerTypeYieldMultiWaffer(**data, **kwargs)

# END PLOTS PER TYPE AND WAFFER #################################################################################
#"""
#
plt.show()
#os.system("read")
#
#
#
