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

# MULTIPROCESSING INITIALIZATION #################################################
if multithrds:
    search = mp.SearchDB_MP
    getparams = mp.GetParams_MP
else:
    search = mp.SearchDB
    getparams = mp.GetParams


def DBSearchPerWaferAndType(GrBase, args):
    GrWs = DbSe.GenGroups(GrBase, 'Wafers.Name', LongName=False)
    ResultsParams = {}
    for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
        GrTypes = DbSe.GenGroups(Grwc, 'TrtTypes.Name', LongName=False)
        ResultsDB = search(GrTypes)
        ResultsParams[Grwn] = getparams(ResultsDB, GrTypes, args)
    return GrWs, ResultsParams


def DBSearchPerType(GrBase, args):
    GrTypes = DbSe.GenGroups(GrBase, 'TrtTypes.Name', LongName=False)
    ResultsParams = {}
    for iType, (nType, cType) in enumerate(GrTypes.items()):
        GrWfs = DbSe.GenGroups(cType, 'Wafers.Name', LongName=False)
        ResultsDB = search(GrWfs)
        ResultsParams[nType] = getparams(ResultsDB, GrWfs, args)
    return GrTypes, ResultsParams


def DataClassification(GrWs, ResultsParams):
    clssfResults = {}
    Results = {}
    for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
        for iarg, (narg, arg) in enumerate(sorted(ResultsParams.get(Grwn).items())):  # Params
            for iType, (TGrn, TGrc) in enumerate(arg.items()):
                quantities = TGrc  # Param 0
                if qty.isActive():
                    Vals = pq.Quantity(quantities)  # np.array also works fine
                else:
                    Vals = quantities * 1e6
                if Vals is not None:
                    g.updateDictOfDicts(Results, Grwn, TGrn, Vals)
            clssfResults[narg] = Results
    return clssfResults


def PlotPerTypeNoise(Results, handles=None, xlabel=None, legendtitle=None, Colors=None, perType="", **kwargs):
    # PLOT 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig, ax = plt.subplots()
    xLab = []
    xPos = []
    pos = 0
    for iWf, (Grwn, Grwc) in enumerate(Results.items()):
        Col = Colors[iWf]
        for nt, (typename, vals) in enumerate(Grwc.items()):
            if vals is not None:
                xPos.append(pos)
                xLab.append(typename)
                g._BoxplotValsGroup(ax, Col, pos, vals.transpose())
                pos += 1
    g._closeBoxplotValsGroup(ax, xPos, xLab, xlabel, "[uVrms]", "Noise (10Hz-1kHz) " + perType, **kwargs)
    g.Legend(ax, legendtitle, handles)


def PlotPerTypeYield(Results, title=None, handles=None, xlabel=None, Colors=None, legendtitle=None, **kwargs):
    # PLOTS 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fig, ax2 = plt.subplots()
    xLab = []
    xPos = []
    pos = 0

    for iWf, (wn, dd) in enumerate(Results.items()):
        Col = Colors[iWf]
        temp = []
        for iType, (Grn, Grc) in enumerate(sorted(dd.items())):  # Param 0
            temp.append(Grc.shape[1])
        totalWf = np.sum(temp)

        work = []
        for iType, (Grn, Grc) in enumerate(sorted(dd.items())):  # Param 0
            work.append((Grc.shape[1] / totalWf)*100)
            xLab.append(Grn)
            xPos.append(pos)
            g._BoxplotValsGroup(ax2, Col, pos, work)
            pos += 1
    g._closeBoxplotValsGroup(ax2, xPos, xLab, xlabel, "Yield [%]", title, **kwargs)
    g.Legend(ax2, legendtitle, handles)


def PlotPerTypeYieldTotal(Results, title=None, Colors=None, xlabel=None,
                          **kwargs):
    # PLOT 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig, ax2 = plt.subplots()
    xLab = []
    xPos = []

    temp = []
    for nt, (typename, vtype) in enumerate(Results.items()):
        for iWf, (nWf, cWf) in enumerate(vtype.items()):
            temp.append(cWf.shape[1])
    total = np.sum(temp)

    work = []
    for nt, (typename, vtype) in enumerate(Results.items()):
        for iWf, (nWf, cWf) in enumerate(vtype.items()):
            Col = Colors[nt]
            xLab.append(typename)
            xPos.append(nt)
            work.append((cWf.shape[1] / total) * 100)
        g._BoxplotValsGroup(ax2, Col, nt, work, **kwargs)
    g._closeBoxplotValsGroup(ax2, xPos, xLab, xlabel, "Yield [%]", title, **kwargs)


#############################
# PLOTS PER WAFER AND TYPE
############################
def PlotsPerWaferAnType(GrBase, **arguments):
    # Devices = DbSe.FindCommonValues(Parameter='Devices.Name', Conditions=Conditions)
    # GrWs = DbSe.GenGroups(GrBase, 'Wafers.Name', LongName=False)
    # GrDevs = DbSe.GenGroups(GrBase, 'Devices.Name', LongName=False)
    GrTypes = DbSe.GenGroups(GrBase, 'TrtTypes.Name', LongName=False)
    ResultsDB = search(GrTypes)
    argParams = {'ResultsDB': dict(ResultsDB), 'GrWfs': GrTypes, 'arguments': arguments, 'args': arguments}
    ResultsParams = getparams(**argParams)
    Vals = g.PlotGroup(ResultsParams, GrTypes, arguments)
    return Vals


# END PLOTS PER WAFER ################################################################
# """

# """ # Add '#' at the beginning of this line to active
############################
# PLOTS PER 1 WAFER AND TYPES
###########################
def PlotsPer1WaferAndTypes(GrBase, **kwargs):
    # DATABASE SEARCH ####################################################################################
    GrWs, ResultsParams = DBSearchPerWaferAndType(GrBase, kwargs)
    # DATA CLASSIFICATION ################################################################################
    Results = DataClassification(GrWs, ResultsParams)

    #######
    arg = kwargs['arg1']
    data = Results['arg1']
    # PLOT 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PlotPerTypeNoise(data, perType="x Wafer", **arg)
    # PLOT 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PlotPerTypeYield(data, title="Working SGFETs x Wafer", **arg)
    # PLOT 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PlotPerTypeYieldTotal(data, title="Working SGFETs x Type", **arg)


def PlotsPerTypes(GrBase, Colors=None, **kwargs):
    # ARGUMENT SPECIFICATION #############################################################################

    # DATABASE SEARCH ####################################################################################
    GrTypes, ResultsParams = DBSearchPerType(GrBase, kwargs)

    types = []
    for item in GrTypes.items():
        types.append(item)

    handles = list((Patch(color=Colors[i], label=sorted(list(GrTypes.keys()))[i])
                    for i in range(0, len(list(GrTypes.keys())))))

    # DATA CLASSIFICATION ################################################################################
    Results = DataClassification(GrTypes, ResultsParams)

    #######
    arg = kwargs['arg1']
    data = Results['arg1']
    # PLOT 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PlotPerTypeNoise(data, handles=handles, perType="x Type", **arg)
    # PLOT 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PlotPerTypeYield(data, title="Working SGFETs x Wafer", handles=handles, **arg)
    # PLOT 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PlotPerTypeYieldTotal(data, title="Working SGFETs x Type", handles=handles, **arg)


def main():
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
        # 'B11870W8',  # (IDIBAPS implants)
        # 'B11601W4',
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

    # PLOT GLOBALS ####################################################################
    Colors = ('r', 'g', 'b', 'm', 'y', 'k')
    legendtitle = 'Wafers'
    xlabel = 'Types'
    handles1 = list((Patch(color=Colors[i], label=sorted(Wafers1)[i]) for i in range(0, len(Wafers1))))
    handles2 = list((Patch(color=Colors[i], label=sorted(Wafers2)[i]) for i in range(0, len(Wafers2))))
    arguments1 = {
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

    arguments2 = {
        'arg1': {
            'Param': 'Vrms',
            'Vgs': -0.1,
            'Ud0Norm': True,
            'Units': 'uV',
            'legendTitle': legendtitle,
            'handles': handles2,
            'Colors': Colors
        }
    }

    arguments3 = {
        'arg1': {
            'Param': 'Vrms',
            'Vgs': -0.1,
            'Ud0Norm': True,
            'Units': 'uV',
            'legendTitle': 'Types',
            'Colors': Colors
        },
        'Colors': Colors
    }

    PlotsPerWaferAnType(GrBase1, **arguments1)
    PlotsPer1WaferAndTypes(GrBase2, **arguments2)
    PlotsPerTypes(GrBase3, **arguments3)


# """"""""""""""""""""""""""""""""""""""""""""""
#
main()
plt.show()
#os.system("read")
#
#
#
