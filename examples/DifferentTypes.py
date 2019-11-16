#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: dragc

"""

import matplotlib.pyplot as plt
import quantities as pq
from matplotlib.patches import Patch

import PyGFETdb.DBSearch as DbSe
import PyGFETdb.GlobalFunctions as g
import PyGFETdb.SearchFunctions as s
from PyGFETdb import PlotFunctions as plot
from PyGFETdb import multithrds, Multiprocessing as mp

# MULTIPROCESSING INITIALIZATION #################################################
if multithrds:
    search = mp.SearchDB_MP
    getparams = mp.GetParams_MP
else:
    search = mp.SearchDB
    getparams = mp.GetParams

#############################
# PLOTS PER WAFER AND TYPE
############################

def PlotsParams(GrBase, arguments, **kwargs):
    GrTypes = DbSe.GenGroups(GrBase, 'TrtTypes.Name', LongName=False)
    ResultsDB = search(GrTypes)
    argParams = {'ResultsDB': dict(ResultsDB), 'GrWfs': GrTypes, 'arguments': arguments, 'args': arguments}
    ResultsParams = getparams(**argParams)
    Vals = plot.PlotGroup(ResultsParams, GrTypes, arguments, **kwargs)
    return Vals


# END PLOTS PER WAFER ################################################################
# """

############################
# PLOTS PER WAFER AND TYPES
###########################
def PlotsPerWaferAndTypes(GrBase, arguments, Colors=None, legendTitle=None, xlabel=None, **kwargs):
    # DATABASE SEARCH ####################################################################################
    GrWs, ResultsParams = s.DBSearchPerWaferAndType(GrBase, arguments, **kwargs)
    # DATA CLASSIFICATION ################################################################################
    Results = g.DataClassification(GrWs, arguments, ResultsParams)
    handles = list((Patch(color=Colors[i], label=sorted(list(GrWs.keys()))[i])
                    for i in range(0, len(list(GrWs.keys())))))

    plot.PlotResults(Results, arguments, Colors=Colors, handles=handles, xlabel=xlabel,
                     legendTitle=legendTitle, **kwargs)

    data = Results['arg5']
    # PLOT 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot.PlotPerTypeNoise(data, legendTitle=legendTitle, handles=handles, Colors=Colors, perType="x Wafer",
                          xlabel=xlabel, **kwargs)
    # PLOT 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot.PlotPerTypeYield(data, Colors=Colors, title="Working SGFETs x Wafer", xlabel=xlabel,
                          legendTitle=legendTitle, perType=" x Wafer", handles=handles, **kwargs)
    # PLOT 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot.PlotPerTypeYieldTotal(data, Colors=Colors, title="Working SGFETs x Wafer", xlabel=xlabel,
                               legendTitle=legendTitle, perType="Overall", handles=handles, **kwargs)


############################
# PLOTS PER TYPES
###########################
def PlotsPerTypes(GrBase, arguments, Colors=None, legendTitle=None, xlabel=None, **kwargs):
    # DATABASE SEARCH ####################################################################################
    GrTypes, ResultsParams = s.DBSearchPerType(GrBase, arguments, **kwargs)
    # DATA CLASSIFICATION ################################################################################
    Results = g.DataClassification(GrTypes, arguments, ResultsParams)
    # PLOTTING ######
    handles = list((Patch(color=Colors[i], label=sorted(list(GrTypes.keys()))[i])
                    for i in range(0, len(list(GrTypes.keys())))))

    plot.PlotResults(Results, arguments, Colors=Colors, handles=handles, xlabel=xlabel,
                     legendTitle=legendTitle, **kwargs)

    data = Results['arg5']
    # PLOT 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot.PlotPerTypeNoise(data, legendTitle=legendTitle, Colors=Colors, xlabel=xlabel, handles=handles,
                          perType="x Type", **kwargs)
    # PLOT 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot.PlotPerTypeYield(data, legendTitle=legendTitle, Colors=Colors, xlabel=xlabel, title="Working SGFETs x Type",
                          perType="x Type", handles=handles, **kwargs)
    # PLOT 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot.PlotPerTypeYieldTotal(data, legendTitle=legendTitle, Colors=Colors, xlabel=xlabel,
                               title="Working SGFETs x Type",
                               perType="Overall", handles=handles, **kwargs)


def main():
    plt.close('all')

    Wafers1 = (
        'B12708W2',  # (in vivo Rob, slices Mavi) Very good
        #'B12142W46',  # (in vivo Rob) # High doping
        # 'B12142W15',  # (in vivo Rob)
        # 'B11870W8',  # (IDIBAPS implants)
        # 'B11601W4',
    )
    Wafers2 = (
        'B12708W2',  # (in vivo Rob, slices Mavi) Very good
        'B12142W46',  # (in vivo Rob) # High doping
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
         'Range': (200, 500),  # Range of allowed values, (Min, Max)
         'Name': 'UD0y',
         'Units': 'mV'},

        {'Param': 'GMV',  # Parameter to evaluate
         'Range': (1e-1, 10),  # Range of allowed values, (Min, Max)
         'ParArgs': {'Vgs': -0.1,  # Bias point to evaluate
                     'Vds': None,
                     'Ud0Norm': True,
                     'Units': "mS/V",
                     },
         'Name': 'GMVy'},

        {'Param': 'Vrms',  # Parameter to evaluate
         'ParArgs': {
             'Vgs': -0.1,  # Bias point to evaluate
             'Vds': None,
             'Ud0Norm': True,
             'NFmin': 10,
             'NFmax': 1000,
             'Units': "uV",
         },
         'Range': (5, 60),  # Range of allowed values, (Min, Max)
         'Units': pq.microvolt / pq.S
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

    arguments1 = {
        'arg0': {
            'Param': 'Vrms',
            'Vgs': -0.1,
            'Ud0Norm': True,
            'yscale': 'log',
            'NFmin': 10,
            'NFmax': 1000,
            'Units': pq.microvolt,
            'title': 'Vrms',
        },
        'arg1': {
            'Param': 'Rds',
            'Vgs': 0,
            'Ud0Norm': True,
            'yscale': 'log',
            'Units': pq.ohm,
            'title': 'Rds',
        },
        'arg2': {
            'Param': 'GMV',
            'Vgs': -0.1,
            'Ud0Norm': True,
            'yscale': 'log',
            'Units': 'uS/V',
            'title': 'GMV',
        },
        'arg3': {
            'Param': 'Ud0',
            'Units': pq.microvolt,
            'title': 'Ud0',
        },
        'arg4': {
            'Param': 'Ids',
            'Vgs=': -0.1,
            'Ud0Norm': True,
            'yscale': 'log',
            'Units': pq.microampere,
            'title': 'Ids',
        },
        'arg5': {
            'Param': 'Vrms',
            'Vgs': -0.1,
            'Ud0Norm': True,
            'Units': pq.microvolt,
            'title': 'Vrms',
        },
    }

    kwargs1 = {
        'arguments': arguments1,
        'Colors': Colors,
        'legendTitle': "Wafers",
        'xlabel': "Types",
        'remove50Hz': True,
    }
    kwargs2 = {
        'arguments': arguments1,
        'Colors': Colors,
        'legendTitle': "Types",
        'xlabel': "Wafers",
        'remove50Hz': True,
    }

    # PLOTS ####################################################################
    # PlotsParams(GrBase3, **kwargs2)
    # PlotsPerWaferAndTypes(GrBase1, **kwargs1)
    # PlotsPerWaferAndTypes(GrBase2, **kwargs1)
    PlotsPerWaferAndTypes(GrBase3, **kwargs1)

    # PlotsPerTypes(GrBase1, **kwargs2)
    # PlotsPerTypes(GrBase2, **kwargs2)
    PlotsPerTypes(GrBase3, **kwargs2)


# """"""""""""""""""""""""""""""""""""""""""""""
#

# qty.setActive(False)
main()
plt.show()
#os.system("read")
#
#
#
