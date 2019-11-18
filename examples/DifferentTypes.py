#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: dragc

"""
import os

import matplotlib.pyplot as plt
import numpy as np
import quantities as pq
from matplotlib.patches import Patch

import PyGFETdb.DBSearch as DbSe
import PyGFETdb.GlobalFunctions as g
import PyGFETdb.SearchFunctions as s
from PyGFETdb import PlotFunctions as plot
from PyGFETdb import qty, multithrds, Multiprocessing as mp
from PyGFETdb.NoiseModel import Fnoise

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

def PlotsPSDperType(GrBase, **kwargs):
    """

    :param GrBase: Conditions to search in the database
    :param kwargs: {remove50Hz: bool}
    :return: None
    """
    arguments = {
        'Fpsd': {
            'Param': 'Fpsd',
        },
        'Vds': {
            'Param': 'Vds',
        }
    }
    qty.setActive(False)
    GrTypes, ResultsParams = s.DBSearchPerWaferAndType(GrBase, arguments, **kwargs)
    for nWf, vWf in GrTypes.items():
        Fpsd = ResultsParams[nWf].get('Fpsd')
        Vds = ResultsParams[nWf].get('Vds')
        for nType, vType in Fpsd.items():
            arguments2 = {
                'PSD': {
                    'Param': 'PSD',
                    'Vds': Vds[nType][0][0],
                },
                'NoA': {'Param': 'NoA'},
                'NoB': {'Param': 'NoB'},
            }
            grPSD, rPSD = s.DBSearchPerWaferAndType(GrBase, arguments2, **kwargs)

            PSD = rPSD[nWf]['PSD'][nType]
            NoA = np.array(rPSD[nWf]['NoA'][nType])
            NoB = np.array(rPSD[nWf]['NoB'][nType])

            Fpsd = vType[0:len(PSD)]
            Fpsd2 = np.array(Fpsd).reshape((1, len(PSD)))
            NoA = np.mean(NoA.transpose(), 1)
            NoA = NoA.reshape((1, NoA.size))
            NoB = np.mean(NoB.transpose(), 1)
            NoB = NoB.reshape((1, NoB.size))

            noise = Fnoise(Fpsd2, NoA[:, len(PSD)], NoB[:, len(PSD)])

            fig, ax = plt.subplots()
            plot.PlotMeanStd(Fpsd, PSD, ax, xscale='log', yscale='log')

            ax.loglog(Fpsd2.transpose(), noise.transpose(), '--')
            ax.set_xlabel("Frequency [Hz]")
            ax.set_ylabel('PSD [A^2/Hz]')

            title = "PSDs for Type {}".format(nType)
            plt.title(title)
    qty.setActive(True)


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
         'Units': "mS/V",
         'ParArgs': {'Vgs': -0.1,  # Bias point to evaluate
                     'Vds': None,
                     'Ud0Norm': True,
                     'Units': "mS/V",  #
                     },
         'Name': 'GMVy'},

        {'Param': 'Vrms',  # Parameter to evaluate
         'ParArgs': {
             'Vgs': -0.1,  # Bias point to evaluate
             'Vds': None,
             'Ud0Norm': True,
             'NFmin': 10,
             'NFmax': 1000,
             # 'Units': "uV",
         },
         'Range': (5e-6, 60e-6),  # Range of allowed values, (Min, Max)
         }
    ]

    DataSelectionConfig2 = [
        {'Param': 'Ud0',  # Parameter to evaluate
         'Range': (200e-3, 500e-3),  # Range of allowed values, (Min, Max)
         'Name': 'UD0y', },
        {'Param': 'GMV',  # Parameter to evaluate
         'Range': (1e-4, 10e-3),  # Range of allowed values, (Min, Max)
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
         'Range': (5e-6, 60e-6),  # Range of allowed values, (Min, Max)
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

    GrBase4 = {'Conditions': Conditions1,
               'Table': CharTable,
               'Last': True,
               'DataSelectionConfig': DataSelectionConfig2
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
            'Units': "kohm",
            'title': 'Rds',
        },
        'arg2': {
            'Param': 'GMV',
            'Vgs': -0.1,
            'Ud0Norm': True,
            'yscale': 'log',
            'Units': 'mS/V',
            'title': 'GMV',
        },
        'arg3': {
            'Param': 'Ud0',
            'Units': pq.millivolt,
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
            'title': 'Vrms',
            'Units': 'uV'
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
    kwargs3 = {
        'remove50Hz': True,
    }
    # PLOTS ####################################################################
    # PlotsParams(GrBase3, **kwargs2)

    # PlotsPerWaferAndTypes(GrBase1, **kwargs1)
    # PlotsPerWaferAndTypes(GrBase2, **kwargs1)
    # PlotsPerWaferAndTypes(GrBase3, **kwargs1)

    # PlotsPerTypes(GrBase1, **kwargs2)
    # PlotsPerTypes(GrBase2, **kwargs2)
    #PlotsPerTypes(GrBase3, **kwargs2)

    PlotsPSDperType(GrBase4, **kwargs3)



# """"""""""""""""""""""""""""""""""""""""""""""
#

main()
plt.show()
os.system("read")
#
#
#
