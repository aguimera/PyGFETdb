#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: dragc

"""
import gc
import sys

import matplotlib.pyplot as plt
import quantities as pq
from matplotlib.patches import Patch

import PyGFETdb.AnalysisFunctions as analysis
import PyGFETdb.SearchFunctions as search
from PyGFETdb import PlotFunctions as plot


############################
# PLOTS PER WAFER AND TYPES
###########################
def PlotsPerWaferAndTypes(GrBase, arguments, Colors=None, legendTitle=None, xlabel=None, **kwargs):
    print(' ')
    print('******************************************************************************')
    print('******* PLOTS PER WAFER AND TYPE *********************************************')
    print('******************************************************************************')
    print(' ')
    # DATABASE SEARCH ####################################################################################
    GrWs, ResultsParams = search.DBSearchPerWaferAndType(GrBase, arguments, **kwargs)
    # DATA CLASSIFICATION ################################################################################
    Results = search.DataClassification(GrWs, arguments, ResultsParams)
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
    plot.PlotPerTypeYieldTotal(data, Colors=Colors, title="Working SGFETs x Wafer", xlabel="Wafers",
                               legendTitle=legendTitle, perType="Overall", handles=handles, **kwargs)

    print('Collect->', gc.collect())


############################
# PLOTS PER TYPES
###########################
def PlotsPerTypes(GrBase, arguments, Colors=None, legendTitle=None, xlabel=None, **kwargs):
    print(' ')
    print('******************************************************************************')
    print('******* PLOTS PER TYPE *******************************************************')
    print('******************************************************************************')
    print(' ')
    # DATABASE SEARCH ####################################################################################
    GrTypes, ResultsParams = search.DBSearchPerType(GrBase, arguments, **kwargs)
    # DATA CLASSIFICATION ################################################################################
    Results = search.DataClassification(GrTypes, arguments, ResultsParams)
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
    plot.PlotPerTypeYieldTotal(data, legendTitle=legendTitle, Colors=Colors, xlabel="Types",
                               title="Working SGFETs x Type",
                               perType="Overall", handles=handles, **kwargs)

    print('Collect->', gc.collect())


#############################
# PLOTS PSD
############################
def PlotsPSDperTypeAndWafer(GrBase, Plot=False, **kwargs):
    """

    :param GrBase: Conditions to search in the database
    :param kwargs: {remove50Hz: bool}
    :return: None
    """
    arguments = {
        'Fpsd': {'Param': 'Fpsd', },
        'PSD': {'Param': 'PSD', 'Vds': 0.05, },
        'NoA': {'Param': 'NoA'},
        'NoB': {'Param': 'NoB'},
    }
    print(' ')
    print('******************************************************************************')
    print('******* NOISE ANALYSIS *******************************************************')
    print('******************************************************************************')
    print(' ')
    GrTypes, rPSD = search.DBSearchPerType(GrBase, arguments, **kwargs.get('db'))
    results = analysis.processAllPSDs(GrTypes, rPSD, **kwargs.get('noise'))
    if Plot:
        plot.PlotResultsPSDPerType(GrTypes, results, rPSD, PlotStd=True)

    print('Collect->', gc.collect())
    return results


def PlotsPSDperWaferAndDevice(GrBase, Plot=False, **kwargs):
    """

    :param GrBase: Conditions to search in the database
    :param kwargs: {remove50Hz: bool}
    :return: None
    """
    arguments = {
        'Fpsd': {'Param': 'Fpsd', },
        'PSD': {'Param': 'PSD', 'Vds': 0.05, },
        'NoA': {'Param': 'NoA'},
        'NoB': {'Param': 'NoB'},
    }
    print(' ')
    print('******************************************************************************')
    print('******* NOISE ANALYSIS *******************************************************')
    print('******************************************************************************')
    print(' ')
    GrTypes, rPSD = search.DBSearchPerWaferAndDevice(GrBase, arguments, **kwargs.get('db'))
    results = analysis.processAllPSDs(GrTypes, rPSD, **kwargs.get('noise'))
    if Plot:
        plot.PlotResultsPSDPerType(GrTypes, results, rPSD)

    print('Collect->', gc.collect())
    return results


def PlotsPSDperDevice(GrBase, Plot=False, **kwargs):
    """

    :param GrBase: Conditions to search in the database
    :param kwargs: {remove50Hz: bool}
    :return: None
    """
    arguments = {
        'Fpsd': {'Param': 'Fpsd', },
        'PSD': {'Param': 'PSD', 'Vds': 0.05, },
        'NoA': {'Param': 'NoA'},
        'NoB': {'Param': 'NoB'},
    }
    print(' ')
    print('******************************************************************************')
    print('******* NOISE ANALYSIS *******************************************************')
    print('******************************************************************************')
    print(' ')
    GrTypes, rPSD = search.DBSearchPerDevice(GrBase, arguments, **kwargs.get('db'))
    results = analysis.processAllPSDs(GrTypes, rPSD, **kwargs.get('noise'))
    if Plot:
        plot.PlotResultsPSDPerType(GrTypes, results, rPSD)

    print('Collect->', gc.collect())
    return results


def PlotsPSDperDeviceAndTrt(GrBase, Plot=False, **kwargs):
    """

    :param GrBase: Conditions to search in the database
    :param kwargs: {remove50Hz: bool}
    :return: None
    """
    arguments = {
        'Fpsd': {'Param': 'Fpsd', },
        'PSD': {'Param': 'PSD', 'Vds': 0.05, },
        'NoA': {'Param': 'NoA'},
        'NoB': {'Param': 'NoB'},
    }
    print(' ')
    print('******************************************************************************')
    print('******* NOISE ANALYSIS *******************************************************')
    print('******************************************************************************')
    print(' ')
    GrTypes, rPSD = search.DBSearchPerDeviceAndTrt(GrBase, arguments, **kwargs.get('db'))
    results = analysis.processPSDsPerTrt(GrTypes, rPSD, **kwargs.get('noise'))
    if Plot:
        plot.PlotResultsPSDPerType(GrTypes, results, rPSD)

    print('Collect->', gc.collect())
    return results

############################
# MAIN
###########################

def main():
    plt.close('all')
    plt.rcParams['figure.max_open_warning'] = False
    gc.enable()

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
    kwargs3 = {  # Per Group
        'db': {
            'remove50Hz': True
        },
        'noise': {
            'tolerance': 0.75,
            'errortolerance': 0.75,
            'gradtolerance': 0.09
        }
    }
    kwargs4 = {  # Per Trt
        'db': {
            'remove50Hz': True
        },
        'noise': {
            'tolerance': 0.75,
            'errortolerance': 0.75,
            'gradtolerance': 0.09
        }
    }
    # PLOTS ####################################################################

    ####### ONE WAFER ##########################
    # PlotsPerWaferAndTypes(GrBase1, **kwargs1)
    # PlotsPerTypes(GrBase1, **kwargs2)

    # PlotsPSDperTypeAndWafer(GrBase1, Plot=True, **kwargs3)
    # PlotsPSDperWaferAndDevice(GrBase1, Plot=True, **kwargs3)
    # PlotsPSDperDeviceAndTrt(GrBase1, Plot=True, **kwargs4)

    ######## 2 WAFERS #####################
    # PlotsPerWaferAndTypes(GrBase2, **kwargs1)
    # PlotsPerTypes(GrBase2, **kwargs2)

    # PlotsPSDperTypeAndWafer(GrBase2, Plot=True, **kwargs3)
    PlotsPSDperWaferAndDevice(GrBase2, Plot=True, **kwargs3)
    # PlotsPSDperDeviceAndTrt(GrBase2, Plot=True, **kwargs4)

    ######## ALL THE WAFERS #####################
    # PlotsPerWaferAndTypes(GrBase3, **kwargs1)
    # PlotsPerTypes(GrBase3, **kwargs2)

    # PlotsPSDperTypeAndWafer(GrBase3, Plot=True, **kwargs3)
    # PlotsPSDperWaferAndDevice(GrBase3, Plot=True, **kwargs3)
    # PlotsPSDperDeviceAndTrt(GrBase3, **kwargs4)


# """"""""""""""""""""""""""""""""""""""""""""""
# END MAIN


main()
plt.show()
# os.system("read")
#
#
#
sys.exit()
