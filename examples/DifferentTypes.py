#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: dragc

"""
import gc
import sys

import matplotlib.pyplot as plt
import quantities as pq

import PyGFETdb.PlotFunctions as plot


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

    kwargs3 = {
        'db': {
            'remove50Hz': True
        },
        'noise': {
            'fluctuation': 38,  # >
            'peak': 58.95,  # >
            'gradient': 2e5,  # <=
            'fiterror': 90,  # >
            'fitgradient': 1e3,  # <=
            'normalization': 1e-22
        },
        'PlotMean': True,
        'PlotStd': True,
        'PlotNoise': True,
    }

    # PLOTS ####################################################################

    ####### ONE WAFER ##########################
    # plot.PlotsPerWaferAndTypes(GrBase1, **kwargs1)
    # plot.PlotsPerTypes(GrBase1, **kwargs2)

    # plot.PlotsPSDperTypeAndWafer(GrBase1, Plot=True, **kwargs3)
    # plot.PlotsPSDperWaferAndDevice(GrBase1, Plot=True, **kwargs3)
    # plot.PlotsPSDperDeviceAndTrt(GrBase1, Plot=True, **kwargs3)

    ######## 2 WAFERS #####################
    # plot.PlotsPerWaferAndTypes(GrBase2, **kwargs1)
    # plot.PlotsPerTypes(GrBase2, **kwargs2)

    # plot.PlotsPSDperTypeAndWafer(GrBase2, Plot=True, **kwargs3)
    # plot.PlotsPSDperWaferAndDevice(GrBase2, Plot=True, **kwargs3)
    # plot.PlotsPSDperDeviceAndTrt(GrBase2, Plot=True, **kwargs3)

    ######## ALL THE WAFERS #####################
    # plot.PlotsPerWaferAndTypes(GrBase3, **kwargs1)
    # plot.PlotsPerTypes(GrBase3, **kwargs2)

    # plot.PlotsPSDperTypeAndWafer(GrBase3, Plot=True, **kwargs3)
    # plot.PlotsPSDperWaferAndDevice(GrBase3, Plot=True, **kwargs3)
    plot.PlotsPSDperDeviceAndTrt(GrBase3, Plot=True, **kwargs3)


# """"""""""""""""""""""""""""""""""""""""""""""
# END MAIN


main()
plt.show()
# os.system("read")
#
#
#
sys.exit()
