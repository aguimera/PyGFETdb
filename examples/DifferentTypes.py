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
    )
    Wafers2 = (
        'B12708W2',  # (in vivo Rob, slices Mavi) Very good
        'B12142W46',  # (in vivo Rob) # High doping
    )

    Wafers3 = (
        'B12708W2',  # (in vivo Rob, slices Mavi) Very good
        'B12142W46',  # (in vivo Rob) # High doping
        'B12142W15',  # (in vivo Rob)
        'B11870W8',  # (IDIBAPS implants)
        'B11601W4',
    )

    # With Quantities
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

    # Without Quantities
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

    Conditions2 = Conditions1
    Conditions2.update({'Wafers.Name = ': Wafers2})

    Conditions3 = Conditions1
    Conditions3.update({'Wafers.Name = ': Wafers3})

    CharTable = 'ACcharacts'

    GrBase1 = {'Conditions': Conditions1,
               'Table': CharTable,
               'Last': True,
               'DataSelectionConfig': DataSelectionConfig
               }

    GrBase2 = GrBase1
    GrBase2.update({'Conditions': Conditions2})

    GrBase3 = GrBase1
    GrBase3.update({'Conditions': Conditions3})

    GrBase4 = GrBase1
    GrBase4.update({'DataSelectionConfig': DataSelectionConfig2})

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

    # PLOT GLOBALS FOR PLOTS PER WAFER
    kwargs1 = {
        'arguments': arguments1,
        'Colors': Colors,
        'legendTitle': "Wafers",
        'xlabel': "Types",
        'remove50Hz': True,
    }

    # PLOT GLOBALS FOR PLOTS PER TYPES
    kwargs2 = kwargs1
    kwargs2.update(
        {
            'legendTitle': "Types",
            'xlabel': "Wafers",
        }
    )

    # PLOT GLOBALS FOR PLOTS PSD PER GROUP AND SUBGROUP
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
        'PlotStd': False,
        'PlotNoise': True,
    }

    # PLOT GLOBALS FOR PLOTS PSD PER GROUP ONLY
    kwargs4 = kwargs3
    kwargs4.update(
        {'PlotStd': True}
    )

    # PLOTS ####################################################################

    # ###### ONE WAFER ##########################
    # plot.PlotsPerWaferAndTypes(GrBase1, **kwargs1)
    # plot.PlotsPerTypes(GrBase1, **kwargs2)

    plot.PlotsPSDperWafer(GrBase1, Plot=True, **kwargs4)
    # plot.PlotsPSDperDevice(GrBase1, Plot=True, **kwargs4)
    # plot.PlotsPSDperType(GrBase1,Plot=True,**kwargs4)

    # plot.PlotsPSDperTypeAndWafer(GrBase1, Plot=True, **kwargs3)
    # plot.PlotsPSDperWaferAndDevice(GrBase1, Plot=True, **kwargs3)
    # plot.PlotsPSDperDeviceAndTrt(GrBase1, Plot=True, **kwargs3)

    # ####### 2 WAFERS #####################
    # plot.PlotsPerWaferAndTypes(GrBase2, **kwargs1)
    # plot.PlotsPerTypes(GrBase2, **kwargs2)

    # plot.PlotsPSDperWafer(GrBase2, Plot=True, **kwargs3)
    # plot.PlotsPSDperDevice(GrBase2,Plot=True,**kwargs4)
    # plot.PlotsPSDperType(GrBase2,Plot=True,**kwargs3)

    # plot.PlotsPSDperTypeAndWafer(GrBase2, Plot=True, **kwargs3)
    # plot.PlotsPSDperWaferAndDevice(GrBase2, Plot=True, **kwargs3)
    # plot.PlotsPSDperDeviceAndTrt(GrBase2, Plot=True, **kwargs3)

    # ####### ALL THE WAFERS #####################
    # plot.PlotsPerWaferAndTypes(GrBase3, **kwargs1)
    # plot.PlotsPerTypes(GrBase3, **kwargs2)

    # plot.PlotsPSDperWafer(GrBase3, Plot=True, **kwargs3)
    # plot.PlotsPSDperDevice(GrBase3, Plot=True, **kwargs4)
    # plot.PlotsPSDperType(GrBase3,Plot=True,**kwargs3)

    # plot.PlotsPSDperTypeAndWafer(GrBase3, Plot=True, **kwargs3)
    # plot.PlotsPSDperWaferAndDevice(GrBase3, Plot=True, **kwargs3)
    # plot.PlotsPSDperDeviceAndTrt(GrBase3, **kwargs3)


# """"""""""""""""""""""""""""""""""""""""""""""
# END MAIN


main()
plt.show()
# os.system("read")
#
#
#
sys.exit()
