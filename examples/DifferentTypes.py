#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: dragc

"""
import gc

import matplotlib.pyplot as plt
import quantities as pq

import PyGFETdb.SearchAndPlots as plot


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

    kwargs3 = {  # Per Trt
        'db': {
            'remove50Hz': True
        },
        'noise': {  # Important parameters
            # PSD   # --------------------
            'fluctuation': 43e-3,
            'gradient': 5e-19,  # <-- max width of the deviation
            'peak': 0.58,
            'gradientmean': 0.5,  # <-- max slope of the PSD

            # Mean PSD
            'meanpsd':False,
            'meanfluctuation': 0.14,
            'meangradient': 2.65e-22,  # <-- max peak of the mean PSD
            'meanpeak': 0.4,
            'meangradientmean': 5e-18,  # <-- max slope of the mean PSD

            # Noise fit
            'fiterror': 0.5,
            'fitgradient': 5e-18,

            # Debug
            'printbad': False,
            'printok': False,

            # Optimization
            'HaltOnFail': False # TODO: Give correct results when halting
        },
        # 'PlotSuperMean': True
        'PlotMean': True,
        'PlotStd': True,
        'PlotNoise': True,
        'PlotOnlyWorking': True,  # All the OK
        # 'PlotOnlyFit': True,      # All the Fit even if not OK
        # 'PlotOnlyPerfect': True  # Only the Working and Fit
    }

    # PLOTS ####################################################################

    # ####### INTERESTING PLOTS #######################
    plot.PlotWorkingDevices(GrBase3, Plot=True, PlotSuperMean=True, **kwargs3)
    # plot.PlotWorkingTrts(GrBase1, Plot=True, **kwargs3)

    # ####### ALL THE PLOTS FOR ALL THE WAFERS #######
    # plot.PlotWafersPerType(GrBase3,**kwargs1)
    # plot.PlotTypesPerWafer(GrBase3,**kwargs2)

    # plot.PlotWorkingTrts(GrBase3, Plot=True, **kwargs3)

    # plot.PlotWorkingWafers(GrBase3, Plot=True, PlotSuperMean=True, **kwargs3)
    # plot.PlotWorkingWafersAndDevices(GrBase3, Plot=True, **kwargs3)

    # plot.PlotWorkingDevices(GrBase3, Plot=True, **kwargs3)
    # plot.PlotWorkingDevicesAndTrts(GrBase3, Plot=True, PlotSuperMean=True, **kwargs3)

    # plot.PlotWorkingTypes(GrBase3, Plot=True, PlotSuperMean=True, **kwargs3)
    # plot.PlotWorkingTypesPerWafer(GrBase3, Plot=True, **kwargs3)
    # plot.PlotWorkingTypesPerTrt(GrBase3,Plot=True,**kwargs3)

# """"""""""""""""""""""""""""""""""""""""""""""
# END MAIN


main()
plt.show()
