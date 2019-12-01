import gc

from matplotlib.patches import Patch

from PyGFETdb import SearchFunctions as s, AnalysisFunctions as analysis, PlotFunctions as plot


########################################################################
#
# SEARCH, GETPARAMS AND PLOT RESULTS COMBOS
#
########################################################################

############################
# SEARCH AND PLOTS PER WAFER AND TYPE
###########################
def PlotWafersPerType(GrBase, arguments, Colors=None, legendTitle=None, xlabel=None, **kwargs):
    print(' ')
    print('******************************************************************************')
    print('******* PLOTS PER WAFER AND TYPE *********************************************')
    print('******************************************************************************')
    print(' ')
    # DATABASE SEARCH ####################################################################################
    GrWs, ResultsParams = s.DBSearchPerWaferAndType(GrBase, arguments, **kwargs)
    # DATA CLASSIFICATION ################################################################################
    Results = s.DataClassification(GrWs, arguments, ResultsParams)
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


####################################
# SEARCH AND PLOT PER TYPE AND WAFER
###################################
def PlotTypesPerWafer(GrBase, arguments, Colors=None, legendTitle=None, xlabel=None, **kwargs):
    print(' ')
    print('******************************************************************************')
    print('******* PLOTS PER TYPE *******************************************************')
    print('******************************************************************************')
    print(' ')
    # DATABASE SEARCH ####################################################################################
    GrTypes, ResultsParams = s.DBSearchPerTypeAndWafer(GrBase, arguments, **kwargs)
    # DATA CLASSIFICATION ################################################################################
    Results = s.DataClassification(GrTypes, arguments, ResultsParams)
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
# SEARCH AND PLOT PSD
############################
def PlotWorkingTypes(GrBase, Plot=False, **kwargs):
    """

    :param GrBase: Conditions to search in the database
    :param Plot: if True Plots the results
    :param kwargs: {remove50Hz: bool}
    :return: None
    """
    arguments = {
        'Fpsd': {'Param': 'Fpsd', },
        'PSD': {'Param': 'PSD', 'Vds': 0.05},
        'NoA': {'Param': 'NoA'},
        'NoB': {'Param': 'NoB'},
    }
    print(' ')
    print('******************************************************************************')
    print('******* TYPE NOISE ANALYSIS **************************************************')
    print('******************************************************************************')
    print(' ')
    GrTypes, rPSD = s.DBSearchPerType(GrBase, arguments, **kwargs.get('db'))
    results = analysis.processAllPSDsPerGroup(rPSD, **kwargs.get('noise'))
    if Plot:
        plot.PlotResultsPSDPerGroup(GrTypes, results, **kwargs)

    print('Collect->', gc.collect())
    return results


def PlotWorkingWafers(GrBase, Plot=False, **kwargs):
    """

    :param GrBase: Conditions to search in the database
    :param Plot: if True Plots the results
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
    print('******* WAFER NOISE ANALYSIS *************************************************')
    print('******************************************************************************')
    print(' ')
    GrTypes, rPSD = s.DBSearchPerWafer(GrBase, arguments, **kwargs.get('db'))
    results = analysis.processAllPSDsPerGroup(rPSD, **kwargs.get('noise'))
    if Plot:
        plot.PlotResultsPSDPerGroup(GrTypes, results, **kwargs)

    print('Collect->', gc.collect())
    return results


def PlotWorkingDevices(GrBase, Plot=False, **kwargs):
    """

    :param GrBase: Conditions to search in the database
    :param Plot: if True Plots the results
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
    print('******* DEVICE NOISE ANALYSIS ************************************************')
    print('******************************************************************************')
    print(' ')
    GrTypes, rPSD = s.DBSearchPerDevice(GrBase, arguments, **kwargs.get('db'))
    results = analysis.processAllPSDsPerGroup(rPSD, **kwargs.get('noise'))
    if Plot:
        plot.PlotResultsPSDPerGroup(GrTypes, results, **kwargs)

    print('Collect->', gc.collect())
    return results


def PlotWorkingTypesPerWafer(GrBase, Plot=False, **kwargs):
    """

    :param GrBase: Conditions to search in the database
    :param Plot: if True Plots the results
    :param kwargs: {remove50Hz: bool}
    :return: None
    """
    arguments = {
        'Fpsd': {'Param': 'Fpsd', },
        'PSD': {'Param': 'PSD', 'Vds': 0.05},
        'NoA': {'Param': 'NoA'},
        'NoB': {'Param': 'NoB'},
    }
    print(' ')
    print('******************************************************************************')
    print('******* NOISE ANALYSIS *******************************************************')
    print('******************************************************************************')
    print(' ')
    GrTypes, rPSD = s.DBSearchPerTypeAndWafer(GrBase, arguments, **kwargs.get('db'))
    results = analysis.processAllPSDsPerSubgroup(GrTypes, rPSD, **kwargs.get('noise'))
    if Plot:
        plot.PlotResultsPSDPerSubgroup(GrTypes, results, rPSD, **kwargs)

    print('Collect->', gc.collect())
    return results


def PlotWorkingWafersAndDevices(GrBase, Plot=False, **kwargs):
    """

    :param GrBase: Conditions to search in the database
    :param Plot: if True Plots the results
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
    GrTypes, rPSD = s.DBSearchPerWaferAndDevice(GrBase, arguments, **kwargs.get('db'))
    results = analysis.processAllPSDsPerSubgroup(GrTypes, rPSD, **kwargs.get('noise'))
    if Plot:
        plot.PlotResultsPSDPerSubgroup(GrTypes, results, rPSD, **kwargs)

    print('Collect->', gc.collect())
    return results


def PlotWorkingDevicesAndTrts(GrBase, Plot=False, **kwargs):
    """

    :param GrBase: Conditions to search in the database
    :param Plot: if True Plots the results
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
    GrTypes, rPSD = s.DBSearchPerDeviceAndTrt(GrBase, arguments, **kwargs.get('db'))
    results = analysis.processAllPSDsPerSubgroup(GrTypes, rPSD, **kwargs.get('noise'))
    if Plot:
        plot.PlotResultsPSDPerSubgroup(GrTypes, results, rPSD, **kwargs)

    print('Collect->', gc.collect())
    return results
