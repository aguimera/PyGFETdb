# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Plot Functions that do not fit in the previous files.

"""

import matplotlib.pyplot as plt
import numpy as np
import quantities as pq

from PyGFETdb import qty


########################################################################
#
# PLOT UTILITY FUNCTIONS
#
########################################################################

def _BoxplotValsGroup(ax, Col, iGroup, Vals, **kwargs):
    """

    :param ax:
    :param Col: Color
    :param iGroup: Position to plot
    :param Vals: Values to plot
    :param kwargs: Aesthetics arguments for the plot
    :return: None
    """
    bplt = ax.boxplot(Vals,
                      positions=(iGroup,),
                      patch_artist=True,  # fill with color
                      widths=0.75,
                      sym='+',
                      labels=('',),
                      #                      notch=True,
                      )

    for element in ('boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps',):
        plt.setp(bplt[element], color=Col)

    for fl in bplt['fliers']:
        fl.set_markeredgecolor(Col)

    for patch in bplt['boxes']:
        patch.set(facecolor=Col)
        patch.set(alpha=0.5)

    return bplt


def _closeBoxplotValsGroup(ax, xPos, xLab, xlabel=None, ylabel=None, title='', legendtitle='', handles=None, **kargs):
    """

    :param ax:
    :param xPos:
    :param xLab:
    :param xlabel: a label for the x axis
    :param ylabel: a label for the y axis
    :param title: a label for the title of the plot
    :param legendtitle: a label for the title of the legend
    :param handles: handles to update the legend
    :param kargs: Aesthetics arguments for the plot

    :return: None
    """
    plt.xticks(xPos, xLab, rotation=75, fontsize='small')
    ax.set_ylabel(ylabel=ylabel, fontsize='large')

    ax.set_xlabel(xlabel, fontsize='large')
    ax.set_title(label=title, fontsize='large')
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0 + chartBox.y0 * 1.65, chartBox.width, chartBox.height * 0.77])


def _PlotValsGroup(Ax, xLab, xPos, iGr, Grn, vals,
                   Col=None, Boxplot=False, ParamUnits=None, **kwargs):
    """
        Auxiliary function to plot the values of a group.
    :param Ax:
    :param xLab:
    :param xPos:
    :param iGr: Number of the group
    :param Grn: Name of the group
    :param vals: The values to plot
    :param Col: Color
    :param Boxplot: if True plots a Boxplot
    :param ParamUnits: Units of the parameters to plot
    :param kwargs:Aesthetics arguments for the plot
    :return: None
    """
    if vals is not None:  # and len(vals) >0:
        if Boxplot:
            Ax.boxplot(vals.transpose(), positions=(iGr + 1,))
            xPos.append(iGr + 1)
        else:
            Ax.plot(np.ones(len(vals)) * iGr, vals, '.', label=(Grn,), color=Col)
            xPos.append(iGr)
            xLab.append(Grn)
    else:
        print('Empty data for: ', Grn)


def _closePlotValsGroup(Ax, xLab, xPos, qtys=None, ParamUnits=None,
                        title=None, legendTitle=None, handles=None, xlabel='', **kwargs):
    """

    :param Ax:
    :param xLab:
    :param xPos:
    :param qtys: Quantities to extract the units automatically
    :param ParamUnits: Units of the parameters to plot
    :param title: Title of the plot
    :param legendTitle: Title of the Legend
    :param handles: Matplotlib handles to plot the Legend
    :param xlabel: String of the title for the x axis
    :param kwargs: Aesthetics arguments for the plot
    :return: None
    """
    units = None
    qtyunits = None
    if qty.isActive():
        units = kwargs.get('Units')
        if qtys is not None:
            qtys = qty.rescaleFromKey(qtys, units)  # Consistency Check
            qtyunits = qty.getQuantityUnits(qtys)
            if qtyunits is not None and (type(qtyunits) is pq.UnitQuantity or type(qtyunits) is pq.Quantity):
                qtyunits = qtyunits.dimensionality.latex
        elif units is not None and (type(units) is pq.UnitQuantity or type(units) is pq.Quantity):
            units = units.dimensionality.latex

    plt.xticks(xPos, xLab, rotation=75, fontsize='small')

    Ax.set_xlabel(xlabel, fontsize='large')

    param = kwargs.get('Param')

    if qtyunits is not None:
        Ax.set_ylabel(param + '[' + qtyunits + ']')
    elif units is not None:
            Ax.set_ylabel(param + '[' + units + ']')
    elif ParamUnits is not None:
        Ax.set_ylabel(param + '[' + ParamUnits + ']')
    else:
        Ax.set_ylabel(param)

    Ax.grid()
    Ax.ticklabel_format(axis='y', style='sci', scilimits=(2, 2))
    if len(xPos) > 1:
        Ax.set_xlim(min(xPos) - 0.5, max(xPos) + 0.5)
    if 'Vgs' in kwargs and 'Vds' in kwargs:
        title = 'Vgs {} Vds {}'.format(kwargs['Vgs'], kwargs['Vds'])
        plt.title(title)
        plt.tight_layout()
    if 'xscale' in list(kwargs.keys()):
        Ax.set_xscale(kwargs['xscale'])
    if 'yscale' in list(kwargs.keys()):
        Ax.set_yscale(kwargs['yscale'])

    Ax.set_title(title, fontsize='large')
    chartBox = Ax.get_position()
    Ax.set_position([chartBox.x0, chartBox.y0 + chartBox.y0 * 1.65, chartBox.width, chartBox.height * 0.77])
    if legendTitle or handles:
        if legendTitle:
            legend = legendTitle
        else:
            legend = ''
        if handles:
            colors = handles
        else:
            colors = None
        Legend(Ax, legend, colors)


def Legend(Ax, legend, handles):
    """

    :param Ax:
    :param legend: Title of the Legend
    :param handles: Matplotlib handles to plot the Legend
    :return: None
    """
    chartBox = Ax.get_position()
    Ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.8, chartBox.height])
    Ax.legend(title=legend, handles=handles, loc='upper right', fontsize='x-small',
              bbox_to_anchor=(0.86, 0.5, 0.5, 0.5), shadow=True, ncol=1)


def PlotMeanStd(Valx, Valy, Ax=None,
                PlotStd=True,
                PlotOverlap=False,
                label=None,
                **kwargs):
    """

    :param Valx: Data for the x axis
    :param Valy: Data for the y axis
    :param Ax: Matplotib.Axis
    :param Color: Color to start plotting
    :param PlotStd: default true, if true plots the standard deviation
    :param PlotOverlap: if true plot over the previous plot
    :param label: label of the Data
    :param kwargs:
    :return:
    """
    scilimits = (-2, 2)
    if PlotStd:
        Color = 'r'
    if Ax is None:
        fig, Ax = plt.subplots()

    if PlotOverlap:
        if Valy is not None:
            Ax.plot(Valx, Valy, alpha=0.2)

    Valy = np.array(Valy)
    if Valy is not None and Valy.size and Valy.ndim == 2:
        try:
            avg = np.mean(Valy, axis=0)
            std = np.std(Valy, axis=0)
            Ax.plot(Valx, avg, label=label)
            if PlotStd:
                Ax.fill_between(Valx, avg - std, avg + std,
                                linewidth=0.0,
                                alpha=0.3)
        except:
            print('Error averaging in PlotStdMean().')

    if 'xscale' in list(kwargs.keys()):
        Ax.set_xscale(kwargs['xscale'])
    else:
        Ax.ticklabel_format(axis='x', style='sci', scilimits=scilimits)

    if 'yscale' in list(kwargs.keys()):
        Ax.set_yscale(kwargs['yscale'])
    else:
        Ax.ticklabel_format(axis='y', style='sci', scilimits=scilimits)

    if 'yscale' in list(kwargs.keys()):
        Ax.set_yscale(kwargs['yscale'])
    else:
        Ax.ticklabel_format(axis='y', style='sci', scilimits=scilimits)


def PlotGroup(ResultsParams, Group, arguments, handles=None, **kwargs):
    """

    :param ResultsParams: The results of a search in the database
    :param Group: A group of conditions to analyse
    :param args: Arguments for getting the parameters
    :param kwargs:
    :return: A dict of args of dicts of groupnames and parameter found in a previous search
    """
    Results = {}
    for iarg, (karg, arg) in enumerate(arguments.items()):
        Results[karg] = {}
        fig, Ax = plt.subplots()
        xLab = []
        xPos = []
        qtys = []
        for iGr, (Grn, Grc) in enumerate(sorted(Group.items())):
            argRes = ResultsParams.get(karg)
            if argRes is not None:
                ParamData = argRes.get(Grn)
                if ParamData is not None:
                    Results[karg][Grn] = ParamData
                    if qty.isActive():
                        qtys.append(ParamData)
                        ParamData = qty.flatten(ParamData)
                    ParamData = np.array(ParamData)
                    _PlotValsGroup(Ax, xLab, xPos, iGr, Grn, ParamData, **arg, **kwargs)
        _closePlotValsGroup(Ax, xLab, xPos, qtys, handles=handles, **arg, **kwargs)

    return Results


########################################################################
#
# PLOTS PER TYPE
#
########################################################################

def PlotPerTypeNoise(Results, handles=None, xlabel=None, legendTitle=None, Colors=None, perType="", **kwargs):
    """
        **Plots the noise analysis**

    :param Results: The values to plot
    :param handles: Matplotlib handles to plot the Legend
    :param xlabel: The title for the x axis
    :param legendTitle: The title for the Legend
    :param Colors: Colors to use in the plot
    :param perType: String to adjust the title of the plot
    :param kwargs: Keyword arguments
    :return: None
    """
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
                vals = qty.flatten(vals)
                _BoxplotValsGroup(ax, Col, pos, np.array(vals).transpose())
                pos += 1
    _closeBoxplotValsGroup(ax, xPos, xLab, xlabel, "[uVrms]", "Noise (10Hz-1kHz) " + perType, **kwargs)
    Legend(ax, legendTitle, handles)


def PlotPerTypeYield(Results, title=None, handles=None, xlabel=None, perType=None, Colors=None, legendTitle=None,
                     **kwargs):
    """
           **Plots the percentage of transistors yield per secondary group**

    :param Results: The values to plot
    :param title: The title of the plot
    :param handles: Matplotlib handles to plot the Legend
    :param xlabel: The title for the x axis
    :param perType: String to adjust the title of the plot
    :param Colors: Colors to use in the plot
    :param legendTitle: The title for the Legend
    :param kwargs: Keyword arguments
    :return: None
    """

    # PLOTS 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fig, ax2 = plt.subplots()
    xLab = []
    xPos = []
    pos = 0

    for iWf, (wn, dd) in enumerate(Results.items()):
        Col = Colors[iWf]
        temp = []
        for iType, (Grn, Grc) in enumerate(sorted(dd.items())):  # Param 0
            temp.append((np.array(Grc)).size)
        totalWf = np.sum(temp)

        work = []
        for iType, (Grn, Grc) in enumerate(sorted(dd.items())):  # Param 0
            work.append(qty.Divide((np.array(Grc).size), totalWf) * 100)
            xLab.append(Grn)
            xPos.append(pos)
            _BoxplotValsGroup(ax2, Col, pos, work)
            pos += 1
    _closeBoxplotValsGroup(ax2, xPos, xLab, xlabel, "Yield [%{}]".format(perType), title, **kwargs)
    Legend(ax2, legendTitle, handles)


def PlotPerTypeYieldTotal(Results, title=None, Colors=None, xlabel=None, perType=None,
                          **kwargs):
    """

           **Plots the percentage of transistors yield per primary group**

    :param Results: The values to plot
    :param title: The title of the plot
    :param Colors: Colors to use in the plot
    :param xlabel: The title for the x axis
    :param kwargs: Keyword arguments
    :return: None
    """

    # PLOT 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig, ax2 = plt.subplots()
    xLab = []
    xPos = []

    temp = []
    for nt, (typename, vtype) in enumerate(Results.items()):
        for iWf, (nWf, cWf) in enumerate(vtype.items()):
            temp.append((np.array(cWf)).size)
    total = np.sum(temp)

    work = []
    for nt, (typename, vtype) in enumerate(Results.items()):
        Col = Colors[nt]
        for iWf, (nWf, cWf) in enumerate(vtype.items()):
            xLab.append(typename)
            xPos.append(nt)
            work.append(qty.Divide((np.array(cWf)).size, total) * 100)
        _BoxplotValsGroup(ax2, Col, nt, work, **kwargs)
    _closeBoxplotValsGroup(ax2, xPos, xLab, xlabel, "Yield [% {}]".format(perType), title, **kwargs)


########################################################################
#
#  PLOT RESULTS
#
########################################################################
def PlotResults(Results, arguments, Colors=None, handles=None, xlabel=None, legendTitle=None, **kwargs):
    """

    :param Results: The results of a search in the database
    :param Group: A group of conditions to analyse
    :param args: Arguments for getting the parameters
    :param kwargs:
    :return: None
    """
    for iarg, (karg, arg) in enumerate(Results.items()):
        fig, Ax = plt.subplots()
        xLab = []
        xPos = []
        qtys = []
        icolor = 0
        pos = 0
        for iGr, (Grn, Grc) in enumerate(sorted(arg.items())):
            for iGr2, (Grn2, ParamData) in enumerate(sorted(Grc.items())):
                if ParamData is not None:
                    if qty.isActive():
                        qtys.append(ParamData)
                    ParamData = qty.flatten(ParamData)
                    ParamData = np.array(ParamData)
                    _PlotValsGroup(Ax, xLab, xPos, pos, Grn2, ParamData, Colors[icolor], **arguments[karg])
                    pos += 1
            icolor += 1
        _closePlotValsGroup(Ax, xLab, xPos, qtys, handles=handles, legendTitle=legendTitle,
                            xlabel=xlabel, **arguments[karg])


def PlotResultsPSDPerSubgroup(GrTypes, results, rPSD, **kwargs):
    """
        **Plots the results of the Noise Analysis**

    :param GrTypes: Group to plot
    :param results: results of Noise Analysis
    :param rPSD: results of a PSD search in the database
    :param PlotStd: Plot Standard Deviation and Noise
    :param PlotMean: Plot PSD Mean, if False Plot all the PSDs
    :param PlotNoise: Plot Noise Mean
    :return: None
    """
    title = ""
    for nType, vType in GrTypes.items():
        Fpsd = rPSD[nType].get('Fpsd')
        temp0 = []
        temp1 = []
        temp2 = []
        temp3 = []
        temp4 = []
        temp5 = []

        for nWf, vWf in Fpsd.items():
            r = results[nType].get(nWf)
            if r is not None:
                temp0.append(r[0])
                temp1.append(r[1])
                temp2.append(r[2])
                temp3.append(r[3])
                temp4.append(r[4])
                temp5.append(r[5])

        temp5 = np.all(temp5)  # perfect
        temp4 = np.all(temp4)  # ok

        if temp1 is not None and len(temp1) > 0:
            PlotPSDMean(
                        temp0[0],  # Fpsd
                        temp1,  # PSD
                        temp3,  # noise
                temp5,  # perfect
                temp4,  # ok
                nType,
                        title=title,
                        **kwargs)


def PlotPSD(ax, Fpsd, PSD, noisefit, perfect, ok, nType, PlotStd=True, PlotMean=True,
            PlotNoise=False, **kwargs):
    """

    :param Fpsd: Data of the x axis
    :param PSD:  Data of the y axis
    :param noise: Data of the y axis for noise fitting
    :param perfect: Boolean to approve the analysis
    :param nType: Name of the Type of Trt
    :param PlotStd: Plot Standard Deviation
    :param PlotMean: Plot PSD Mean, if False Plot all the PSDs
    :param PlotNoise: Plot Noise Mean
    :return: None
    """

    noise = np.array(noisefit)

    for i, item in enumerate(PSD):
        item = np.array(item)
        for n in range(0, item.ndim - 3):
            item = reduce(item)
        for n in range(0, item.ndim - 3):
            item = reduce(item)
        if item is not None and item.ndim == 3 and noise.ndim == 3:
            item = item.reshape((noise.shape[1], item.shape[2], noise.shape[2]))
        else:
            item = item.reshape((item.shape[1], item.shape[0], item.shape[2]))
        for i2, item0 in enumerate(item):
            if item0 is not None and item0.ndim == 2:
                if PlotMean:
                    PlotMeanStd(Fpsd, item0, ax, xscale='log', yscale='log', PlotStd=PlotStd)
                else:
                    for i3, item2 in enumerate(item0):
                        PlotMeanStd(Fpsd, item2, ax, PlotOverlap=True, xscale='log', yscale='log',
                                    PlotStd=PlotStd)
        if PlotNoise and noise is not None:
            if PlotMean or len(noise) != len(PSD):
                noisei = processNoiseForPlotPSD(noise)
            else:
                noisei = processNoiseForPlotPSD(noise[i])

            ax.loglog(Fpsd, noisei, '--')

    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel('PSD [A^2/Hz]')

    title = "PSDs {}/{} for {}".format("OK" if perfect else "NOK", "FIT" if ok else "NOK", nType)
    plt.title(title)


def PlotResultsPSDPerGroup(GrTypes, results, **kwargs):
    """
        **Plots the results of the Noise Analysis**

    :param GrTypes: Group to plot
    :param results: results of Noise Analysis
    :param rPSD: results of a PSD search in the database
    :param PlotStd: Plot Standard Deviation and Noise
    :param PlotMean: Plot PSD Mean, if False Plot all the PSDs
    :param PlotNoise: Plot Noise Mean
    :return: None
    """
    for nType, vType in GrTypes.items():
        r = results.get(nType)
        if r is not None:
            PlotPSDMean(r[0][0],  # Fpsd
                        r[1],  # PSD
                        r[3],  # noise
                        r[5],  # perfect
                        r[4],  # ok
                        nType,
                        **kwargs)


def PlotPSDMean(Fpsd, PSD, noisefit, perfect, ok, nType, PlotSuperMean=None, PlotOnlyWorking=False,
                PlotOnlyFit=False, PlotOnlyPerfect=False, **kwargs):
    """

    :param Fpsd: Data of the x axis
    :param PSD:  Data of the y axis
    :param noisefit: Data of the y axis for noise fitting
    :param perfect: Boolean to approve the analysis
    :param nType: Name of the Type of Trt
    :param PlotSuperMean: Plot PSD Mean, if False Plot all the PSDs
    :param plotonlyworking: Plots only the working items (Devices, Trts, etc).
    :return: None
    """
    if np.all(perfect) and np.all(ok) and PlotOnlyPerfect or \
            (not PlotOnlyPerfect and
             ((np.all(perfect) and PlotOnlyWorking or PlotOnlyFit) and (
                     np.all(ok) and PlotOnlyFit or not PlotOnlyFit))):

        fig, ax = plt.subplots()

        mPSD = PSD
        if PlotSuperMean:
            try:
                mPSD = np.array(PSD)
                mPSD = [np.mean(mPSD.transpose(), mPSD.ndim - 1).transpose()]
            except:
                pass

        PlotPSD(ax, Fpsd, mPSD, noisefit, perfect, ok, nType, **kwargs)


def processNoiseForPlotPSD(noise):
    noisei = np.array(noise)
    if noisei.ndim == 2:
        noisei = np.mean(noisei, 0)
    if noisei.ndim > 2:
        tn = []
        for n in noisei:
            nm = processNoiseForPlotPSD(n)
            tn.append(nm)
        tn = np.array(tn)
        noisei = np.mean(tn, 0)
    return noisei


def reduce(item):
    item = np.array(item)
    try:
        item = np.mean(item, 0)
    except ValueError:
        tl = []
        for i in item:
            t = reduce(i)
            tl.append(t)
        item = tl
    item = np.array(item)
    return item
