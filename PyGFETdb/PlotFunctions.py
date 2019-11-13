# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Plot Functions that do not fit in the previous files.
"""
import matplotlib.pyplot as plt
import numpy as np
import quantities as pq

from PyGFETdb import qty


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
    :param xlabel:
    :param ylabel:
    :param title:
    :param legendtitle:
    :param handles:
    :param kargs: Aesthetics arguments for the plot

    :return: None
    """
    plt.xticks(xPos, xLab, rotation=45, fontsize='small')
    ax.set_ylabel(ylabel=ylabel, fontsize='large')

    ax.set_xlabel(xlabel, fontsize='large')
    ax.set_title(label=title, fontsize='large')


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

    if type(units) is pq.UnitQuantity or type(units) is pq.Quantity:
        units = qty.getQuantityUnits(units)
        if units is not None:
            units = units.dimensionality.latex

    plt.xticks(xPos, xLab, rotation=45, fontsize='small')

    Ax.set_xlabel(xlabel, fontsize='large')

    param = kwargs.get('Param')

    if qty.isActive() and qtys is not None:
        qtyunits = qty.getQuantityUnits(qtys)
        if qtyunits is not None:
            qtyunits = qtyunits.dimensionality.latex
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
    Ax.set_position([chartBox.x0, chartBox.y0 + chartBox.y0 * 1.3, chartBox.width * 0.8, chartBox.height * 0.8])
    Ax.legend(title=legend, handles=handles, loc='upper right', fontsize='x-small',
              bbox_to_anchor=(0.89, 0.5, 0.5, 0.5), shadow=True, ncol=1)


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


def PlotPerTypeNoise(Results, handles=None, xlabel=None, legendTitle=None, Colors=None, perType="", **kwargs):
    """
        Plots the noise analysis

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
           Plots the percentage of transistors yield per secondary group

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
            work.append((np.array(Grc).size / totalWf) * 100)
            xLab.append(Grn)
            xPos.append(pos)
            _BoxplotValsGroup(ax2, Col, pos, work)
            pos += 1
    _closeBoxplotValsGroup(ax2, xPos, xLab, xlabel, "Yield [%{}]".format(perType), title, **kwargs)
    Legend(ax2, legendTitle, handles)


def PlotPerTypeYieldTotal(Results, title=None, Colors=None, xlabel=None,
                          **kwargs):
    """

           Plots the percentage of transistors yield per primary group

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
            work.append(((np.array(cWf)).size / total) * 100)
        _BoxplotValsGroup(ax2, Col, nt, work, **kwargs)
    _closeBoxplotValsGroup(ax2, xPos, xLab, xlabel, "Yield [% x Type]", title, **kwargs)
