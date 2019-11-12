# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.
"""
import matplotlib.pyplot as plt
import numpy as np

from PyGFETdb import qty


def updateDictOfLists(dict, key, value):
    """
    Modifies a dictionary of lists, appending the value at the list obtained
    of applying the key to the dictionary
    :param dict: A dictionary to update
    :param key: The key to update
    :param value: The value to append
    :return: None
    """
    k = dict.get(key)
    if k is None:
        dict[key] = [value]
    else:
        k.append(value)


def updateDictOfDicts(dict, key1, key2, value):
    """
    Modifies a dictionary of dictionaries, updating the value at the dictionary obtained
    of applying the key to the dictionary
    :param dict: A dictionary to update
    :param key1: The key to search in the dictionary
    :param key2: The key to update in the result of searching the first key
    :param value: The value to update
    :return: None
    """
    k = dict.get(key1)
    if k is None:
        dict[key1] = {key2: value}
    else:
        k[key2] = value

def _BoxplotValsGroup(ax, Col, iGroup, Vals, **kwargs):
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
    plt.xticks(xPos, xLab, rotation=45, fontsize='small')
    ax.set_ylabel(ylabel=ylabel, fontsize='large')

    ax.set_xlabel(xlabel, fontsize='large')
    ax.set_title(label=title, fontsize='large')


def _PlotValsGroup(Ax, xLab, xPos, iGr, Grn, vals,
                   Boxplot=False, ParamUnits=None, **kwargs):

    if vals is not None:  # and len(vals) >0:
        if Boxplot:
            Ax.boxplot(vals.transpose(), positions=(iGr + 1,))
            xPos.append(iGr + 1)
        else:
            Ax.plot(np.ones(len(vals)) * iGr, vals, '.', label=(Grn,))
            xPos.append(iGr)
            xLab.append(Grn)
    else:
        print('Empty data for: ', Grn)


def _closePlotValsGroup(Ax, xLab, xPos, qtys=None, ParamUnits=None,
                        title=None, legendTitle=None, handles=None, xlabel='', **kwargs):
    units = kwargs.get('Units')
    plt.xticks(xPos, xLab, rotation=45, fontsize='small')

    Ax.set_xlabel(xlabel, fontsize='large')

    if ParamUnits is not None:
        Ax.set_ylabel(kwargs['Param'] + '[' + ParamUnits + ']')
    else:
        Ax.set_ylabel(kwargs['Param'])

    if qty.isActive() and qtys is not None:
        qtyunits = qty.getQuantityUnits(qtys)
        if qtyunits:
            Ax.set_ylabel(kwargs['Param'] + '[' + qtyunits + ']')
        elif units is not None:
            Ax.set_ylabel(kwargs['Param'] + '[' + units + ']')

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
    chartBox = Ax.get_position()
    Ax.set_position([chartBox.x0, chartBox.y0 + chartBox.y0 * 1.3, chartBox.width * 0.8, chartBox.height * 0.8])
    Ax.legend(title=legend, handles=handles, markerscale=10, loc='upper right',
              bbox_to_anchor=(0.89, 0.5, 0.5, 0.5), shadow=True, ncol=1)


def PlotGroup(ResultsParams, Group, args, **kwargs):
    """

    :param ResultsParams: The results of a search in the database
    :param Group: A group of conditions to analyse
    :param args: Arguments for getting the parameters
    :param kwargs:
    :return: A dict of args of dicts of groupnames and parameter found in a previous search
    """
    Results = {}
    for iarg, (karg, arg) in enumerate(args.items()):
        Results[karg] = {}
        fig, Ax = plt.subplots()
        xLab = []
        xPos = []
        qtys = None
        for iGr, (Grn, Grc) in enumerate(sorted(Group.items())):
            argRes = ResultsParams.get(karg)
            if argRes is not None:
                ParamData = argRes.get(Grn)
                if ParamData is not None:
                    Results[karg][Grn] = ParamData
                    if qty.isActive():
                        qtys = np.array(ParamData)
                        ParamData = qty.flatten(ParamData)
                    ParamData = np.array(ParamData)
                    _PlotValsGroup(Ax, xLab, xPos, iGr, Grn, ParamData, **arg)
        _closePlotValsGroup(Ax, xLab, xPos, qtys, **arg)
    return Results


