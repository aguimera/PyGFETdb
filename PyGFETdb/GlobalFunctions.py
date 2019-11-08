# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.
"""
import matplotlib.pyplot as plt
import numpy as np

import PyGFETdb.DBAnalyze as DbAn
import PyGFETdb.DBSearch as DbSe
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
        dict[key] = value
    else:
        k.append(value)


def _PlotValsGroup(Ax, xLab, xPos, iGr, Grn, vals, Boxplot=False, ParamUnits=None, **kwargs):
    if vals is not None:  # and len(vals) >0:
        if Boxplot:
            Ax.boxplot(vals.transpose(), positions=(iGr + 1,))
            xPos.append(iGr + 1)
        else:
            Ax.plot(np.ones(len(vals)) * iGr, vals, '*')
            xPos.append(iGr)
            xLab.append(Grn)
    else:
        print('Empty data for: ', Grn)


def _closePlotValsGroup(Ax, xLab, xPos, qtys=None, ParamUnits=None, **kwargs):
    units = kwargs.get('Units')
    plt.xticks(xPos, xLab, rotation=45)

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


def GetParams(ResultsDB, GrWfs, args: dict, **kwargs):
    """

    :param ResultsDB: The results of a search in the database
    :param Group: A group of conditions to analyse
    :param args: Arguments for getting the parameters
    :param Plot: Bool that activates plotting
    :param kwargs:
    :return: A dict of args of dicts of groupnames and parameter found in a previous search
    """
    Results = {}
    for karg, arg in args.items():
        Results[karg] = {}
        for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
            Results[karg][Wfn] = {}
            if type(Wfc) is dict and Wfc.get('Conditions') is None:
                for iGr, (Grn, Grc) in enumerate(sorted(Wfc.items())):
                    Data = ResultsDB.get(Grn)
                    if Data is not None:
                        ParamData = DbAn.GetParam(Data, **arg)
                        if ParamData is not None:
                            Results[karg][Wfn][Grn] = ParamData
            else:
                Data = ResultsDB.get(Wfn)
                if Data is not None:
                    ParamData = DbAn.GetParam(Data, **arg)
                    if ParamData is not None:
                        Results[karg][Wfn] = ParamData
    return Results


def PlotGroup(ResultsParams, Group, args, **kwargs):
    """

    :param ResultsParams: The results of a search in the database
    :param Group: A group of conditions to analyse
    :param args: Arguments for getting the parameters
    :param kwargs:
    :return: A dict of args of dicts of groupnames and parameter found in a previous search
    """
    Results = {}
    for karg, arg in args.items():
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


def SearchDB(GrWfs, **kwargs):
    """

    :param Group: A group of conditions
    :return: The results of the search in the database
    """
    ResultsDB = {}
    for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
        ResultsDB[Wfn] = {}
        if type(Wfc) is dict and Wfc.get('Conditions') is None:
            for iGr, (Grn, Grc) in enumerate(Wfc.items()):
                Data, Trts = DbSe.GetFromDB(**Grc)
                ResultsDB[Wfn][Grn] = dict(Data)
        else:
            Data, Trts = DbSe.GetFromDB(**Wfc)
            ResultsDB[Wfn] = dict(Data)

    return dict(ResultsDB)
