# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.
"""
import matplotlib.pyplot as plt
import numpy as np

import PyGFETdb.DBAnalyze as DbAn
import PyGFETdb.DBSearch as DbSe
import PyGFETdb.Thread as Thread
from PyGFETdb import multithrds
from PyGFETdb import qty

"""
def getFromDB(Groups):
    PyGFETdb.multithrds = False
    pool = Thread.PyFETdb(DbSe)
    ret = list([])
    for iGr, (Grn, Grc) in enumerate(sorted(Groups.items())):
        print('Getting data for ', Grn)
        if PyGFETdb.multithrds:
            pool.call('GetFromDB', [Grc])
            ret = pool.getResults()
            del pool
        else:
            ret.append([Grn, DbSe.GetFromDB(**Grc)])
    PyGFETdb.multithrds = True
    return ret
"""


def GetParamsThread(args, ResultsDB, Group, **kwargs):
    # GetParamsThread
    if multithrds:  # is not None:
        pool = Thread.PyFETdb(DbAn)
        for iarg, arg in enumerate(args):
            for iGr, (Grn, Grc) in enumerate(sorted(Group.items())):
                pool.call('GetParamThread', [Grn, ResultsDB[Grn], iarg], **arg)
        tlist = pool.getResults()
        listres = processResults(args, tlist, **kwargs)
        del pool
    else:
        listres = GetParams(ResultsDB, **kwargs)
    return listres


def processResults(args, ResultsDict, **kwargs):
    Results = {}
    for iarg, (arg) in enumerate(args):
        Results[iarg] = {}
        for r, rd in ResultsDict.items():
            dict = rd.get(iarg)
            if dict is not None:
                for kn, kv in dict.items():
                    for iGr, (Grn, Grc) in enumerate(kv.items()):
                        Results[iarg][Grn] = Grc
        return Results


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

    if qty.Quantities and qtys is not None:
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


def GetParams(ResultsDB, Group, args, **kwargs):
    """

    :param ResultsDB: The results of a search in the database
    :param Group: A group of conditions to analyse
    :param args: Arguments for getting the parameters
    :param Plot: Bool that activates plotting
    :param kwargs:
    :return: A dict of args of dicts of groupnames and parameter found in a previous search
    """
    Results = {}
    for iarg, arg in enumerate(args):
        Results[iarg] = {}
        for iGr, (Grn, Grc) in enumerate(sorted(Group.items())):
            Data = ResultsDB[Grn]
            ParamData = DbAn.GetParam(Data, **arg)
            if ParamData is not None:
                Results[iarg][Grn] = ParamData
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
    for iarg, arg in enumerate(args):
        Results[iarg] = {}
        fig, Ax = plt.subplots()
        xLab = []
        xPos = []
        for iGr, (Grn, Grc) in enumerate(sorted(Group.items())):
            ParamData = ResultsParams.get(iarg).get(Grn)
            if ParamData is not None:
                Results[iarg][Grn] = ParamData
                if qty.isActive():
                    qtys = np.array(ParamData)
                    ParamData = qty.flatten(ParamData)
                ParamData = np.array(ParamData)
                _PlotValsGroup(Ax, xLab, xPos, iGr, Grn, ParamData, **arg)
        _closePlotValsGroup(Ax, xLab, xPos, qtys, **arg)
    return Results


def SearchDB(Group):
    """

    :param Group: A group of conditions
    :return: The results of the search in the database
    """
    ResultsDB = {}
    for iGr, (Grn, Grc) in enumerate(sorted(Group.items())):
        Data, Trts = DbSe.GetFromDB(**Grc)
        ResultsDB[Grn] = dict(Data)
    return ResultsDB
