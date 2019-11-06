# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.
"""
import matplotlib.pyplot as plt
import numpy as np

import PyGFETdb
import PyGFETdb.DBAnalyze as DbAn
import PyGFETdb.DBSearch as DbSe
import PyGFETdb.Thread as Thread
from PyGFETdb import multithrds
from PyGFETdb import qty

"""
pool = Thread.PyFETdb(DbSe)
for iGr, (Grn, Grc) in enumerate(sorted(GrDevs.items())):
    print('Getting data for ', Grn)
    #Data, Trts = DbSe.GetFromDB(**Grc)
    pool.call('GetFromDB',[Grc])
ResultsDB = pool.getResults()
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


def GetParamsThread(args, ResultsDB, Multidata=None, **kwargs):
    # GetParamsThread
    pool = Thread.PyFETdb(DbAn)
    if Multidata is not None:
        for r in ResultsDB:
            for a in args:
                _GetParamsThread(pool, a, r)
    else:
        _GetParamsThread(pool, args, ResultsDB)
    listres = pool.getResults()
    del pool
    return listres


def _GetParamsThread(pool, args, ResultsDB):
    for Grn, (Data, Trts) in ResultsDB:  #
        if len(Data) > 0:
            for i, arg in enumerate(args):
                if multithrds and pool:
                    pool.call('GetParamThread', [Grn, Data, i], **arg)
                else:
                    return DbAn.GetParamThread(Grn, Data, i, **arg)


def processResults(args, ResultsDB, MultiData=None, **kwargs):
    if MultiData is not None:
        Results = []
        for r in ResultsDB:
            for a in args:
                Results.append(_processResults(r, a, **kwargs))
        return Results
    else:
        return _processResults(ResultsDB, args, **kwargs)


def _processResults(Grc, args, groupcallback=None, argcallback=None, **kwargs):
    args = np.array(args)
    for i, arg in enumerate(args):
        for res, results in Grc:
            if i == results[0]:  # argument number
                for grn, [Param, Data] in results[1].items():  # group dict
                    groupcallback(grn, Param, Data, **kwargs)
        argcallback(grn, Param, Data, **kwargs)
    return results[1]


def Plot(GrDevs, Grc, args, Boxplot=False, ParamUnits=None, **kwargs):
    if multithrds:
        return _PlotThreads(GrDevs, Grc, args, **kwargs)
    else:
        _Plot(GrDevs, Grc, **args)


def updateDictOfLists(dict, key, value):
    k = dict.get(key)
    if k is None:
        dict[key] = value
    else:
        k.append(value)


def _Plot(Groups, Grc, Boxplot=False, ParamUnits=None, **kwargs):
    Vals = {}
    fig, Ax = plt.subplots()
    for iGr, (Grn) in enumerate(Groups.items()):
        xLab = []
        xPos = []
        for (Grn, data) in Grc.items():
            qtys = data
            vals = np.array(data)
            if vals is None:
                continue
            if vals.size == 0:
                continue
            vals = qty.flatten(vals)
            # vals = np.array(vals)
            PlotValsGroup(Ax, xLab, xPos, iGr, Grn, vals, **kwargs)
            updateDictOfLists(Vals, Grn, qtys)
    if vals is not None and len(vals) > 0:
        closePlotValsGroup(Ax, xLab, xPos, qtys, **kwargs)
    return Vals


def _PlotThreads(GrDevs, Grc, args, Boxplot=False, ParamUnits=None, **kwargs):
    args = np.array(args)
    Vals = {}
    for iarg, arg in enumerate(args):
        fig, Ax = plt.subplots()
        qtys = None
        param = arg.get('Param')
        xLab = []
        xPos = []
        iGr = 0
        units = None
        for gr in GrDevs.items():
            tdata = []  # np.array()
            tqtys = []
            for results in Grc:
                for iarg2, a2 in results.items():
                    if iarg2 == iarg:
                        if a2 is not None:
                            p = a2.get(param)
                            if p is not None:
                                for grn, data in p.items():
                                    # tdata.append(data)
                                    qtys = data
                                    vals = np.array(data)
                                    if vals is None:
                                        continue
                                    if vals.size == 0:
                                        continue
                                    tdata.append(vals)
                                    tqtys.append(qtys)
            tdata = qty.flatten(tdata)
            tdata = np.array(tdata)
            if tdata is not None and len(tdata) > 0:
                PlotValsGroup(Ax, xLab, xPos, iGr, gr, tdata, **args[iarg])
                iGr += 1
                Vals[grn] = tqtys
        if tdata is not None and len(tdata) > 0:
            closePlotValsGroup(Ax, xLab, xPos, units, **args[iarg])
            plt.show()
    return Vals


def PlotValsGroup(Ax, xLab, xPos, iGr, Grn, vals, Boxplot=False, ParamUnits=None, **kwargs):
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


def closePlotValsGroup(Ax, xLab, xPos, qtys=None, ParamUnits=None, **kwargs):
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
