# -*- coding: utf-8 -*-
"""

@author: dragc

"""

import itertools
import sys

import numpy as np

import PyGFETdb.DBAnalyze as DbAn
import PyGFETdb.DBSearch as DbSe
from PyGFETdb import multithrds, Thread, AnalysisFunctions as analysis

superthreading = False

if not superthreading:
    getParam = 'GetParam'
else:
    getParam = '_GetParam'


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
    return processResults(Results, args)


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
                kwargs.update(**Grc)
                Data, Trts = DbSe.GetFromDB(**kwargs)
                if Data is not None:
                    ResultsDB[Wfn][Grn] = dict(Data)
        else:
            kwargs.update(**Wfc)
            Data, Trts = DbSe.GetFromDB(**kwargs)
            if Data is not None:
                ResultsDB[Wfn] = dict(Data)

    return dict(ResultsDB)


def processGetFromDB(results):
    """
        **Gathers the results of different processes after a search in the database**

    :param results: Results obtained with multi-processing from the database
    :return: All the results from the search gathered in a tuple
    """
    ret = (None, None)
    if not multithrds:
        ret = results
    else:
        if type(results) is dict:
            Data = {}
            Trts = []
            try:
                for i, (data, trts) in results.items():
                    Data.update(data)
                    Trts.append(trts)
                ret = (Data, Trts)
            except ValueError:
                print(sys.exc_info())
        elif type(results) is list:
            Data = {}
            Trts = []
            try:
                for data, trts in results:
                    Data.update(data)
                    Trts.append(trts)
                ret = (Data, Trts)
            except ValueError:
                print(sys.exc_info())

    return ret


def processResults(Results, args):
    """

    :param Results: results to classify
    :param args: arguments to search
    :return: A dict of arguments with each corresponding result
    """
    Ret = {}
    if args is not None:
        for karg, arg in args.items():
            Ret[karg] = {}
            #            for rk,rd in Results.items():
            if type(Results) is dict:
                tdict = Results.get(karg)
                if tdict is not None:
                        for iWf, (Wfn, Wfc) in enumerate(tdict.items()):
                            Ret[karg][Wfn] = {}
                            if type(Wfc) is dict and Wfc.get('Conditions') is None:
                                for iGr, (Grn, Grc) in enumerate(tdict.items()):
                                    Ret[karg][Wfn][Grn] = Grc
                            else:
                                Ret[karg][Wfn] = Wfc
    else:
        Ret = Results

    return Ret


def SearchDB_MP(GrWfs, **kwargs):
    """

    :param Group: A group of conditions
    :return: The results of the search in the database
    """
    ResultsDB = {}
    thread = Thread.MultiProcess(DbSe)
    for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
        if type(Wfc) is dict and Wfc.get('Conditions') is None:
            for iGr, (Grn, Grc) in enumerate(Wfc.items()):
                key = Wfn + ' ' + Grn
                thread.initcall(key, DbSe)
                kwargs.update(**Grc)
                thread.call(key, DbSe, 'GetFromDB', {}, **kwargs)
        else:
            thread.initcall(Wfn, DbSe)
            kwargs.update(**Wfc)
            thread.call(Wfn, DbSe, 'GetFromDB', {}, **kwargs)

    for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
        ResultsDB[Wfn] = {}
        if type(Wfc) is dict and Wfc.get('Conditions') is None:
            for iGr, (Grn, Grc) in enumerate(Wfc.items()):
                ResultsDB[Grn] = {}
                key = Wfn + ' ' + Grn
                ResultsDB[Wfn][Grn] = processGetFromDB(thread.getResults(key))[0]
        else:
            ResultsDB[Wfn] = processGetFromDB(thread.getResults(Wfn))[0]

    return ResultsDB


def GetParams_MP(ResultsDB, GrWfs, arguments, **kwargs):
    """

    :param ResultsDB: The results of a search in the database
    :param Group: A group of conditions to analyse
    :param args: Arguments for getting the parameters
    :param Plot: Bool that activates plotting
    :param kwargs:
    :return: A dict of args of dicts of groupnames and parameter found in a previous search
    """
    thread = Thread.MultiProcess(DbAn)
    for karg, arg in arguments.items():
        for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
            if type(Wfc) is dict and Wfc.get('Conditions') is None:
                for iGr, (Grn, Grc) in enumerate(sorted(Wfc.items())):
                    Data = ResultsDB.get(Grn)
                    if Data is not None:
                        key = (karg + ' ' + Wfn + ' ' + Grn)
                        kargs = {'Data': Data, 'args': arg}
                        thread.initcall(key, DbAn)
                        thread.call(key, DbAn, 'GetParam', kargs, **kargs)

            else:
                Data = ResultsDB.get(Wfn)
                if Data is not None:
                    key = karg + ' ' + Wfn
                    thread.initcall(key, DbAn)
                    kargs = {'Data': Data, 'args': arg}
                    thread.call(key, DbAn, 'GetParam', kargs, **kargs)

    # """""
    res = {}
    for karg, arg in arguments.items():
        for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
            if type(Wfc) is dict and Wfc.get('Conditions') is None:
                for iGr, (Grn, Grc) in enumerate(sorted(Wfc.items())):
                    key = (karg + ' ' + Wfn + ' ' + Grn)
                    res[key] = thread.getResults(key)
            else:
                key = karg + ' ' + Wfn
                res[key] = thread.getResults(key)

    Results = processGetParams_MP(GrWfs, res, arguments)
    # """

    return Results


def processGetParams_MP(GrWfs, Results, args):
    """

    :param GrWfs: A group of conditions
    :param Results: The results of the parameter search
    :param args: a dict with the parameters to be searched
    :return: The results ordered by parameter
    """
    Ret = {}
    if args is not None:
        for karg, arg in args.items():
            Ret[karg] = {}
            if type(Results) is dict:
                processOneArg(GrWfs, Results, Ret, karg)
    else:
        Ret = Results

    return Ret


def processOneArg(GrWfs, Results, Ret, karg):
    """

    :param GrWfs: A group of conditions
    :param Results: Results of the parameter search
    :param Ret: A dict to update with Eesults
    :param karg: Key-name of the argument to process
    :return: None but modifies Ret dict
    """
    for iWf, (Wfn, Wfc) in enumerate(GrWfs.items()):
        Ret[karg][Wfn] = {}
        if type(Wfc) is dict and Wfc.get('Conditions') is None:
            for iGr, (Grn, Grc) in enumerate(Wfc.items()):
                for key, result in Results.items():
                    if key.startswith(karg + ' ' + Wfn + ' ' + Grn):
                        k = result.get(key)
                        if k is not None:
                            l = list(itertools.chain(k))
                            Ret[karg][Wfn][Grn] = l
        else:
            for key, result in Results.items():
                if key.startswith(karg + ' ' + Wfn):
                    k = result.get(key)
                    if k is not None:
                        l = list(itertools.chain(k))
                        Ret[karg][Wfn] = l


def processPSDs_MP(GrTypes, rPSD, tolerance=1.5e-22, noisetolerance=1.5e-25):
    """

    :param GrTypes: Group to process
    :param rPSD: Results of a param search of PSD, Fpsd, NoA and NoB
    :param tolerance:
    :param noisetolerance:
    :return: A dict with the results of the processing
    """
    results = {}
    keys = {}
    i = 0
    print('******************************************************************************')
    print('******* RESULTS OF THE ANALYSIS **********************************************')
    print('******************************************************************************')
    print(' ')
    thread = Thread.MultiProcess(analysis)

    for nType, vType in GrTypes.items():
        Fpsd = rPSD[nType].get('Fpsd')
        PSD = rPSD[nType]['PSD']
        NoA = rPSD[nType]['NoA']
        NoB = rPSD[nType]['NoB']
        keys[nType] = {}
        for nWf, vWf in Fpsd.items():
            PSDt = PSD[nWf]
            NoAt = NoA[nWf]
            NoBt = NoB[nWf]

            Fpsdt = vWf[0:len(PSDt)]
            Fpsd2t = np.array(Fpsdt).reshape((1, len(PSDt)))
            args = {
                'PSD': PSDt,
                'Fpsd': Fpsd2t,
                'NoA': NoAt,
                'NoB': NoBt,
                'tolerance': tolerance,
                'noisetolerance': noisetolerance
            }
            key = thread.initcall(Thread.key(), analysis)
            keys[nType][nWf] = key
            thread.call(key, analysis, 'processNoise', args, **args)

    for nType, vType in GrTypes.items():
        Fpsd = rPSD[nType].get('Fpsd')
        PSD = rPSD[nType]['PSD']
        results[nType] = {}
        for nWf, vWf in Fpsd.items():
            PSDt = PSD[nWf]
            Fpsdt = vWf[0:len(PSDt)]
            Fpsd2t = np.array(Fpsdt).reshape((1, len(PSDt)))
            key = keys[nType][nWf]
            r = thread.getResults(key)
            [noise, ok, perfect, grad, gradnoise] = r[key]
            results[nType][nWf] = [Fpsdt, PSDt, Fpsd2t, noise, ok, perfect, grad, gradnoise]

    return results
