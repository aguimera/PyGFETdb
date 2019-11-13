# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.
"""

import PyGFETdb.DBSearch as DbSe
from PyGFETdb import multithrds, Multiprocessing as mp

# MULTIPROCESSING INITIALIZATION #################################################
if multithrds:
    search = mp.SearchDB_MP
    getparams = mp.GetParams_MP
else:
    search = mp.SearchDB
    getparams = mp.GetParams

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


def DBSearchPerWaferAndType(GrBase, args):
    """

    :param GrBase:
    :param args:
    :return: a group of conditions and the results of the search
    """
    GrWs = DbSe.GenGroups(GrBase, 'Wafers.Name', LongName=False)
    ResultsParams = {}
    for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
        GrTypes = DbSe.GenGroups(Grwc, 'TrtTypes.Name', LongName=False)
        ResultsDB = search(GrTypes)
        ResultsParams[Grwn] = getparams(ResultsDB, GrTypes, args)
    return GrWs, ResultsParams


def DBSearchPerType(GrBase, args):
    """

    :param GrBase:
    :param args:
    :return: a group of conditions and the results of the search
    """
    GrTypes = DbSe.GenGroups(GrBase, 'TrtTypes.Name', LongName=False)
    ResultsParams = {}
    for iType, (nType, cType) in enumerate(GrTypes.items()):
        GrWfs = DbSe.GenGroups(cType, 'Wafers.Name', LongName=False)
        ResultsDB = search(GrWfs)
        ResultsParams[nType] = getparams(ResultsDB, GrWfs, args)
    return GrTypes, ResultsParams


def DataClassification(GrWs, arguments, ResultsParams):
    """

    :param GrWs:
    :param ResultsParams: the results obtained of a search of parameters
    :return: A dict with the results classified by argument to plot
    """
    clssfResults = {}
    Results = {}
    for iarg, (narg, carg) in enumerate(arguments.items()):
        clssfResults[narg] = {}
        for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
            Results[Grwn] = {}
            k1 = ResultsParams.get(Grwn)
            k2 = k1.get(narg)
            for iType, (TGrn, TGrc) in enumerate(k2.items()):
                if TGrc is not None:
                    k3 = Results[Grwn].get(TGrn)
                    if k3 is None:
                        Results[Grwn][TGrn] = []
                    updateDictOfLists(Results[Grwn], TGrn, TGrc)
        clssfResults[narg].update(Results)
    return clssfResults

