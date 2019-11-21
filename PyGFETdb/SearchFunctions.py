# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Search Functions that do not fit in the previous files.

"""

import PyGFETdb.DBSearch as DbSe
from PyGFETdb import multithrds, Multiprocessing as mp, GlobalFunctions as g

# MULTIPROCESSING INITIALIZATION #################################################
if multithrds:
    search = mp.SearchDB_MP
    getparams = mp.GetParams_MP
else:
    search = mp.SearchDB
    getparams = mp.GetParams


def DBSearchPerWaferAndType(GrBase, args, **kwargs):
    """

    :param GrBase: A group of conditions
    :param args: a dict with the parameters to search
    :return: a group of conditions and the results of the search
    """
    GrWs = DbSe.GenGroups(GrBase, 'Wafers.Name', LongName=False)
    ResultsParams = {}
    for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
        print('')
        print('Searching Waffer {}...'.format(Grwn))
        GrTypes = DbSe.GenGroups(Grwc, 'TrtTypes.Name', LongName=False)
        ResultsDB = search(GrTypes, **kwargs)
        ResultsParams[Grwn] = getparams(ResultsDB, GrTypes, args)
    return GrWs, ResultsParams


def DBSearchPerWaferAndDevice(GrBase, args, **kwargs):
    """

    :param GrBase: A group of conditions
    :param args: a dict with the parameters to search
    :return: a group of conditions and the results of the search
    """
    GrWs = DbSe.GenGroups(GrBase, 'Wafers.Name', LongName=False)
    ResultsParams = {}
    for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
        print('')
        print('Searching Waffer {}...'.format(Grwn))
        GrTypes = DbSe.GenGroups(Grwc, 'Devices.Name', LongName=False)
        ResultsDB = search(GrTypes, **kwargs)
        ResultsParams[Grwn] = getparams(ResultsDB, GrTypes, args)
    return GrWs, ResultsParams


def DBSearchPerDevice(GrBase, args, **kwargs):
    """

    :param GrBase: A group of conditions
    :param args: a dict with the parameters to search
    :return: a group of conditions and the results of the search
    """
    GrWs = DbSe.GenGroups(GrBase, 'Devices.Name', LongName=False)
    ResultsParams = {}
    for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
        print('')
        print('Searching Device {}...'.format(Grwn))
        ResultsDB = search(GrWs, **kwargs)
        ResultsParams[Grwn] = getparams(ResultsDB, GrWs, args)
    return GrWs, ResultsParams


def DBSearchPerType(GrBase, args, **kwargs):
    """

    :param GrBase: a group of conditions
    :param args: a dict with the parameters to search
    :return: a group of conditions and the results of the search
    """
    GrTypes = DbSe.GenGroups(GrBase, 'TrtTypes.Name', LongName=False)
    ResultsParams = {}
    for iType, (nType, cType) in enumerate(GrTypes.items()):
        print('')
        print('Searching Type {}...'.format(nType))
        GrWfs = DbSe.GenGroups(cType, 'Wafers.Name', LongName=False)
        ResultsDB = search(GrWfs, **kwargs)
        ResultsParams[nType] = getparams(ResultsDB, GrWfs, args)
    return GrTypes, ResultsParams


def DataClassification(GrWs, arguments, ResultsParams):
    """

    :param GrWs: A group of conditions
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
                    g.updateDictOfLists(Results[Grwn], TGrn, TGrc)
        clssfResults[narg].update(Results)
    return clssfResults


def DBSearchPerWafer(GrBase, args):
    """

    :param GrBase: a group of conditions
    :param args: a dict with the parameters to search
    :return: a group of conditions and the results of the search
    """
    GrTypes = DbSe.GenGroups(GrBase, 'TrtTypes.Name', LongName=False)
    ResultsDB = search(GrTypes)
    argParams = {'ResultsDB': dict(ResultsDB), 'GrWfs': GrTypes, 'arguments': args, 'args': args}
    ResultsParams = getparams(**argParams)
    return GrTypes, ResultsParams
