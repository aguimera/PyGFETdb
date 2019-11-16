# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Search Functions that do not fit in the previous files.

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


def DBSearchPerWaferAndType(GrBase, args, **kwargs):
    """

    :param GrBase: A group of conditions
    :param args: a dict with the parameters to search
    :return: a group of conditions and the results of the search
    """
    GrWs = DbSe.GenGroups(GrBase, 'Wafers.Name', LongName=False)
    ResultsParams = {}
    for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
        GrTypes = DbSe.GenGroups(Grwc, 'TrtTypes.Name', LongName=False)
        ResultsDB = search(GrTypes, **kwargs)
        ResultsParams[Grwn] = getparams(ResultsDB, GrTypes, args)
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
        GrWfs = DbSe.GenGroups(cType, 'Wafers.Name', LongName=False)
        ResultsDB = search(GrWfs, **kwargs)
        ResultsParams[nType] = getparams(ResultsDB, GrWfs, args)
    return GrTypes, ResultsParams
