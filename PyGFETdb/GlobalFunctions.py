# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.
"""

import quantities as pq

import PyGFETdb.DBSearch as DbSe
from PyGFETdb import qty, multithrds, Multiprocessing as mp

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


def DBSearchPerWaferAndType(GrBase, args):
    GrWs = DbSe.GenGroups(GrBase, 'Wafers.Name', LongName=False)
    ResultsParams = {}
    for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
        GrTypes = DbSe.GenGroups(Grwc, 'TrtTypes.Name', LongName=False)
        ResultsDB = search(GrTypes)
        ResultsParams[Grwn] = getparams(ResultsDB, GrTypes, args)
    return GrWs, ResultsParams


def DBSearchPerType(GrBase, args):
    GrTypes = DbSe.GenGroups(GrBase, 'TrtTypes.Name', LongName=False)
    ResultsParams = {}
    for iType, (nType, cType) in enumerate(GrTypes.items()):
        GrWfs = DbSe.GenGroups(cType, 'Wafers.Name', LongName=False)
        ResultsDB = search(GrWfs)
        ResultsParams[nType] = getparams(ResultsDB, GrWfs, args)
    return GrTypes, ResultsParams


def DataClassification(GrWs, ResultsParams):
    clssfResults = {}
    Results = {}
    for iWf, (Grwn, Grwc) in enumerate(GrWs.items()):
        for iarg, (narg, arg) in enumerate(sorted(ResultsParams.get(Grwn).items())):  # Params
            for iType, (TGrn, TGrc) in enumerate(arg.items()):
                quantities = TGrc  # Param 0
                if qty.isActive():
                    Vals = pq.Quantity(quantities)  # np.array also works fine
                else:
                    Vals = quantities * 1e6
                if Vals is not None:
                    updateDictOfDicts(Results, Grwn, TGrn, Vals)
            clssfResults[narg] = Results
    return clssfResults
