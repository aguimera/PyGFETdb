# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.

"""

import numpy as np


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
                    updateDictOfLists(Results[Grwn], TGrn, TGrc)
        clssfResults[narg].update(Results)
    return clssfResults


def remove(Valx: np.array, index):
    Valx: list = Valx.tolist()
    Valx.remove(Valx[index])
    Valx = np.array(Valx)
    return Valx


def process50Hz(Array, process: bool):
    if process:
        for i in range(1, 2):
            Array = remove(Array, 48)
    return Array
