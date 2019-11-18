# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.

"""

import numpy as np


def updateDictOfLists(dict, key, value):
    """
        **Modifies a dictionary of lists, appending the value at the list obtained
        of applying the key to the dictionary**

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


def remove(Valx, index):
    """
        **Removes the value at the index specified from an array**

    :param Valx: array to remove the value
    :param index: index to remove the value
    :return: the array without the value at the index specified
    """
    if type(Valx) is list:
        Valx.remove(Valx[index])
    else:
        Valx = Valx.tolist()
        Valx.remove(Valx[index])
        Valx = np.array(Valx)
    return Valx


def process50Hz(Array, process):
    """
        **Removes the frequency 50Hz**

    :param Array: array of frequencies
    :param process: bool that activates the function
    :return: the array without the frequency 50Hz
    """
    if process:
        # remove 50Hz
        for i in range(1, 2):  # To widen the effect increase the 2
            Array = remove(Array, 48)
    return Array


def processBelow1Hz(Array, process):
    """
        **Removes the frequencies below 1Hz**

    :param Array: array of frequencies
    :param process: bool that activates the function
    :return: the array without the frequencies below 1Hz
    """
    if process:
        #  remove below 1Hz
        for i in range(1, 15):  # To widen the effect increase the 15
            Array = remove(Array, 0)
    return Array


def processHiFreqs(Array, process):
    """
        **Removes the higher frequencies**

    :param Array: array of frequencies
    :param process: bool that activates the function
    :return: the array without the higher frequencies
    """
    if process:
        # remove the highest frequencies
        for i in range(1, 22):  # To widen the effect increase the 22
            Array = remove(Array, Array.size - 1)
    return Array


def processFreqs(Array, process):
    """
            **Removes all the unwanted frequencies**

        :param Array: array of frequencies
        :param process: bool that activates the function
        :return: the array without the unwanted frequencies
        """
    Array = process50Hz(Array, process)
    Array = processBelow1Hz(Array, process)
    Array = processHiFreqs(Array, process)
    return Array
