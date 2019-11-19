# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.

"""

import numpy as np

from PyGFETdb import qty
from PyGFETdb.NoiseModel import Fnoise


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
        for i in range(1, 25):  # To widen the effect increase the 22
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


def processNoise(PSD, Fpsd, NoA, NoB):
    Fpsd2 = Fpsd

    if NoA is not None and len(NoA[0].shape) == 1:  # Only a wafer
        NoA = np.array(NoA)
        NoB = np.array(NoB)

        NoA = np.mean(NoA.transpose(), 1)
        NoA = NoA.reshape((1, NoA.size))

        NoB = np.mean(NoB.transpose(), 1)
        NoB = NoB.reshape((1, NoB.size))
        noise = Fnoise(Fpsd2, NoA[:, len(PSD)], NoB[:, len(PSD)])
    else:  # more than one wafer
        rA = []
        for item in NoA:
            rA.append(np.mean(item))

        rB = []
        for item in NoB:
            rB.append(np.mean(item))

        NoA = np.array(rA)
        NoB = np.array(rB)

        noise = Fnoise(Fpsd2, NoA[:len(Fpsd2)], NoB[:len(Fpsd2)])

    ok, perfect, grad, noisegrad = isPSDok(PSD, Fpsd2, noise)

    return noise, ok, perfect, grad, noisegrad


def processPSDs(GrTypes, rPSD):
    results = {}
    for nType, vType in GrTypes.items():
        Fpsd = rPSD[nType].get('Fpsd')
        results[nType] = {}
        PSD = rPSD[nType]['PSD']
        NoA = rPSD[nType]['NoA']
        NoB = rPSD[nType]['NoB']
        for nWf, vWf in Fpsd.items():
            PSDt = PSD[nWf]
            NoAt = NoA[nWf]
            NoBt = NoB[nWf]

            Fpsdt = vWf[0:len(PSDt)]
            Fpsd2t = np.array(Fpsdt).reshape((1, len(PSDt)))

            noise, ok, perfect, grad, noisegrad = processNoise(PSDt, Fpsd2t, NoAt, NoBt)
            results[nType][nWf] = (Fpsdt, PSDt, Fpsd2t, noise, ok, perfect, grad, noisegrad)
    return results


def isPSDok(PSD, Fpsd, noise):
    mPSD = np.mean(PSD, 1)

    dx = np.diff(Fpsd)
    y = np.diff(mPSD)

    grad = qty.Divide(y, dx)
    perfect = np.all(grad < 0)  #

    y2 = np.diff(noise)
    noisegrad = qty.Divide(y2, dx)

    ok = np.all(noisegrad < 0)

    if perfect:
        print('PSD Noise PERFECT -> {}'.format(np.max(grad)))
    else:
        print('PSD Noise NOK -> {}'.format(np.max(grad)))

    if ok:
        print('Noise Fitted OK -> {}'.format(np.mean(noisegrad)))

    return ok, perfect, grad, noisegrad
