# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Analysis Functions that do not fit in the previous files.

"""

import numpy as np

from PyGFETdb import qty, GlobalFunctions as g
from PyGFETdb.NoiseModel import Fnoise


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
            Array = g.remove(Array, 48)
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
            Array = g.remove(Array, 0)
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
        for i in range(1, 27):  # To widen the effect increase the 25
            Array = g.remove(Array, Array.size - 1)
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


def processNoise(PSD, Fpsd, NoA, NoB, tolerance=1.5e-22, errortolerance=1.3e-19, gradtolerance=-2.7e-18):
    """

    :param PSD: PSD of a Group
    :param Fpsd: Fpsd of a Group
    :param NoA: NoA of a Group
    :param NoB: NoB of a Group
    :param tolerance: margin of error allowed on the analysis of the Mean PSD
    :return: noise:
    :return: ok: True if the noise is fitted well
    :return: perfect: True if the Mean PSD gradient is acceptable
    :return: grad: Mean PSD gradient
    :return: noisegrad: Noise gradient

    """
    Fpsd2 = Fpsd
    noise = None
    if NoA is not None and len(NoA) > 0 and len(NoA[0].shape) == 1:  # Only a wafer
        NoA = np.array(NoA)
        NoB = np.array(NoB)

        NoA = np.mean(NoA.transpose(), 1)
        NoA = NoA.reshape((1, NoA.size))

        NoB = np.mean(NoB.transpose(), 1)
        NoB = NoB.reshape((1, NoB.size))
        noise = Fnoise(Fpsd2, NoA[:, len(Fpsd2)], NoB[:, len(Fpsd2)])
    elif NoA is not None:  # more than one wafer
        rA = []
        for item in NoA:
            rA.append(np.mean(item))

        rB = []
        for item in NoB:
            rB.append(np.mean(item))

        NoA = np.array(rA)
        NoB = np.array(rB)

        noise = Fnoise(Fpsd2, NoA[:len(Fpsd2)], NoB[:len(Fpsd2)])

    ok, perfect, grad, noisegrad = isMeanPSDok(PSD, Fpsd2, noise, tolerance, errortolerance, gradtolerance)

    return [noise, ok, perfect, grad, noisegrad]


def processPSDs(GrTypes, rPSD, tolerance=1.5e-22, errortolerance=1.3e-19, gradtolerance=-2.7e-18):
    """

    :param GrTypes: Group to process
    :param rPSD: Results of a param search of PSD, Fpsd, NoA and NoB
    :param tolerance:
    :return: A dict with the results of the processing
    """
    results = {}
    i = 0
    okc = 0
    perfectc = 0
    print(' ')
    print('******************************************************************************')
    print('******* RESULTS OF THE NOISE ANALYSIS ****************************************')
    print('******************************************************************************')
    print(' ')
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
            if NoAt is not None and len(NoAt) > 0:
                Fpsdt = vWf[0:len(PSDt)]
                Fpsd2t = np.array(Fpsdt).reshape((1, len(PSDt)))
                i += 1
                print('***************************************************')
                print('{}) Group:{}, Subgroup:{}'.format(i, nType, nWf))
                print('***************************************************')

                [noise, ok, perfect, grad, noisegrad] = processNoise(PSDt, Fpsd2t, NoAt, NoBt,
                                                                 tolerance, errortolerance, gradtolerance)
                print(' ')
                if ok:
                    okc += 1
                if perfect:
                    perfectc += 1
                results[nType][nWf] = [Fpsdt, PSDt, Fpsd2t, noise, ok, perfect, grad, noisegrad]

    print('******************************************************************************')
    print('*******  NOISE ANALYSIS SUMMARY ********************************************')
    print('******************************************************************************')
    print('Perfect PSDs -> {} of {} : {} %'.format(perfectc, i, (perfectc / i) * 100))
    print('Noise Fitted OK -> {} of {} : {} %'.format(okc, i, (okc / i) * 100))
    print('******************************************************************************')
    print('******* END OF THE NOISE ANALYSIS ********************************************')
    print('******************************************************************************')
    return results


def isMeanPSDok(PSD, Fpsd, noise, tolerance=2.5e-2, errortolerance=0.8, gradtolerance=0.08):
    """

       :param PSD: PSD of a Group
       :param Fpsd: Fpsd of a Group
       :param noise: Mean noise
       :param tolerance: margin of error allowed on the analysis of the Mean PSD
       :param tolerance: margin of error allowed on the analysis of the Mean Noise
       :return: ok: True if the noise is fitted well
       :return: perfect: True if the Mean PSD gradient is acceptable
       :return: grad: Mean PSD gradient
       :return: noisegrad: Noise gradient

       """
    mPSD = np.mean(PSD, 1)

    dx = np.diff(Fpsd)

    y = np.diff(mPSD)
    y2 = np.diff(noise)

    grad = qty.Divide(y, dx) / np.max(mPSD)  # Gradient of the mean PSD
    perfect = np.all(grad <= tolerance)

    noisegrad = qty.Divide(y2, dx) / np.max(mPSD)  # Gradient of the noise fitting

    graderror = grad - noisegrad
    error = (mPSD - noise) / np.max(mPSD)
    minerr = np.min(error)
    meangraderr = np.mean(graderror)

    # ok1 = minerr > 0 and minerr < errortolerance
    # ok2 = meangraderr < 0 and np.abs(meangraderr) < gradtolerance

    ok1 = minerr < errortolerance
    ok2 = np.abs(meangraderr) < gradtolerance

    ok = ok1 and ok2

    if perfect:
        print('Mean PSD Noise PERFECT -> {}'.format(np.max(grad)))
    else:
        print('Mean PSD Noise BAD -> {}'.format(np.max(grad)))

    if ok:
        print('Noise Fitted OK -> error:{} grad-error:{}'.format(minerr, meangraderr))
    else:
        print('Noise Fitted BAD -> error:{} grad-error:{}'.format(minerr, meangraderr))
        if not ok1:
            print('     Noise Error BAD -> error:{}'.format(minerr))
        if not ok2:
            print('     Noise GradError BAD -> grad-error:{}'.format(meangraderr))


    return ok, perfect, grad, noisegrad


def isPSDok(PSD, Fpsd, noise, tolerance=0.5, errortolerance=-1.3, gradtolerance=0.05):
    """

       :param PSD: PSD of a Group
       :param Fpsd: Fpsd of a Group
       :param noise: Mean noise
       :param tolerance: margin of error allowed on the analysis of the Mean PSD
       :param tolerance: margin of error allowed on the analysis of the Mean Noise
       :return: ok: True if the noise is fitted well
       :return: perfect: True if the Mean PSD gradient is acceptable
       :return: grad: Mean PSD gradient
       :return: noisegrad: Noise gradient

       """
    mPSD = PSD[1:]

    dx = np.diff(Fpsd.transpose())[1:]

    y = np.diff(mPSD)
    y2 = np.diff(noise.transpose())

    grad = qty.Divide(y, dx) / np.max(mPSD)  # Gradient of the mean PSD
    perfect = np.all(grad <= tolerance)

    noisegrad = qty.Divide(y2, dx) / np.max(mPSD)  # Gradient of the noise fitting

    graderror = grad - noisegrad
    error = (mPSD - noise) / np.max(mPSD)
    minerr = np.min(error)
    meangraderr = np.mean(graderror)

    ok1 = minerr > errortolerance
    ok2 = np.abs(meangraderr) <= gradtolerance

    ok = ok1 and ok2




    if perfect:
        print('PSD Noise OK -> {}'.format(np.max(grad)))
    else:
        print('PSD Noise BAD -> {}'.format(np.max(grad)))

    if ok:
        print('Noise Fitted OK -> error:{} grad-error:{}'.format(minerr, meangraderr))
    else:
        print('Noise Fitted BAD -> error:{} grad-error:{}'.format(minerr, meangraderr))
        if not ok1:
            print('     Noise Error BAD -> error:{}'.format(minerr))
        if not ok2:
            print('     Noise GradError BAD -> grad-error:{}'.format(meangraderr))

    return ok, perfect, grad, noisegrad
