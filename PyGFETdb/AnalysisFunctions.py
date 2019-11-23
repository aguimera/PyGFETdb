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
        for i in range(1, 20):  # To widen the effect increase the 15
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


def processNoise(PSD, Fpsd, NoA, NoB, tolerance=1, errortolerance=1, gradtolerance=1):
    """

    :param PSD: PSD of a Group
    :param Fpsd: Fpsd of a Group
    :param NoA: NoA of a Group
    :param NoB: NoB of a Group
    :param tolerance: margin of error allowed on the analysis of the Mean PSD
    :param errortolerance: Maximum error allowed in the Noise fit
    :param gradtolerance: Maximum error of the gradient in the Noise fit

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


def isMeanPSDok(PSD, Fpsd, noise, tolerance=3.1e-2, errortolerance=-0.1, gradtolerance=0.09):
    """

        :param PSD: PSD of a Group
        :param Fpsd: Fpsd of a Group
        :param noise: Mean noise
        :param tolerance: margin of error allowed on the analysis of the Mean PSD
        :param errortolerance: Maximum error allowed in the Noise fit
        :param gradtolerance: Maximum error of the gradient in the Noise fit

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

    ok1 = minerr <= errortolerance
    ok2 = np.abs(meangraderr) <= gradtolerance

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


def processAllNoise(PSD, Fpsd, NoA, NoB, tolerance=1.0, errortolerance=1.0, gradtolerance=1.0):
    """

    :param PSD: PSD of a Group
    :param Fpsd: Fpsd of a Group
    :param NoA: Noise A of a Group
    :param NoB: Noise B of a Group
    :param tolerance: margin of error allowed on the analysis of the  PSD
    :param errortolerance: Maximum error allowed in the Noise fit
    :param gradtolerance: Maximum error of the gradient in the Noise fit
    :return: noise:
    :return: ok: True if the noise is fitted well
    :return: perfect: True if the PSD gradient is acceptable
    :return: grad: PSD gradient
    :return: noisegrad: Noise gradient

    """
    noise = np.array([])
    ok = False
    perfect = False
    grad = 0
    noisegrad = np.array([])
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    if NoA is not None and len(NoA) > 0:
        f = np.array(Fpsd).reshape((1, len(PSD)))
        noise = Fnoise(f, NoA[:, :len(f)], NoB[:, :len(f)])
        mPSD = np.mean(PSD, 1)
        for item in noise:
            [ok, perfect, grad, noisegrad] = isPSDok(mPSD, Fpsd, item.transpose()[1:],
                                                     tolerance, errortolerance, gradtolerance)
            temp1.append(ok)
            temp2.append(perfect)
            temp3.append(grad)
            temp4.append(noisegrad)

        noise = np.mean(noise.transpose(), 1)
        noise = noise.reshape(len(noise))

        ok = np.all(temp1)
        perfect = np.all(temp2)

        temp3 = np.array(temp3)
        grad = np.mean(temp3.transpose(), 1)
        grad = grad.reshape(len(grad))

        temp4 = np.array(temp4)
        noisegrad = np.mean(temp4.transpose(), 1)
        noisegrad = noisegrad.reshape(len(noisegrad))

    return [noise, ok, perfect, grad, noisegrad]


def processAllPSDs(GrTypes, rPSD, tolerance=1, errortolerance=1, gradtolerance=1):
    """

    :param GrTypes: Group to process
    :param rPSD: Results of a param search of PSD, Fpsd, NoA and NoB
    :param tolerance: Maximum error allowed in the PSD
    :param errortolerance: Maximum error allowed in the Noise fit
    :param gradtolerance: Maximum error of the gradient in the Noise fit
    :return: A dict with the results of the processing
    """
    results = {}
    i = 0
    itype = []
    itypeo = []
    iwf = []
    iwfo = []
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
        perfectct = 0
        okct = 0
        it = 0
        for nWf, vWf in Fpsd.items():
            PSDt = PSD[nWf]
            NoAt = NoA[nWf]
            NoBt = NoB[nWf]
            perfectcw = 0
            okcw = 0
            iw = 0
            if NoAt is not None and len(NoAt) > 0:
                Fpsdt = vWf[0:len(PSDt)]
                i += 1
                it += 1
                iw += 1
                print('***************************************************')
                print('{}) Group:{}, Subgroup:{}'.format(i, nType, nWf))
                print('***************************************************')

                if type(NoAt) is list:
                    temp0 = []
                    temp1 = []
                    temp2 = []
                    temp3 = []
                    temp4 = []
                    for i, item in enumerate(NoAt):
                        [noise, ok, perfect, grad, noisegrad] = processAllNoise(PSDt, Fpsdt, NoAt[i], NoBt[i],
                                                                                tolerance, errortolerance,
                                                                                gradtolerance)
                        temp0.append(noise)
                        temp1.append(ok)
                        temp2.append(perfect)
                        temp3.append(grad)
                        temp4.append(noisegrad)

                    temp0 = np.array(temp0)
                    temp0 = temp0.transpose()
                    noise = np.mean(temp0, 1)
                    ok = np.all(temp1)
                    perfect = np.all(temp2)
                    grad = np.mean(temp3)
                    noisegrad = np.mean(temp4, 1)
                else:
                    [noise, ok, perfect, grad, noisegrad] = processAllNoise(PSDt, Fpsdt, NoAt, NoBt,
                                                                            tolerance, errortolerance, gradtolerance)
                Fpsd2t = np.array(Fpsdt).reshape((1, len(PSDt)))

                print(' ')
                if ok:
                    okc += 1
                    okct += 1
                    okcw += 1
                if perfect:
                    perfectc += 1
                    perfectct += 1
                    perfectcw += 1
                results[nType][nWf] = [Fpsdt, PSDt, Fpsd2t, noise, ok, perfect, grad, noisegrad]
            if iw > 0:
                iwf.append((nWf, perfectcw, iw))
                iwfo.append((nWf, okcw, iw))
        if it > 0:
            itype.append((nType, perfectct, it))
            itypeo.append((nType, okct, it))
    print('******************************************************************************')
    print('*******  NOISE ANALYSIS SUMMARY ********************************************')
    print('******************************************************************************')
    print(' ')
    print('******************************************************************************')
    print('*** PER SUBGROUP *************************************************************')
    print('******************************************************************************')
    for n, (nw, p, t) in enumerate(iwf):
        o = iwfo[n][1]
        print('{}:'.format(nw))
        print('     Perfect PSDs -> {} of {} : {} %'.format(p, t, (p / t) * 100 if t > 0 else 0))
        print('     Noise Fitted OK -> {} of {} : {} %'.format(o, t, (o / t) * 100 if t > 0 else 0))
    print('******************************************************************************')
    print('*** PER GROUP ****************************************************************')
    print('******************************************************************************')
    for n, (nt, p, t) in enumerate(itype):
        o = itypeo[n][1]
        print('{}:'.format(nt))
        print('     Perfect PSDs -> {} of {} : {} %'.format(p, t, (p / t) * 100 if t > 0 else 0))
        print('     Noise Fitted OK -> {} of {} : {} %'.format(o, t, (o / t) * 100 if t > 0 else 0))
    print('******************************************************************************')
    print('*** TOTAL*********************************************************************')
    print('******************************************************************************')
    print('Perfect PSDs -> {} of {} : {} %'.format(perfectc, i, (perfectc / i) * 100 if i > 0 else 0))
    print('Noise Fitted OK -> {} of {} : {} %'.format(okc, i, (okc / i) * 100 if i > 0 else 0))
    print('******************************************************************************')
    print('******* END OF THE NOISE ANALYSIS ********************************************')
    print('******************************************************************************')

    return results


def isPSDok(PSD, Fpsd, noise, tolerance=0.75, errortolerance=0.75, gradtolerance=0.09):
    """

       :param PSD: PSD of a Group
       :param Fpsd: Fpsd of a Group
       :param noise: Noise
       :param tolerance: margin of error allowed on the analysis of the PSD
       :param tolerance: margin of error allowed on the analysis of the Noise
       :return: ok: True if the noise is fitted well
       :return: perfect: True if the PSD gradient is acceptable
       :return: grad: PSD gradient
       :return: noisegrad: Noise gradient

       """
    mPSD = PSD[1:]
    f = Fpsd.transpose()
    dx = np.diff(f)[1:]

    y = np.diff(mPSD)
    y2 = np.diff(noise.transpose())

    grad = qty.Divide(y, dx) / np.max(mPSD)  # Gradient of the mean PSD
    perfect = np.all(np.abs(grad) <= tolerance)

    noisegrad = qty.Divide(y2, dx) / np.max(mPSD)  # Gradient of the noise fitting

    graderror = grad - noisegrad
    error = (mPSD - noise) / np.max(mPSD)
    minerr = np.abs(np.min(error))
    meangraderr = np.mean(graderror)

    ok1 = minerr <= errortolerance
    ok2 = np.abs(meangraderr) <= gradtolerance

    ok = ok1 and ok2

    if perfect:
        print('PSD Noise OK -> {}'.format(np.max(np.abs(grad))))
    else:
        print('PSD Noise BAD -> {}'.format(np.max(np.abs(grad))))

    if ok:
        print('    Noise Fitted OK -> error:{} grad-error:{}'.format(minerr, meangraderr))
    else:
        print('    Noise Fitted BAD -> error:{} grad-error:{}'.format(minerr, meangraderr))
        if not ok1:
            print('         Noise Error BAD -> error:{}'.format(minerr))
        if not ok2:
            print('         Noise GradError BAD -> grad-error:{}'.format(meangraderr))

    return ok, perfect, grad, noisegrad

