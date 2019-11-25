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
        for i in range(1, 28):  # To widen the effect increase the 25
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


def processNoise(PSD, Fpsd, NoA, NoB, tolerance=0.75, errortolerance=0.75, gradtolerance=0.09):
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


def isMeanPSDok(PSD, Fpsd, noise, tolerance=0.75, errortolerance=0.75, gradtolerance=0.09):
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
    perfect = np.all(np.abs(grad) <= tolerance)

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
        print('Mean PSD Noise PERFECT -> {}'.format(np.max(np.abs(grad))))
    else:
        print('Mean PSD Noise BAD -> {}'.format(np.max(np.abs(grad))))

    if ok:
        print('Noise Fitted OK -> error:{} grad-error:{}'.format(minerr, meangraderr))
    else:
        print('Noise Fitted BAD -> error:{} grad-error:{}'.format(minerr, meangraderr))
        if not ok1:
            print('     Noise Error BAD -> error:{}'.format(minerr))
        if not ok2:
            print('     Noise GradError BAD -> grad-error:{}'.format(meangraderr))

    return ok, perfect, grad, noisegrad


def processAllNoise(PSD, Fpsd, NoA, NoB, fluctuation=0.905, peak=0.355, gradient=0.94, fiterror=0.31, fitgradient=0.09):
    """

    :param PSD: PSD of a Group
    :param Fpsd: Fpsd of a Group
    :param NoA: Noise A of a Group
    :param NoB: Noise B of a Group
    :param fluctuation: Maximum change allowed between the maximum value and the minimum value
    :param peak: Maximum change allowed between the maximum value and the mean value
    :param gradient: Maximum gradient allowed
    :param fiterror: Maximum error allowed in the fit
    :param fitgradient: Maximum error allowed in the gradient of the fit
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
    mPSD = None
    if NoA is not None and len(NoA) > 0:
        f = np.array(Fpsd).reshape((1, len(PSD)))
        noise = Fnoise(f, NoA[:, -len(f):], NoB[:, -len(f):])
        if PSD.ndim == 3:
            PSD = PSD.reshape(PSD.shape[0], NoA.shape[0], NoA.shape[1])
            mPSD = PSD[:, :, -noise.shape[1]:]
        elif PSD.ndim == 2:
            mPSD = np.array([PSD.transpose()])
        else:
            mPSD = np.array([PSD])
            mPSD = mPSD.reshape(PSD.shape[0], NoA.shape[0], NoA.shape[1])
            mPSD = mPSD[:, :, -noise.shape[1]:]
        for psd in mPSD:
            for i, item in enumerate(noise):
                [ok, perfect, grad, noisegrad] = isPSDok(psd[i], Fpsd, item.transpose()[1:],
                                                     fluctuation, peak, gradient, fiterror, fitgradient)
                temp1.append(ok)
                temp2.append(perfect)
                temp3.append(grad)
                temp4.append(noisegrad)

        noise = noise.reshape(len(noise), len(PSD))

        ok = np.all(temp1)
        perfect = np.all(temp2)
        grad = temp3  # grad.reshape(len(PSD))
        noisegrad = temp4  #noisegrad.reshape(len(PSD))

    return [mPSD, noise, ok, perfect, grad, noisegrad]


def processAllPSDs(GrTypes, rPSD, fluctuation=0.905, peak=0.35, gradient=0.94, fiterror=0.31, fitgradient=0.09):

    """

    :param GrTypes: Group to process
    :param rPSD: Results of a param search of PSD, Fpsd, NoA and NoB
    :param fluctuation: Maximum change allowed between the maximum value and the minimum value
    :param peak: Maximum change allowed between the maximum value and the mean value
    :param gradient: Maximum gradient allowed
    :param fiterror: Maximum error allowed in the fit
    :param fitgradient: Maximum error allowed in the gradient of the fit
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
                    temp5 = []
                    for i, item in enumerate(NoAt):
                        [mPSD, noise, ok, perfect, grad, noisegrad] = processAllNoise(PSDt, Fpsdt, NoAt[i], NoBt[i],
                                                                                      fluctuation, peak, gradient,
                                                                                      fiterror,
                                                                                      fitgradient)
                        temp0.append(noise)
                        temp1.append(ok)
                        temp2.append(perfect)
                        temp3.append(grad)
                        temp4.append(noisegrad)
                        temp5.append(mPSD)

                    noise = temp0
                    ok = np.all(temp1)
                    perfect = np.all(temp2)
                    grad = temp3
                    noisegrad = temp4
                    mPSD = temp5
                else:
                    [mPSD, noise, ok, perfect, grad, noisegrad] = processAllNoise(PSDt, Fpsdt, NoAt, NoBt,
                                                                                  fluctuation, peak, gradient,
                                                                                  fiterror, fitgradient)
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
                results[nType][nWf] = [Fpsdt, mPSD, Fpsd2t, noise, ok, perfect, grad, noisegrad]
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


def isPSDok(PSD, Fpsd, noise, fluctuation=0.188, peak=0.256, gradient=0.354, fiterror=0.5, fitgradient=0.009):
    """
       :param PSD: PSD of a Group
       :param Fpsd: Fpsd of a Group
       :param noise: Noise
       :param fluctuation: Maximum change allowed between the maximum value and the minimum value
       :param peak: Maximum change allowed between the maximum value and the mean value
       :param gradient: Maximum gradient allowed
       :param fiterror: Maximum error allowed in the fit
       :param fitgradient: Maximum error allowed in the gradient of the fit

       :return: ok: True if the noise is fitted well
       :return: perfect: True if the PSD gradient is acceptable
       :return: grad: PSD gradient
       :return: noisegrad: Noise gradient

       """
    mPSD = PSD[1:]
    f = Fpsd.transpose()
    dx = np.diff(f)[1:]

    # Max fluctuation of the PSD
    fluct = np.abs((mPSD - np.abs(np.min(mPSD))) / np.max(mPSD))
    maxfluct = (1 - np.max(fluct))

    # Max peak of the PSD
    pk = np.abs((mPSD - np.abs(np.mean(mPSD))) / np.max(mPSD))
    maxpeak = 1 - np.max(pk)

    # Gradient of the PSD
    y = np.diff(mPSD)
    grad = qty.Divide(y, dx)
    absgrad = np.abs(grad)
    maxgrad = np.max(absgrad)
    graderror = maxgrad

    # Error of the noise fitting
    fitpeak = np.abs((mPSD - np.abs(np.mean(noise))) / np.max(mPSD))
    fitmaxerr = 1 - np.max(fitpeak)

    # Gradient of the noise fitting
    y2 = np.diff(noise.transpose())
    noisegrad = qty.Divide(y2, dx)
    absnoisegrad = np.abs(noisegrad)
    fitgraderror = absgrad - absnoisegrad
    absgraderror = np.abs(fitgraderror)
    fitmaxgraderror = np.mean(absgraderror)

    # Consistency check
    if maxfluct > 1 or maxpeak > 1 or graderror > 1 or fitmaxerr > 1 or fitmaxgraderror > 1:
        raise ValueError(
            'Value>1 -> f:{} p:{} g:{} fe:{} fg:{}'.format(maxfluct, maxpeak, graderror, fitmaxerr, fitmaxgraderror))

    # Conditions of the PSD
    perfect1 = np.all(maxfluct <= fluctuation)
    perfect2 = np.all(maxpeak <= peak)
    perfect3 = np.all(graderror <= gradient)
    perfect = perfect1 and perfect2 and perfect3

    # Conditions of the noise fit
    ok1 = fitmaxerr <= fiterror
    ok2 = fitmaxgraderror <= fitgradient
    ok = ok1 and ok2

    # Print output
    if perfect:
        pass
        #print('PSD Noise OK -> fluctuation:{} peak:{} gradient:{}'.format(maxfluct, maxpeak, graderror))
    else:
        #print('PSD Noise BAD -> fluctuation:{} peak:{} gradient:{}'.format(maxfluct, maxpeak, graderror))
        if not perfect1:
            print('         PSD FLUCTUATION BAD -> fluctuation:{}'.format(maxfluct))
        if not perfect2:
            print('         PSD PEAK BAD -> peak:{}'.format(maxpeak))
        if not perfect3:
            print('         PSD GRADIENT BAD -> gradient:{}'.format(graderror))

    if ok:
        pass
        #print('    Noise Fitted OK -> error:{} grad-error:{}'.format(fitmaxerr, fitmaxgraderror))
    else:
        #print('    Noise Fitted BAD -> error:{} grad-error:{}'.format(fitmaxerr, fitmaxgraderror))
        if not ok1:
            print('         Noise Error BAD -> error:{}'.format(fitmaxerr))
        if not ok2:
            print('         Noise GradError BAD -> grad-error:{}'.format(fitmaxgraderror))

    return ok, perfect, grad, noisegrad

