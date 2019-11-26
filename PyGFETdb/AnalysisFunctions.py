# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Analysis Functions that do not fit in the previous files.

"""

import numpy as np

from PyGFETdb import qty, GlobalFunctions as g
from PyGFETdb.NoiseModel import Fnoise


########################################################################
#
#  FREQUENCY FILTERS
#
########################################################################
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


########################################################################
#
#  PSD ANALYSIS
#
########################################################################

def processAllNoise(PSD, Fpsd, NoA, NoB, fluctuation=0.905, peak=0.355, gradient=0.94, fiterror=0.31, fitgradient=0.09,
                    normalization=None):
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
    if NoA is not None and len(NoA) > 0:
        f = np.array(Fpsd[:NoA.shape[0]])
        noise = Fnoise(f, NoA[:len(f), :], NoB[:len(f), :])
        fpsd = np.array(Fpsd[:1])
        for psd in PSD:
            tnoise = noise.transpose()[1:]
            [ok, perfect, grad, noisegrad] = isPSDok(psd, fpsd, tnoise.transpose(),
                                                     fluctuation, peak, gradient, fiterror, fitgradient,
                                                     normalization)
            temp1.append(ok)
            temp2.append(perfect)
            temp3.append(grad)
            temp4.append(noisegrad)

        noise = np.mean(noise.transpose(), 1)
        noise = noise.reshape(noise.size)

        ok = np.all(temp1)
        perfect = np.all(temp2)
        grad = temp3
        noisegrad = temp4

    return [PSD, noise, ok, perfect, grad, noisegrad]


def processAllPSDsPerDevice(rPSD, fluctuation=0.905, peak=0.35, gradient=0.94, fiterror=0.31, fitgradient=0.09,
                            normalization=None):
    """

    :param rPSD: Results of a param search of PSD, Fpsd, NoA and NoB
    :param fluctuation: Maximum change allowed between the maximum value and the minimum value
    :param peak: Maximum change allowed between the maximum value and the mean value
    :param gradient: Maximum gradient allowed
    :param fiterror: Maximum error allowed in the fit
    :param fitgradient: Maximum error allowed in the gradient of the fit
    :return: A dict with the results of the processing
    """
    results = {}
    ic = 0
    itypes = 0
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

    if normalization is not None:
        maxpsd = normalization
    else:
        maxpsd = 1

    Fpsd = rPSD.get('Fpsd')
    results = {}
    PSD = rPSD['PSD']
    NoA = rPSD['NoA']
    NoB = rPSD['NoB']
    for nWf, vWf in Fpsd.items():
        PSDt = PSD[nWf]
        NoAt = NoA[nWf]
        NoBt = NoB[nWf]
        perfectcw = 0
        okcw = 0
        iw = 0
        if NoAt is not None and len(NoAt) > 0:
            Fpsdt = vWf[0:len(PSDt)]
            ic += 1
            iw += 1
            print('***************************************************')
            print('{}) Group:{}'.format(ic, nWf))
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
                                                                                  fitgradient, maxpsd)
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
                                                                              fiterror, fitgradient, maxpsd)
            Fpsd2t = np.array(Fpsdt)

            print(' ')
            if ok:
                okc += 1
                okcw += 1
            if perfect:
                perfectc += 1
                perfectcw += 1
            results[nWf] = [Fpsdt, mPSD, Fpsd2t, noise, ok, perfect, grad, noisegrad]
        if iw > 0:
            iwf.append((nWf, perfectcw, iw))
            iwfo.append((nWf, okcw, iw))
    print('******************************************************************************')
    print('*******  NOISE ANALYSIS SUMMARY ********************************************')
    print('******************************************************************************')
    print(' ')
    print('******************************************************************************')
    print('*** PER GROUP *************************************************************')
    print('******************************************************************************')
    for iw, (nw, p, t) in enumerate(iwf):
        o = iwfo[iw][1]
        print('{}:'.format(nw))
        print('     Perfect PSDs -> {} of {} : {} %'.format(p, t, (p / t) * 100 if t > 0 else 0))
        print('     Noise Fitted OK -> {} of {} : {} %'.format(o, t, (o / t) * 100 if t > 0 else 0))
    print('******************************************************************************')
    print('*** TOTAL*********************************************************************')
    print('******************************************************************************')
    print('Perfect PSDs -> {} of {} : {} %'.format(perfectc, ic, (perfectc / ic) * 100 if ic > 0 else 0))
    print('Noise Fitted OK -> {} of {} : {} %'.format(okc, ic, (okc / ic) * 100 if ic > 0 else 0))
    print('******************************************************************************')
    print('******* END OF THE NOISE ANALYSIS ********************************************')
    print('******************************************************************************')

    return results



def processAllPSDs(GrTypes, rPSD, fluctuation=0.905, peak=0.35, gradient=0.94, fiterror=0.31, fitgradient=0.09,
                   normalization=None):

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
    ic = 0
    itypes = 0
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

    if normalization is not None:
        maxpsd = normalization
    else:
        maxpsd = 1

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
                ic += 1
                it += 1
                iw += 1
                print('***************************************************')
                print('{}) Group:{}, Subgroup:{}'.format(ic, nType, nWf))
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
                                                                                      fitgradient, maxpsd)
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
                    mPSD = PSDt
                else:
                    [mPSD, noise, ok, perfect, grad, noisegrad] = processAllNoise(PSDt, Fpsdt, NoAt, NoBt,
                                                                                  fluctuation, peak, gradient,
                                                                                  fiterror, fitgradient, maxpsd)
                Fpsd2t = np.array(Fpsdt)

                print(' ')
                if ok:
                    okc += 1
                    okct += 1
                    okcw += 1
                if perfect:
                    perfectc += 1
                    perfectct += 1
                    perfectcw += 1
                results[nType][nWf] = [Fpsdt[0], mPSD, Fpsd2t[0], noise, ok, perfect, grad, noisegrad]
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
    for iw, (nw, p, t) in enumerate(iwf):
        o = iwfo[iw][1]
        print('{}:'.format(nw))
        print('     Perfect PSDs -> {} of {} : {} %'.format(p, t, (p / t) * 100 if t > 0 else 0))
        print('     Noise Fitted OK -> {} of {} : {} %'.format(o, t, (o / t) * 100 if t > 0 else 0))
    print('******************************************************************************')
    print('*** PER GROUP ****************************************************************')
    print('******************************************************************************')
    for it, (nt, p, t) in enumerate(itype):
        o = itypeo[it][1]
        print('{}:'.format(nt))
        print('     Perfect PSDs -> {} of {} : {} %'.format(p, t, (p / t) * 100 if t > 0 else 0))
        print('     Noise Fitted OK -> {} of {} : {} %'.format(o, t, (o / t) * 100 if t > 0 else 0))
    print('******************************************************************************')
    print('*** TOTAL*********************************************************************')
    print('******************************************************************************')
    print('Perfect PSDs -> {} of {} : {} %'.format(perfectc, ic, (perfectc / ic) * 100 if ic > 0 else 0))
    print('Noise Fitted OK -> {} of {} : {} %'.format(okc, ic, (okc / ic) * 100 if ic > 0 else 0))
    print('******************************************************************************')
    print('******* END OF THE NOISE ANALYSIS ********************************************')
    print('******************************************************************************')

    return results


def isPSDok(PSD, Fpsd, noise, fluctuation=38.0, peak=58.95, gradient=2e5, fiterror=10.0, fitgradient=1e3,
            normalization=1e-22):
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

    mPSD = PSD
    dx = np.diff(Fpsd)

    if normalization is None:
        maxmPSD = 1
    else:
        maxmPSD = normalization

    # Max fluctuation of the PSD
    fluct = np.abs((mPSD - np.abs(np.min(mPSD))) / maxmPSD)
    maxfluct = (np.max(fluct))

    # Max peak of the PSD
    pk = np.abs((mPSD - np.abs(np.mean(mPSD))) / maxmPSD)
    maxpeak = np.max(pk)

    # Gradient of the PSD
    y = np.diff(mPSD)

    if y.ndim == 2:
        grad = qty.Divide(y[:, :dx.shape[1]], dx[:, :y.shape[1]])
    else:
        grad = qty.Divide(y, dx)

    absgrad = np.abs(grad)
    maxgrad = np.max(absgrad)
    graderror = maxgrad / maxmPSD

    # Error of the noise fitting
    fitpeak = np.abs((mPSD - np.abs(np.mean(noise))) / maxmPSD)
    fitmaxerr = np.max(fitpeak)

    # Gradient of the noise fitting
    y2 = np.diff(noise.transpose())

    if y2.ndim == 2:
        noisegrad = qty.Divide(y2[:, :dx.shape[1]], dx[:, :y2.shape[1]])
    else:
        noisegrad = qty.Divide(y2, dx[1:])

    absnoisegrad = np.abs(noisegrad)
    if absgrad.ndim == 2:
        fitgraderror = absgrad[1:, :absnoisegrad.shape[1]] - absnoisegrad[:, :absgrad.shape[1]]
    else:
        fitgraderror = absgrad[1:] - absnoisegrad
    absgraderror = np.abs(fitgraderror)
    fitmaxgraderror = np.mean(absgraderror) / maxmPSD

    # Conditions of the PSD
    perfect1 = np.all(maxfluct > fluctuation)
    perfect2 = np.all(maxpeak > peak)
    perfect3 = np.all(graderror <= gradient)
    perfect = perfect1 and perfect2 and perfect3

    # Conditions of the noise fit
    ok1 = fitmaxerr > fiterror
    ok2 = fitmaxgraderror <= fitgradient
    ok = ok1 and ok2

    # Print output
    print(' ')
    if perfect:
        pass

        print('PSD Noise OK ->')
        if perfect1:
            print('         PSD FLUCTUATION OK -> fluctuation:{}'.format(maxfluct))
        if perfect2:
            print('         PSD PEAK OK -> peak:{}'.format(maxpeak))
        if perfect3:
            print('         PSD GRADIENT OK -> gradient:{}'.format(graderror))

    else:
        print('PSD Noise BAD ->')
        if not perfect1:
            print('         PSD FLUCTUATION BAD -> fluctuation:{}'.format(maxfluct))
        if not perfect2:
            print('         PSD PEAK BAD -> peak:{}'.format(maxpeak))
        if not perfect3:
            print('         PSD GRADIENT BAD -> gradient:{}'.format(graderror))

    if ok:
        print('    Noise Fitted OK -> error:{} grad-error:{}'.format(fitmaxerr, fitmaxgraderror))
        if ok1:
            print('         Noise Error OK -> error:{}'.format(fitmaxerr))
        if ok2:
            print('         Noise GradError OK -> grad-error:{}'.format(fitmaxgraderror))
    else:
        print('    Noise Fitted BAD -> error:{} grad-error:{}'.format(fitmaxerr, fitmaxgraderror))
        if not ok1:
            print('         Noise Error BAD -> error:{}'.format(fitmaxerr))
        if not ok2:
            print('         Noise GradError BAD -> grad-error:{}'.format(fitmaxgraderror))
    print(' ')

    return ok, perfect, grad, noisegrad

