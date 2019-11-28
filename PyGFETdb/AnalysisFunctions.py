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

def processAllPSDsPerGroup(rPSD, fluctuation=0.905, peak=0.35, gradient=0.94, fiterror=0.31, fitgradient=0.09,
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

            [mPSD, noise, ok, perfect, grad, noisegrad] = processAllNoise(PSDt, Fpsdt, NoAt, NoBt,
                                                                          fluctuation, peak, gradient,
                                                                          fiterror, fitgradient, maxpsd)
            mPSD = [mPSD]

            Fpsdt = np.array(Fpsdt)

            print(' ')
            if ok:
                okc += 1
                okcw += 1
            if perfect:
                perfectc += 1
                perfectcw += 1
            results[nWf] = [Fpsdt, mPSD, Fpsdt, noise, ok, perfect, grad, noisegrad]
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


def processAllNoise(PSD, Fpsd, NoA, NoB, fluctuation=0.905, peak=0.355, gradient=0.94, fiterror=0.31,
                    fitgradient=0.09,
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
    temp0 = []
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []

    if type(NoA) is list():
        for i, item in enumerate(NoA):
            [noise, ok, perfect, grad, noisegrad] = processNoA(PSD, Fpsd, NoA[i], NoB[i],
                                                               fluctuation, peak, gradient,
                                                               fiterror,
                                                               fitgradient, normalization)

            temp0.append(noise)
            temp1.append(ok)
            temp2.append(perfect)
            temp3.append(grad)
            temp4.append(noisegrad)

        noise = temp0
        ok = np.all(temp1)
        perfect = np.all(temp2)
        grad = temp3
        noisegrad = temp4
    else:
        [noise, ok, perfect, grad, noisegrad] = processNoA(PSD, Fpsd, NoA, NoB,
                                                           fluctuation, peak, gradient,
                                                           fiterror, fitgradient, normalization)

    return [PSD, noise, ok, perfect, grad, noisegrad]


def processNoA(PSD, Fpsd, NoA, NoB,
               fluctuation=0.905, peak=0.355, gradient=0.94, fiterror=0.31,
               fitgradient=0.09, normalization=None):
    noise = np.array([])
    ok = False
    perfect = False
    grad = 0
    noisegrad = np.array([])

    temp0 = []
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    if NoA is not None and len(NoA) > 0:
        if type(NoA) is list:
            for i in range(0, len(NoA) - 1):
                [noise, ok, perfect, grad, noisegrad] = processNoA(PSD, Fpsd, NoA[i], NoB[i],
                                                                   fluctuation, peak, gradient,
                                                                   fiterror, fitgradient, normalization)
                temp0.append(noise)
                temp1.append(ok)
                temp2.append(perfect)
                temp3.append(grad)
                temp4.append(noisegrad)

            noise = temp0
            ok = np.all(temp1)
            perfect = np.all(temp2)
            grad = temp3
            noisegrad = temp4

            return [noise, ok, perfect, grad, noisegrad]

        try:
            f = np.array(Fpsd).reshape((1, len(Fpsd)))
            noise = Fnoise(f, NoA[:, -len(Fpsd):], NoB[:, -len(Fpsd):])

            retnoise = np.mean(noise.transpose(), 1)
            retnoise = retnoise.reshape(len(retnoise))
            temp0.append(retnoise)
        except:
            try:
                f0 = Fpsd[:NoA.shape[1]]
                f = np.array(f0).reshape((1, len(f0)))
                noise = Fnoise(f, NoA, NoB)

                retnoise = np.array([np.zeros(len(Fpsd))])
                temp0.append(retnoise)
            except:
                noise = None  # TODO: Better Fit Calculation? or Bad Fit Data?

                retnoise = np.array([np.zeros(len(Fpsd))])
                temp0.append(retnoise)
        try:
            if PSD.ndim == 3:
                PSD = PSD.reshape(PSD.shape[0], NoA.shape[0], NoA.shape[1])
                mPSD = PSD[:, :, -noise.shape[1]:]
            elif PSD.ndim == 2:
                mPSD = PSD.transpose()
                mPSD = mPSD.reshape(NoA.shape[0], NoA.shape[1], mPSD.shape[1])
                mPSD = np.array([mPSD])
            else:
                mPSD = PSD
                mPSD = np.array([mPSD, mPSD])
                mPSD = np.array([mPSD.transpose()])
        except ValueError:
            print("Warning: Possible wrong analysis. Noise fit data insufficient.")
            mPSD = np.array([PSD.transpose()])

        for tpsd in mPSD:
            if noise is None:
                ok, perfect, grad, noisegrad = processIsPSDok(tpsd, tpsd, None, Fpsd,
                                                              fluctuation, peak, gradient,
                                                              fiterror, fitgradient, normalization)
                temp1.append(ok)
                temp2.append(perfect)
                temp3.append(grad)
                temp4.append(noisegrad)
            else:
                for i, item in enumerate(noise):
                    ok, perfect, grad, noisegrad = processIsPSDok(tpsd, tpsd[i], item, Fpsd,
                                                                  fluctuation, peak, gradient,
                                                                  fiterror, fitgradient, normalization)
                    temp1.append(ok)
                    temp2.append(perfect)
                    temp3.append(grad)
                    temp4.append(noisegrad)
        noise = temp0
        ok = np.all(temp1)
        perfect = np.all(temp2)
        grad = temp3
        noisegrad = temp4

    return [noise, ok, perfect, grad, noisegrad]


def processIsPSDok(psd, psdi, noise, Fpsd,
                   fluctuation=0.905, peak=0.355, gradient=0.94, fiterror=0.31,
                   fitgradient=0.09, normalization=None):
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []

    noisei = None

    if noise is not None:
        noisei = noise.transpose()[1:]

    if psd.ndim > 2:
        for dpsd in psdi:
            [ok, perfect, grad, noisegrad] = isPSDok(dpsd, Fpsd, noisei,
                                                     fluctuation, peak, gradient, fiterror, fitgradient,
                                                     normalization)
            temp1.append(ok)
            temp2.append(perfect)
            temp3.append(grad)
            temp4.append(noisegrad)
    else:
        [ok, perfect, grad, noisegrad] = isPSDok(psdi, Fpsd, noisei,
                                                 fluctuation, peak, gradient, fiterror, fitgradient,
                                                 normalization)
        temp1.append(ok)
        temp2.append(perfect)
        temp3.append(grad)
        temp4.append(noisegrad)

    return temp1, temp2, temp3, temp4


def processAllPSDsPerSubgroup(GrTypes, rPSD, fluctuation=0.905, peak=0.35, gradient=0.94, fiterror=0.31,
                              fitgradient=0.09,
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
                        [mPSD, noise, ok, perfect, grad, noisegrad] = processAllNoise(PSDt, Fpsdt, NoAt[i],
                                                                                      NoBt[i],
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
                Fpsdt = np.array(Fpsdt).reshape(len(PSDt))

                print(' ')
                if ok:
                    okc += 1
                    okct += 1
                    okcw += 1
                if perfect:
                    perfectc += 1
                    perfectct += 1
                    perfectcw += 1
                results[nType][nWf] = [Fpsdt, mPSD, Fpsdt, noise, ok, perfect, grad, noisegrad]
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
    if PSD is None or (type(PSD) is not list and type(PSD) is not np.ndarray):
        print('Wrong PSD measurement.')
        return False, False, [], []

    mPSD = PSD[1:]

    if normalization is None:
        maxmPSD = 1
    else:
        maxmPSD = normalization

    try:
        # Max fluctuation of the PSD
        fluct = np.abs((mPSD - np.abs(np.min(mPSD))) / maxmPSD)
        maxfluct = (np.max(fluct))
        perfect1 = np.all(maxfluct > fluctuation)
    except:
        perfect1 = False
        maxfluct = np.nan

    try:
        # Max peak of the PSD
        pk = np.abs((mPSD - np.abs(np.mean(mPSD))) / maxmPSD)
        maxpeak = np.max(pk)
        perfect2 = np.all(maxpeak > peak)
    except:
        perfect2 = False
        maxpeak = np.nan

    try:
        # Gradient of the PSD
        f = Fpsd.transpose()
        dx = np.diff(f)[1:]

        y = np.diff(mPSD)
        grad = qty.Divide(y, dx)
        absgrad = np.abs(grad)
        maxgrad = np.max(absgrad)
        graderror = maxgrad / maxmPSD
        perfect3 = np.all(graderror <= gradient)
    except:
        perfect3 = False
        graderror = np.nan
        grad = []
        absgrad = np.nan

    try:
        # Error of the noise fitting
        fitpeak = np.abs((mPSD - np.abs(np.mean(noise))) / maxmPSD)
        fitmaxerr = np.max(fitpeak)
        ok1 = fitmaxerr > fiterror
    except:
        ok1 = False
        fitmaxerr = np.nan

    try:
        # Gradient of the noise fitting
        f = Fpsd.transpose()
        dx = np.diff(f)[1:noise.size]

        y2 = np.diff(noise.transpose())
        noisegrad = qty.Divide(y2, dx)
        absnoisegrad = np.abs(noisegrad)

        fitgraderror = absgrad[:absnoisegrad.shape[0]] - absnoisegrad
        absgraderror = np.abs(fitgraderror)
        fitmaxgraderror = np.mean(absgraderror) / maxmPSD
        ok2 = fitmaxgraderror <= fitgradient
    except:
        ok2 = False
        fitmaxgraderror = np.nan
        noisegrad = []

    # Conditions of the PSD
    perfect = perfect1 and perfect2 and perfect3

    # Conditions of the noise fit
    ok = ok1 and ok2

    # Print output
    print(' ')
    if perfect:
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
        print('    Noise Fitted OK ->')
        if ok1:
            print('         Noise Error OK -> error:{}'.format(fitmaxerr))
        if ok2:
            print('         Noise GradError OK -> grad-error:{}'.format(fitmaxgraderror))
    else:
        print('    Noise Fitted BAD ->')
        if not ok1:
            print('         Noise Error BAD -> error:{}'.format(fitmaxerr))
        if not ok2:
            print('         Noise GradError BAD -> grad-error:{}'.format(fitmaxgraderror))
    print(' ')

    return ok, perfect, grad, noisegrad

