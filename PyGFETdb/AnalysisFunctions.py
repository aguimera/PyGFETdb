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

def processAllPSDsPerGroup(rPSD, **kwargs):
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
            Fpsdt = vWf
            ic += 1
            iw += 1
            print('***************************************************')
            print('{}) Group:{}'.format(ic, nWf))
            print('***************************************************')

            [mPSD, noise, ok, perfect, grad, noisegrad] = processAllNoiseGroup(PSDt, Fpsdt, NoAt, NoBt,
                                                                               **kwargs)

            mPSD = mPSD

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


def processAllNoiseGroup(PSD, Fpsd, NoA, NoB, **kwargs):
    temp0 = []
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    temp5 = []
    for i, item in enumerate(PSD):
        [mPSD, noise, ok, perfect, grad, noisegrad] = processAllNoise(item, Fpsd[i], NoA[i], NoB[i],
                                                                      **kwargs)
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
    return [mPSD, noise, ok, perfect, grad, noisegrad]


def processAllNoise(PSD, Fpsd, NoA, NoB, **kwargs):
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
    temp0 = []
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    mPSD = PSD

    if type(NoA) is list():
        for i, item in enumerate(NoA):
            [mPSD, noise, ok, perfect, grad, noisegrad] = processNoA(PSD, Fpsd, NoA[i], NoB[i],
                                                                     **kwargs)

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
        [mPSD, noise, ok, perfect, grad, noisegrad] = processNoA(PSD, Fpsd, NoA, NoB,
                                                                 **kwargs)

    return [mPSD, noise, ok, perfect, grad, noisegrad]


def processNoA(PSD, Fpsd, NoA, NoB,
               **kwargs):
    noise = np.array([])
    ok = False
    perfect = False
    grad = np.array([])
    noisegrad = np.array([])
    mPSD = PSD

    if NoA is not None and len(NoA) > 0:
        if type(NoA) is list:
            return processNoAlist(PSD, Fpsd, NoA, NoB, **kwargs)
        NoA = NoA.transpose()
        NoB = NoB.transpose()

        noise = calculateNoise(Fpsd, NoA, NoB)

        mPSD = reshapePSD(PSD, NoA)

        if mPSD is None:
            return [noise, False, False, [], []]

        if mPSD.ndim == 4:
            ok, perfect, grad, noisegrad = processPSDlist(Fpsd, mPSD, noise, **kwargs)
        elif mPSD.ndim == 5:
            ok, perfect, grad, noisegrad = processListPSDlist(Fpsd, mPSD, noise, **kwargs)
        else:
            raise BaseException("PSD number of dimensions insufficient.")

    return [mPSD, noise, ok, perfect, grad, noisegrad]


def reshapePSD(PSD, NoA):
    if PSD.ndim == 3:
        mPSD = PSD
        n = int(mPSD.size / mPSD.shape[1] / NoA.shape[0] / NoA.shape[1])
        mPSD = mPSD.reshape(n, NoA.shape[0], NoA.shape[1], mPSD.shape[1])
        mPSD = np.array([mPSD])
    elif PSD.ndim == 2:
        try:
            mPSD = PSD.transpose()
            mPSD = mPSD.reshape(NoA.shape[0], NoA.shape[1], mPSD.shape[1])
            mPSD = np.array([[mPSD]])
        except ValueError:
            try:
                mPSD = PSD.transpose()
                n = int(mPSD.shape[0] / NoA.shape[0] / NoA.shape[1])
                mPSD = mPSD.reshape(n, NoA.shape[0], NoA.shape[1], mPSD.shape[1])
                mPSD = np.array([mPSD])
            except ValueError:
                print("Warning: Possible wrong analysis measures. PSDs not equally distributed.")
                mPSD = None
    else:
        mPSD = PSD
        mPSD = np.array([mPSD, mPSD])
        mPSD = np.array([mPSD.transpose()])

    return mPSD


def calculateNoise(Fpsd, NoA, NoB):
    f = np.array([Fpsd])
    tnoise = []
    for i, item in enumerate(NoA):
        NoAi = NoA[i]
        NoBi = NoB[i]
        noise = Fnoise(f.transpose(), [NoAi], [NoBi])
        tnoise.append(noise.transpose())
    noise = np.array(tnoise)
    return noise


def processNoAlist(PSD, Fpsd, NoA, NoB,
                   **kwargs):
    temp0 = []
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    for i in range(0, len(NoA) - 1):
        [noise, ok, perfect, grad, noisegrad] = processNoA(PSD, Fpsd, NoA[i], NoB[i],
                                                           **kwargs)
        temp0.append(noise)
        temp1.append(ok)
        temp2.append(perfect)
        temp3.append(grad)
        temp4.append(noisegrad)

    noise = np.array(temp0)
    ok = np.all(temp1)
    perfect = np.all(temp2)
    grad = temp3
    noisegrad = temp4

    return [noise, ok, perfect, grad, noisegrad]


def processPSDlist(Fpsd, mPSD, noise,
                   **kwargs):
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    for tpsd in mPSD:
        if noise is None:
            ok, perfect, grad, noisegrad = processIsPSDok(tpsd, tpsd, None, Fpsd,
                                                          **kwargs)
            temp1.append(ok)
            temp2.append(perfect)
            temp3.append(grad)
            temp4.append(noisegrad)
        else:
            for i, item in enumerate(noise):
                ok, perfect, grad, noisegrad = processIsPSDok(tpsd, tpsd[i], item, Fpsd,
                                                              **kwargs)
                temp1.append(ok)
                temp2.append(perfect)
                temp3.append(grad)
                temp4.append(noisegrad)
    ok = np.all(temp1)
    perfect = np.all(temp2)
    grad = temp3
    noisegrad = temp4

    return ok, perfect, grad, noisegrad


def processListPSDlist(Fpsd, mPSD, noise, **kwargs):
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    for tpsd in mPSD:
        ok, perfect, grad, noisegrad = processPSDlist(Fpsd, tpsd, noise, **kwargs)
        temp1.append(ok)
        temp2.append(perfect)
        temp3.append(grad)
        temp4.append(noisegrad)
    ok = np.all(temp1)
    perfect = np.all(temp2)
    grad = temp3
    noisegrad = temp4

    return ok, perfect, grad, noisegrad


def processIsPSDok(psd, psdi, noise, Fpsd, **kwargs):
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []

    noisei = None

    if noise is not None:
        noisei = noise.transpose()[1:]

    if psd.ndim > 2:
        for dpsd in psdi:
            ok, perfect, grad, noisegrad = processVgs(dpsd, Fpsd, noisei, **kwargs)
            temp1.append(ok)
            temp2.append(perfect)
            temp3.append(grad)
            temp4.append(noisegrad)
    else:
        ok, perfect, grad, noisegrad = processVgs(psdi, Fpsd, noisei,
                                                  **kwargs)
        temp1.append(ok)
        temp2.append(perfect)
        temp3.append(grad)
        temp4.append(noisegrad)

    return temp1, temp2, temp3, temp4


def processAllPSDsPerSubgroup(GrTypes, rPSD, **kwargs):

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
                Fpsdt = np.array(vWf[0])
                PSDt = np.array(PSDt)
                ic += 1
                it += 1
                iw += 1
                print('***************************************************')
                print('{}) Group:{}, Subgroup:{}'.format(ic, nType, nWf))
                print('***************************************************')

                mPSD = PSDt
                if type(NoAt) is list:
                    temp0 = []
                    temp1 = []
                    temp2 = []
                    temp3 = []
                    temp4 = []
                    for i, item in enumerate(NoAt):
                        [mPSD, noise, ok, perfect, grad, noisegrad] = processAllNoise(PSDt, Fpsdt,
                                                                                      NoAt[i],
                                                                                      NoBt[i],
                                                                                      **kwargs)

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
                    [mPSD, noise, ok, perfect, grad, noisegrad] = processAllNoise(PSDt, Fpsdt, NoAt, NoBt, **kwargs)

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


def processVgs(PSD, Fpsd, noise, **kwargs):
    if noise is None:
        noise = np.array([])
        return isPSDok(PSD, Fpsd, noise, **kwargs)

    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    for item in noise.transpose():
        ok, perfect, grad, noisegrad = isPSDok(PSD, Fpsd, item, **kwargs)

        temp1.append(ok)
        temp2.append(perfect)
        temp3.append(grad)
        temp4.append(noisegrad)
    return temp1, temp2, temp3, temp4


def isPSDok(PSD, Fpsd, noise, fluctuation=38.0, peak=58.95, gradient=2e5, fiterror=10.0, fitgradient=1e3,
            normalization=None, **kwargs):
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
        maxmPSD = np.max(mPSD)
    else:
        maxmPSD = normalization

        # Max fluctuation of the PSD
    fluct = (mPSD - np.min(mPSD)) / maxmPSD
    maxfluct = (np.max(fluct))

    if normalization is None:
        maxfluct = 1 - maxfluct

    perfect1 = np.all(maxfluct < fluctuation)

    # Max peak of the PSD
    pk = (mPSD - np.mean(mPSD)) / maxmPSD
    maxpeak = np.max(pk)
    if normalization is None:
        maxpeak = 1 - maxpeak
    perfect2 = np.all(maxpeak < peak)

    # Gradient of the PSD
    f = Fpsd.transpose()
    dx = np.diff(f)[1:]

    y = np.diff(mPSD)
    grad = qty.Divide(y, dx)
    absgrad = np.abs(grad)
    maxgrad = np.max(absgrad)
    graderror = maxgrad / maxmPSD
    perfect3 = np.all(graderror <= gradient)

    # Error of the noise fitting
    fitpeak = (mPSD - np.mean(noise)) / maxmPSD
    fitmaxerr = np.max(fitpeak)
    if normalization is None:
        fitmaxerr = 1 - fitmaxerr
    ok1 = fitmaxerr < fiterror

    # Gradient of the noise fitting
    f = Fpsd.transpose()
    dx = np.diff(f)[1:noise.shape[0]]

    y2 = np.diff(noise.transpose())
    noisegrad = qty.Divide(y2, dx)
    absnoisegrad = np.abs(noisegrad)

    fitgraderror = absgrad - absnoisegrad
    absgraderror = np.abs(fitgraderror)
    fitmaxgraderror = np.mean(absgraderror) / maxmPSD
    ok2 = fitmaxgraderror <= fitgradient

    # Conditions of the PSD
    perfect = perfect1 and perfect2 and perfect3

    # Conditions of the noise fit
    ok = ok1 and ok2

    printReport(True, perfect, perfect1, perfect2, perfect3, ok, ok1, ok2, maxfluct, maxpeak, graderror, fitmaxerr,
                fitmaxgraderror)

    return ok, perfect, grad, noisegrad


def printReport(debug,
                perfect, perfect1, perfect2, perfect3,
                ok, ok1, ok2,
                maxfluct, maxpeak, graderror,
                fitmaxerr, fitmaxgraderror):
    if debug:
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
