# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Analysis Functions that do not fit in the previous files.

"""

import numpy as np

import PyGFETdb
from PyGFETdb import qty, GlobalFunctions as g, Thread as mp
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


def process60Hz(Array, process):
    """
        **Removes the frequency 60Hz**

    :param Array: array of frequencies
    :param process: bool that activates the function
    :return: the array without the frequency 50Hz
    """
    if process:
        # remove 60Hz
        for i in range(1, 2):  # To widen the effect increase the 2
            Array = g.remove(Array, 50)
    return Array


def processBelow10Hz(Array, process):
    """
        **Removes the frequencies below 1Hz**

    :param Array: array of frequencies
    :param process: bool that activates the function
    :return: the array without the frequencies below 1Hz
    """
    if process:
        #  remove below 1Hz
        for i in range(1, 35):  # To widen the effect increase the 15
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
    Array = process60Hz(Array, process)
    Array = processBelow10Hz(Array, process)
    Array = processHiFreqs(Array, process)
    return Array


########################################################################
#
#  PSD ANALYSIS
#
########################################################################

def processAllPSDsPerGroup(rPSD, HaltOnFail=False, meanpsd=True,**kwargs):
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
    perfect = False
    perfectc = 0
    print(' ')
    print('******************************************************************************')
    print('******* RESULTS OF THE NOISE ANALYSIS ****************************************')
    print('******************************************************************************')
    print(' ')

    #kwargs.update({'HaltOnFail': HaltOnFail})

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
                                                                               HaltOnFail=HaltOnFail,
                                                                               **kwargs)

            Fpsdt = np.array(Fpsdt)

            if np.all(perfect) and meanpsd:
                try:
                    meanPSD = np.mean(PSDt, 0)
                    perfect2, ok2 = isMeanPSDOk(Fpsdt[0], meanPSD, noise, **kwargs)
                    ok = np.all(ok) and ok2
                    perfect = np.all(perfect) and perfect2
                except ValueError:
                    print()
                    print("**WARNING** : Mean PSD not valid. Skipping check.")

            printReportPerSubgroup(perfect, ok)

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

    PrintSummaryPerGroup(iwf, iwfo, perfectc, okc, ic)

    return results


def PrintSummaryPerGroup(iwf, iwfo, perfectc, okc, ic):
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


def processAllNoiseGroup(PSD, Fpsd, NoA, NoB, HaltOnFail=False, **kwargs):
    temp0 = []
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    temp5 = []
    #kwargs.update({'HaltOnFail': HaltOnFail})

    for i, item in enumerate(PSD):
        mPSD, noise, ok, perfect, grad, noisegrad = processAllNoise(item, Fpsd[i], NoA[i], NoB[i],
                                                                    HaltOnFail=HaltOnFail,
                                                                    **kwargs)
        temp0.append(noise)
        temp1.append(ok)
        temp2.append(perfect)
        temp3.append(grad)
        temp4.append(noisegrad)
        temp5.append(mPSD)
        if not np.all(temp2) and HaltOnFail:
            print()
            print('Bad results: Halting analysis.')
            print()
            noise = np.array(temp0)
            noise = np.mean(noise, 0)
            ok = np.all(temp1)
            perfect = np.all(temp2)
            grad = temp3
            noisegrad = temp4
            mPSD = temp5
            return [mPSD, noise, ok, perfect, grad, noisegrad]
    noise = np.array(temp0)
    noise = np.mean(noise, 0)
    ok = np.all(temp1)
    perfect = np.all(temp2)
    grad = temp3
    noisegrad = temp4
    mPSD = temp5

    return [mPSD, noise, ok, perfect, grad, noisegrad]


def processAllNoise(PSD, Fpsd, NoA, NoB, **kwargs):
    try:
        PSD = np.array(PSD)
        n = int(PSD.size / 500 + 5)
    except ValueError:
        n = len(PSD) * 20 + 5
    # print('Processing PSDs... forking {} threads'.format(n))
    thread = mp.MultiProcess(PyGFETdb.AnalysisFunctions, n)
    key = thread.initcall(mp.key(), PyGFETdb.AnalysisFunctions)
    args = {'PSD': PSD, 'Fpsd': Fpsd, 'NoA': NoA, 'NoB': NoB}
    args.update(kwargs)
    thread.call(key, PyGFETdb.AnalysisFunctions, '_processAllNoise', args, **args)
    r = thread.getResults(key)[key][0]
    thread.end()
    del thread
    return r


def _processAllNoise(PSD, Fpsd, NoA, NoB, HaltOnFail=False,
                     **kwargs):
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
    temp5 = []

    ##kwargs.update({'HaltOnFail': HaltOnFail})

    if type(NoA) is list():
        for i, item in enumerate(NoA):
            mPSD, noise, ok, perfect, grad, noisegrad = processNoA(PSD, Fpsd, NoA[i], NoB[i],
                                                                   HaltOnFail=HaltOnFail,
                                                                   **kwargs)

            temp0.append(noise)
            temp1.append(ok)
            temp2.append(perfect)
            temp3.append(grad)
            temp4.append(noisegrad)
            temp5.append(mPSD)
            if not np.all(temp2) and HaltOnFail:
                noise = np.array(temp0)
                noise = np.mean(noise, 0)
                ok = np.all(temp1)
                perfect = np.all(temp2)
                grad = temp3
                noisegrad = temp4
                mPSD = temp5
                return mPSD, noise, ok, perfect, grad, noisegrad
        noise = np.array(temp0)
        noise = np.mean(noise, 0)
        ok = temp1
        perfect = temp2
        grad = temp3
        noisegrad = temp4
        mPSD = temp5
    else:
        mPSD, noise, ok, perfect, grad, noisegrad = processNoA(PSD, Fpsd, NoA, NoB,
                                                               HaltOnFail=HaltOnFail,
                                                               **kwargs)

    return mPSD, noise, ok, perfect, grad, noisegrad


def isMeanPSDOk(Fpsd, PSD, noise, meanfluctuation=None, meanpeak=None, meangradient=None, meangradientmean=None,
                fiterror=None, fitgradient=None,**kwargs):


    meanPSD = np.mean(PSD, 1)
    meannoise = np.mean(noise, 0)
    meannoise = meannoise.reshape(meannoise.size)[1:]

    ok, perfect, grad2, noisegrad2 = isPSDok(meanPSD, Fpsd, meannoise.transpose(),
                                             meanfluctuation, meanpeak, meangradient,
                                             meangradientmean, fiterror, fitgradient,printbad=True,printok=False)

    print()
    if perfect:
        print('Analysing Mean PSD... OK.')
    else:
        print('Analysing Mean PSD... BAD.')
    if ok:
        print('Mean PSD Fit... OK.')
    else:
        print('Mean PSD Fit... BAD.')
    print()
    return perfect, ok


def processNoA(PSD, Fpsd, NoA, NoB, HaltOnFail=False, **kwargs):
    retnoise = np.array([])
    grad = np.array([])
    noisegrad = np.array([])
    mPSD = PSD
    ok = False
    perfect = False

    #kwargs.update({'HaltOnFail': HaltOnFail})

    if NoA is not None and len(NoA) > 0:
        if type(NoA) is list:
            return processNoAlist(PSD, Fpsd, NoA, NoB,
                                  HaltOnFail=HaltOnFail,
                                  **kwargs)

        noise, retnoise = calculateNoise(Fpsd, NoA, NoB)

        mPSD = reshapePSD(PSD, NoA)

        if mPSD is None:
            return PSD, retnoise, False, False, [], []

        try:
            if mPSD.ndim == 4:
                ok, perfect, grad, noisegrad = processPSDlist(Fpsd, mPSD, noise,
                                                              HaltOnFail=HaltOnFail,
                                                              **kwargs)
            elif mPSD.ndim == 5:
                ok, perfect, grad, noisegrad = processListPSDlist(Fpsd, mPSD, noise,
                                                                  HaltOnFail=HaltOnFail,
                                                                  **kwargs)
            else:
                raise BaseException("PSD number of dimensions insufficient.")
        except AttributeError:
            temp0 = []
            # temp1=[]
            temp2 = []
            temp3 = []
            temp4 = []
            temp5 = []
            for item in mPSD:
                tmPSD, retnoise, ok, perfect, grad, noisegrad = processNoA(item, Fpsd, NoA, NoB,
                                                                           HaltOnFail=HaltOnFail,
                                                                           **kwargs)


                temp0.append(tmPSD)
                # temp1.append(retnoise)
                temp2.append(ok)
                temp3.append(perfect)
                temp4.append(grad)
                temp5.append(noisegrad)
                if not np.all(temp3) and HaltOnFail:
                    mPSD = temp0
                    ok = temp2
                    perfect = temp3
                    grad = temp4
                    noisegrad = temp5
                    return mPSD, noise, ok, perfect, grad, noisegrad

            mPSD = temp0
            ok = temp2
            perfect = temp3
            grad = temp4
            noisegrad = temp5

    return mPSD, retnoise, ok, perfect, grad, noisegrad


def reshapePSD(PSD, NoA):
    mPSD = None
    try:
        if PSD is not None:
            mPSD = _reshapePSD(PSD, NoA)
    except AttributeError:
        mPSD = []
        for item in PSD:
            r = reshapePSD(item, NoA)
            mPSD.append(r)
    return mPSD


def _reshapePSD(PSD, NoA):
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
                # print("Warning: Possible wrong analysis measures. PSDs not equally distributed.")
                mPSD = None
    elif PSD.ndim == 4:
        mPSD = PSD
    elif PSD.ndim == 5:
        mPSD = PSD
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
        if noise is not None:
            tnoise.append(noise.transpose())
    if len(tnoise) > 0:
        noise = np.array(tnoise)
        retnoise = np.mean(noise, 0)
    else:
        noise = np.array([])
        retnoise = np.repeat(0, len(Fpsd))
    return noise, retnoise


def processNoAlist(PSD, Fpsd, NoA, NoB, HaltOnFail=False,
                   **kwargs):
    temp0 = []
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    #kwargs.update({'HaltOnFail': HaltOnFail})

    for i in range(0, len(NoA) - 1):
        noise, ok, perfect, grad, noisegrad = processNoA(PSD, Fpsd, NoA[i], NoB[i],
                                                         HaltOnFail=HaltOnFail,
                                                         **kwargs)
        temp0.append(noise)
        temp1.append(ok)
        temp2.append(perfect)
        temp3.append(grad)
        temp4.append(noisegrad)
        if not np.all(temp2) and HaltOnFail:
            noise = np.array(temp0)
            noise = np.mean(noise, 0)
            ok = temp1
            perfect = temp2
            grad = temp3
            noisegrad = temp4
            return noise, ok, perfect, grad, noisegrad

    noise = np.array(temp0)
    noise = np.mean(noise, 0)
    ok = temp1
    perfect = temp2
    grad = temp3
    noisegrad = temp4

    return noise, ok, perfect, grad, noisegrad


def processPSDlist(Fpsd, mPSD, noise, HaltOnFail=False,
                   **kwargs):
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    #kwargs.update({'HaltOnFail': HaltOnFail})

    for tpsd in mPSD:
        if noise is None:
            ok, perfect, grad, noisegrad = processIsPSDok(tpsd, tpsd, None, Fpsd,
                                                          HaltOnFail=HaltOnFail,
                                                          **kwargs)
            temp1.append(ok)
            temp2.append(perfect)
            temp3.append(grad)
            temp4.append(noisegrad)
            if not np.all(temp2) and HaltOnFail:
                return temp1, temp2, temp3, temp4
        else:
            if noise.shape[0] != tpsd.shape[0]:
                tpsd = tpsd.transpose()
            for i, item in enumerate(noise):
                ok, perfect, grad, noisegrad = processIsPSDok(tpsd, tpsd[i], item, Fpsd,
                                                              HaltOnFail=HaltOnFail,
                                                              **kwargs)
                temp1.append(ok)
                temp2.append(perfect)
                temp3.append(grad)
                temp4.append(noisegrad)
                if not np.all(temp2) and HaltOnFail:
                    return temp1, temp2, temp3, temp4

    return temp1, temp2, temp3, temp4


def processListPSDlist(Fpsd, mPSD, noise, HaltOnFail=False, **kwargs):
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    #kwargs.update({'HaltOnFail': HaltOnFail})

    for tpsd in mPSD:
        ok, perfect, grad, noisegrad = processPSDlist(Fpsd, tpsd, noise,                                                                               HaltOnFail=HaltOnFail,
 **kwargs)

        temp1.append(ok)
        temp2.append(perfect)
        temp3.append(grad)
        temp4.append(noisegrad)
        if not np.all(temp2) and HaltOnFail:
            ok = temp1
            perfect = temp2
            grad = temp3
            noisegrad = temp4
            return ok, perfect, grad, noisegrad
    ok = temp1
    perfect = temp2
    grad = temp3
    noisegrad = temp4

    return ok, perfect, grad, noisegrad


def processIsPSDok(psd, psdi, noise, Fpsd, HaltOnFail=False, **kwargs):
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    #kwargs.update({'HaltOnFail': HaltOnFail})

    noisei = None

    if noise is not None:
        noisei = noise.transpose()[1:]

    if psd.ndim > 2:
        if noisei.shape[0] + 1 != psdi.shape[1]:
            psdi = psdi.transpose()
        for dpsd in psdi:
            ok, perfect, grad, noisegrad = processVgs(dpsd, Fpsd, noisei,
                                                      HaltOnFail=HaltOnFail,
                                                      **kwargs)
            temp1.append(ok)
            temp2.append(perfect)
            temp3.append(grad)
            temp4.append(noisegrad)
            if not np.all(temp2) and HaltOnFail:
                ok = temp1
                perfect = temp2
                grad = temp3
                noisegrad = temp4

                return ok, perfect, grad, noisegrad

    else:
        ok, perfect, grad, noisegrad = processVgs(psdi, Fpsd, noisei,
                                                  HaltOnFail=HaltOnFail,
                                                  **kwargs)
        temp1=ok
        temp2=perfect
        temp3=grad
        temp4=noisegrad

    ok = temp1
    perfect = temp2
    grad = temp3
    noisegrad = temp4

    return ok, perfect, grad, noisegrad


def processAllPSDsPerSubgroup(GrTypes, rPSD,  **kwargs):

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

        ic, it, okc, okct, perfectc, perfectct, iwf, iwfo, itype, itypeo, results = \
            processSubGroup(Fpsd, PSD, NoA, NoB, ic, it, okc, okct, perfectc, perfectct, iwf, iwfo,
                            itype, itypeo, results,
                            nType, **kwargs)
        if it > 0:
            itype.append((nType, perfectct, it))
            itypeo.append((nType, okct, it))

    PrintSummaryPerSubGroup(iwf, iwfo, itype, itypeo, perfectc, okc, ic)

    return results


def processSubGroup(Fpsd, PSD, NoA, NoB, ic, it, okc, okct, perfectc, perfectct, iwf, iwfo,
                    itype, itypeo, results,
                    nType, HaltOnFail=False, meanpsd=True, **kwargs):
    #kwargs.update({'HaltOnFail': HaltOnFail})

    for nWf, vWf in Fpsd.items():
        PSDt = PSD[nWf]
        NoAt = NoA[nWf]
        NoBt = NoB[nWf]
        perfectcw = 0
        okcw = 0
        iw = 0
        perfect = False
        if NoAt is not None and len(NoAt) > 0:
            Fpsdt = np.array(vWf[0])
            PSDt = [PSDt]
            ic += 1
            it += 1
            iw += 1
            print('***************************************************')
            print('{}) Group:{}, Subgroup:{}'.format(ic, nType, nWf))
            print('***************************************************')

            #kwargs.update({'HaltOnFail': HaltOnFail})

            mPSD, noise, ok, perfect, grad, noisegrad = processPSDSubgroup(Fpsdt, PSDt,
                                                                           NoAt, NoBt,
                                                                           HaltOnFail=HaltOnFail,
                                                                           **kwargs)
            if np.all(perfect) and meanpsd:
                try:
                    meanPSD = np.mean(PSDt, 0)
                    perfect2, ok2 = isMeanPSDOk(Fpsdt[0], meanPSD, noise,
                                            HaltOnFail=HaltOnFail,
                                            **kwargs)
                    perfect = np.all(perfect) and perfect2
                    ok = np.all(ok) and ok2
                except ValueError:
                    print()
                    print("**WARNING** : Mean PSD not valid. Skipping check.")

            printReportPerSubgroup(perfect, ok)

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
        if not np.all(perfect) and HaltOnFail:
            PrintSummaryPerSubGroup(iwf, iwfo, itype, itypeo, perfectc, okc, ic)
            print()
            print('Bad results: Halting analysis.')
            print()
            return ic, it, okc, okct, perfectc, perfectct, iwf, iwfo, itype, itypeo, results
    return ic, it, okc, okct, perfectc, perfectct, iwf, iwfo, itype, itypeo, results


def PrintSummaryPerSubGroup(iwf, iwfo, itype, itypeo, perfectc, okc, ic):
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


def processPSDSubgroup(Fpsdt, PSDt, NoAt, NoBt, HaltOnFail=None, **kwargs):
    #kwargs.update({'HaltOnFail': HaltOnFail})
    if type(NoAt) is list:
        temp0 = []
        temp1 = []
        temp2 = []
        temp3 = []
        temp4 = []
        temp5 = []
        for i, item in enumerate(NoAt):
            mPSD, noise, ok, perfect, grad, noisegrad = processAllNoise(PSDt, Fpsdt,
                                                                        NoAt[i],
                                                                        NoBt[i],
                                                                        HaltOnFail=HaltOnFail,
                                                                        **kwargs)

            temp0.append(noise)
            temp1.append(ok)
            temp2.append(perfect)
            temp3.append(grad)
            temp4.append(noisegrad)
            temp5.append(mPSD)
            if not np.all(temp2) and HaltOnFail:
                noise = np.array(temp0)
                ok = temp1
                perfect = temp2
                grad = temp3
                noisegrad = temp4
                mPSD = temp5
                return mPSD, noise, ok, perfect, grad, noisegrad
        noise = np.array(temp0)
        try:
            noise = np.mean(noise, 0)
        except ValueError:
            "None value in noise fit calculation."
            pass
        ok = temp1
        perfect = temp2
        grad = temp3
        noisegrad = temp4
        mPSD = temp5
    else:
        mPSD, noise, ok, perfect, grad, noisegrad = processAllNoise(PSDt, Fpsdt, NoAt, NoBt,
                                                                    HaltOnFail=HaltOnFail,
                                                                    **kwargs)
    return mPSD, noise, ok, perfect, grad, noisegrad


def processVgs(PSD, Fpsd, noise, HaltOnFail=False, **kwargs):

    #kwargs.update({'HaltOnFail': HaltOnFail})

    if noise is None:
        noise = np.array([])
        ok, perfect, grad, noisegrad = isPSDok(PSD, Fpsd, noise, **kwargs)
        if not np.all(perfect) and HaltOnFail:
            return ok, perfect, grad, noisegrad

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
        if not np.all(temp2) and HaltOnFail:
            ok = temp1
            perfect = temp2
            grad = temp3
            noisegrad = temp4
            return ok, perfect, grad, noisegrad

    ok = temp1
    perfect = temp2
    grad = temp3
    noisegrad = temp4

    return ok, perfect, grad, noisegrad


def isPSDok(PSD, Fpsd, noise, fluctuation=43e-3, peak=0.58, gradient=5e-19, gradientmean=0.5,
            fiterror=0.3, fitgradient=5e-21,
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
    fluct = (mPSD - np.abs(np.min(mPSD))) / maxmPSD
    maxfluct = (np.max(fluct))

    if normalization is None:
        maxfluct = 1 - maxfluct

    perfect1 = np.all(maxfluct < fluctuation)

    # Max peak of the PSD
    pk = np.abs(mPSD - np.abs(np.mean(mPSD))) / maxmPSD
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
    graderror = maxgrad
    perfect3 = np.all(graderror <= gradient)

    # Mean Gradient of the PSD
    mgrad = np.abs(np.mean(grad))
    perfect4 = np.all(mgrad <= gradientmean)

    # Error of the noise fitting
    fitpeak = np.abs(mPSD - np.abs(np.mean(noise))) / maxmPSD
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

    fitgraderror = np.mean(absgrad) - np.mean(absnoisegrad)
    absgraderror = np.abs(fitgraderror)
    fitmaxgraderror = np.mean(absgraderror)

    ok2 = np.all(fitmaxgraderror <= fitgradient)

    # Conditions of the PSD
    perfect = perfect1 and perfect2 and perfect3 and perfect4

    # Conditions of the noise fit
    ok = ok1 and ok2

    printReport(perfect, perfect1, perfect2, perfect3, perfect4, ok, ok1, ok2,
                maxfluct, maxpeak, graderror, mgrad,
                fitmaxerr, fitmaxgraderror,**kwargs)

    return ok, perfect, grad, noisegrad


def printReport(
        perfect, perfect1, perfect2, perfect3, perfect4, ok, ok1, ok2,
        maxfluct, maxpeak, graderror, mgrad,
        fitmaxerr, fitmaxgraderror, printbad=True, printok=True,**kwargs):
    # Print output
    if perfect and printok:
        print()
        print('PSD Noise OK ->')
        if perfect1:
            print('         PSD FLUCTUATION OK -> fluctuation:{}'.format(maxfluct))
        if perfect2:
            print('         PSD PEAK OK -> peak:{}'.format(maxpeak))
        if perfect3:
            print('         PSD GRADIENT OK -> gradient:{}'.format(graderror))
        if perfect4:
            print('         PSD MEAN GRADIENT OK -> mean-gradient:{}'.format(mgrad))
    if not perfect and printbad:
        print()
        print('PSD Noise BAD ->')
        if not perfect1:
            print('         PSD FLUCTUATION BAD -> fluctuation:{}'.format(maxfluct))
        if not perfect2:
            print('         PSD PEAK BAD -> peak:{}'.format(maxpeak))
        if not perfect3:
            print('         PSD GRADIENT BAD -> gradient:{}'.format(graderror))
        if not perfect4:
            print('         PSD MEAN GRADIENT BAD -> mean-gradient:{}'.format(mgrad))

    if ok and printok:
        print('    Noise Fitted OK ->')
        if ok1:
            print('         Noise Error OK -> error:{}'.format(fitmaxerr))
        if ok2:
            print('         Noise GradError OK -> grad-error:{}'.format(fitmaxgraderror))
            print()
    if not ok and printbad:
        print('    Noise Fitted BAD ->')
        if not ok1:
            print('         Noise Error BAD -> error:{}'.format(fitmaxerr))
        if not ok2:
            print('         Noise GradError BAD -> grad-error:{}'.format(fitmaxgraderror))
            print()


def printReportPerSubgroup(perfect, ok):
    if perfect:
        print()
        print('PSD Noise OK.')
    else:
        print()
        print('PSD Noise BAD.')
    if ok:
        print('    Noise Fitted OK.')
        print()
    else:
        print('    Noise Fitted BAD.')
        print()
