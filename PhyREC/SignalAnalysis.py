#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 11:12:23 2018

@author: aguimera
"""
import numpy as np
import matplotlib.pyplot as plt
import quantities as pq
from scipy import signal
from collections import OrderedDict
import PhyREC.PlotWaves as Rplt
import neo
from PhyREC.NeoInterface import NeoSegment, NeoSignal


def nFFTFMin(Fs, Fmin):
    return int(2**(np.around(np.log2(Fs/Fmin))+1))


def threshold_detection(signal, threshold=0.0 * pq.mV, sign='above',
                        RelaxTime=None):
    """
    Returns the times when the analog signal crosses a threshold.
    Usually used for extracting spike times from a membrane potential.
    Adapted from version in NeuroTools.

    Parameters
    ----------
    signal : neo AnalogSignal object
        'signal' is an analog signal.
    threshold : A quantity, e.g. in mV
        'threshold' contains a value that must be reached
        for an event to be detected. Default: 0.0 * mV.
    sign : 'above' or 'below'
        'sign' determines whether to count thresholding crossings
        that cross above or below the threshold.
    format : None or 'raw'
        Whether to return as SpikeTrain (None)
        or as a plain array of times ('raw').

    Returns
    -------
    result_st : neo SpikeTrain object
        'result_st' contains the spike times of each of the events (spikes)
        extracted from the signal.
    """

    assert threshold is not None, "A threshold must be provided"

    if sign == 'above':
        cutout = np.where(signal > threshold)[0]
    elif sign == 'below':
        cutout = np.where(signal < threshold)[0]

    if len(cutout) <= 0:
        events = np.zeros(0)
    else:
        take = np.where(np.diff(cutout) > 1)[0] + 1
        take = np.append(0, take)

        time = signal.times.rescale('s')
        events = time[cutout][take]

    if RelaxTime:
        outevents = []
        told = 0*pq.s
        for te in events:
            if (te-told) > RelaxTime:
                outevents.append(te)
                told = te
        outevents = np.array(outevents)*pq.s

    else:
        outevents = events

    return outevents


def PlotSpectralCoherence(RefSig, Signals, Time=None, nFFT=2**17, FMin=None, Ax=None,
                          scaling='density', Units=None, Noverlap=2, **LineKwargs):

    if Ax is None:
        Fig, Ax = plt.subplots()

    if FMin is not None:
        nFFT = nFFTFMin(RefSig.sampling_rate, FMin)

    noverlap = nFFT/Noverlap

    csACC = []
    for sl in Signals:
        if not hasattr(sl, 'GetSignal'):
            continue
        
        f, cs = signal.coherence(RefSig,
                                 sl.GetSignal(Time, Units=Units),
                                 fs=RefSig.sampling_rate,
                                 nperseg=nFFT,
                                 noverlap=noverlap,
                                 axis=0)

        if hasattr(sl, 'LineKwargs'):
            lkwargs = sl.LineKwargs.copy()
            lkwargs.update(LineKwargs)
        else:
            lkwargs = LineKwargs

        if 'label' not in lkwargs:
            lkwargs['label'] = sl.name

        Ax.loglog(f, cs, **lkwargs)        
        csACC.append(cs)
    return f, csACC


def PlotPSD(Signals, Time=None, nFFT=2**17, FMin=None, Ax=None,
            scaling='density', Units=None, Noverlap=2, **LineKwargs):

    if Ax is None:
        Fig, Ax = plt.subplots()

    PSD = {}        
    for sl in Signals:
        if hasattr(sl, 'GetSignal'):            
            sig = sl.GetSignal(Time, Units=Units)
        else:
            sig = sl.time_slice(*Time)

        if FMin is not None:
            nFFT = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/FMin))+1))

        noverlap = nFFT/Noverlap

        ff, psd = signal.welch(x=sig, fs=sig.sampling_rate, axis=0,
                               window='hanning',
                               nperseg=nFFT,
                               noverlap=noverlap,
                               scaling=scaling)

        if scaling == 'density':
            units = sig.units**2/pq.Hz
        elif scaling == 'spectrum':
            units = sig.units**2

        slN = sl.name
        PSD[slN] = {}
        PSD[slN]['psd'] = psd * units
        PSD[slN]['ff'] = ff

        if hasattr(sl, 'LineKwargs'):
            lkwargs = sl.LineKwargs.copy()
            lkwargs.update(LineKwargs)
        else:
            lkwargs = {}
            lkwargs.update(LineKwargs)

        if 'label' not in lkwargs:
            lkwargs['label'] = slN
        
        Ax.loglog(ff, psd, **lkwargs)

    Ax.set_xlabel('Frequency [Hz]')
    Ax.set_ylabel('[' + str(units).split(' ')[-1] + ']')

    handles, labels = Ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))

#    nLines = len(by_label)
#    nlc = 4
#    if nLines > nlc:
#        ncol = (nLines / nlc) + ((nLines % nlc) > 0)
#    else:
#        ncol = 1
    Ax.legend(by_label.values(), by_label.keys(),
              loc='best',
#              ncol=1,
              fontsize='x-small')

    return PSD


def PlotEventAvg(Signals, TimesEvent, TimeAvg, Time=None,
                 OverLap=True, Std=True, Units=None,
                 FileOutPrefix=None, SpecSlot=None):

    if len(Signals)>1:
        ft, Axt = plt.subplots()
    else:
        ft = None
        Axt = None

    for sl in Signals:
        avg = np.array([])
        if SpecSlot is None:
            fig, Ax = plt.subplots()
            Ax.set_xlabel('Time [s]')
        else:
            fig, ((Ax, ac),(Axs, Axc)) = plt.subplots(2,
                                          2,
                                          sharex=True,
                                          gridspec_kw={'width_ratios': (10, 1)})    
            ac.axis('off')
            Axc.axis('off')
            Axs.set_xlabel('Time [s]')
            Axs.set_ylabel('Frequency [Hz]')
            
        Ts = sl.Signal.sampling_period
        nSamps = int((TimeAvg[1]-TimeAvg[0])/Ts)
        t = np.arange(nSamps)*Ts + TimeAvg[0]

        for et in TimesEvent:
            start = et+TimeAvg[0]
            stop = et+TimeAvg[1]

            st = sl.GetSignal((start, stop), Units=Units)[:nSamps]
            try:
                avg = np.hstack([avg, st]) if avg.size else st
                if OverLap:
                    Ax.plot(t, st, 'k-', alpha=0.1)
            except:
                print ('Error', nSamps, et, avg.shape, st.shape)

        MeanT = np.mean(avg, axis=1)

        if SpecSlot is not None:
            SpecSlot.Ax = Axs
            SpecSlot.CAx = Axc
            SpecSlot.Signal = NeoSignal(Signal=st.duplicate_with_new_array(signal=MeanT*st.units))
            SpecSlot.Signal.signal.t_start = t[0]
            SpecSlot.DispName = SpecSlot.Signal.Name
            SpecSlot.PlotSignal(None)

        Ax.plot(t, MeanT, 'r-')
        if Std:
            StdT = np.std(avg, axis=1)
            Ax.fill_between(t, MeanT+StdT, MeanT-StdT,
                            facecolor='r', alpha=0.5)

        ylim = Ax.get_ylim()
        Ax.vlines(0, ylim[0], ylim[1],
                  linestyles='dashdot',
                  alpha=0.7,
                  colors='g')
        Ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

#        plt.title(sl.Name)
        
        su = str(st.units).split(' ')[-1]
#        Ax.set_ylabel(sl.Name + ' [' + su + ']')
        Ax.set_ylabel(sl.Name)

        if Axt is not None:
            Axt.plot(t, MeanT, label=sl.Name)

        if FileOutPrefix is not None:
            fig.savefig(FileOutPrefix + sl.Name + '.png')
        

    if Axt is not None:
        ylim = Axt.get_ylim()
        Axt.vlines(0, ylim[0], ylim[1],
                   linestyles='dashdot',
                   alpha=0.7,
                   colors='g')
        Axt.set_ylabel('[' + su + ']')
        Axt.set_xlabel('Time [s]')
        Axt.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        Axt.legend()

    if FileOutPrefix is not None:
        ft.savefig(FileOutPrefix + '.png')


#        if Spect:
#            nFFT = int(2**(np.around(np.log2(Fs/sl.SpecFmin))+1))
#            noverlap = int((Ts*nFFT - sl.SpecTimeRes)/Ts)
#            TWindOff = (nFFT * Ts)/8
#
#            f, tsp, Sxx = signal.spectrogram(MeanT, Fs,
#                                             window='hanning',
#                                             nperseg=nFFT,
#                                             noverlap=noverlap,
#                                             scaling='density',
#                                             axis=0)
#
#            finds = np.where(f < sl.SpecFmax)[0][1:]
#            print Sxx.shape
#            r, c = Sxx.shape
#            S = Sxx.reshape((r, c))[finds][:]
#            pcol = AxS.pcolormesh(tsp + TimeWindow[0].magnitude + TWindOff,
#                                  f[finds],
#                                  np.log10(S),
#                                  vmin=np.log10(np.max(S))+sl.SpecMinPSD,
#                                  vmax=np.log10(np.max(S)),
#                                  cmap=sl.SpecCmap)
#            f, a = plt.subplots(1, 1)
#            f.colorbar(pcol)
#
#        Axt.plot(t, np.mean(avg, axis=1), label=sl.DispName)
#        ft.canvas.draw()
#
#    Axt.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
#    Axt.legend()
#    ft.canvas.draw()
#    plt.show()

#
#    def PlotHist(self, Time, Resamp=False, Binds=250):
#        fig, Ax = plt.subplots()
#
#        for sl in self.Slots:
#            Ax.hist(sl.GetSignal(Time, Resamp=Resamp),
#                    Binds,
#                    alpha=0.5)
#            Ax.set_yscale('log')
#
#    def PlotPSD(self, Time, nFFT=2**17, FMin=None, Resamp=False):
#
#        if not self.FigFFT or not plt.fignum_exists(self.FigFFT.number):
#            self.FigFFT, self.AxFFT = plt.subplots()
#
#        PSD = {}
#        for sl in self.Slots:
#            sig = sl.GetSignal(Time, Resamp=Resamp)
#            if FMin:
#                nFFT = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/FMin))+1)) 
#
#            ff, psd = signal.welch(x=sig, fs=sig.sampling_rate,
#                                   window='hanning',
#                                   nperseg=nFFT, scaling='density', axis=0)
#            slN = sl.SigName
#            PSD[slN] = {}
#            PSD[slN]['psd'] = psd
#            PSD[slN]['ff'] = ff
#            self.AxFFT.loglog(ff, psd, label=sl.DispName)
#
#        self.AxFFT.set_xlabel('Frequency [Hz]')
#        self.AxFFT.set_ylabel('PSD [V^2/Hz]')
#        self.AxFFT.legend()
#        return PSD



#
#    def PlotPSDSNR(self, (EventRec, EventName), TDelay, TEval, DevDCVals,
#                   Time=None, nFFT=2**17):
#
#        etimes = EventRec.GetEventTimes(EventName, Time)
#
#        DevACVals = {}
#
#        for sl in self.Slots:
#            fig, (ax1, ax2)  = plt.subplots(1,2)
#            psd = np.array([])
#            DCVals = DevDCVals[sl.SigName]
#            for ne, et in enumerate(etimes):
#                start = et+TDelay
#                stop = et+TDelay+TEval
#
#                sig = sl.GetSignal((start, stop), Resamp=False)
#                fpsd, npsd = signal.welch(x=sig, fs=sig.sampling_rate,
#                                          window='hanning',
#                                          nperseg=nFFT, scaling='density', axis=0)                
#                psd = np.hstack([psd, npsd]) if psd.size else npsd
#
##                Flin = fpsd[1:].magnitude
##                Flog = np.logspace(np.log10(Flin[0]),
##                                   np.log10(Flin[-1]), 100)
##                Flog = np.round(Flog, 9)
##                Flin = np.round(Flin, 9)
##                intpsd = interpolate.interp1d(Flin, npsd[1:].transpose())(Flog)
##                a, b, _ = FETAna.noise.FitNoise(Flog,
##                                                intpsd, Fmin=150, Fmax=5e3)
#
#                ax1.loglog(fpsd, npsd)
##                ax2.loglog(Flog, intpsd/FETAna.noise(Flog, a, b))
#
#            psdD = {'Vd0': psd.transpose()}
#            ACVals = {'PSD': psdD,
#                      'gm': None,
#                      'Vgs': DCVals['Vgs'],
#                      'Vds': DCVals['Vds'],
#                      'Fpsd': fpsd.magnitude,
#                      'Fgm': None,
#                      'ChName': sl.SigName,
#                      'Name': sl.DispName,
#                      'GMPoly': DCVals['GMPoly'],
#                      'IdsPoly': DCVals['IdsPoly'],
#                      'Ud0': DCVals['Ud0'],
#                      'IsOK': DCVals['IsOK'],
#                      'DateTime': DCVals['DateTime']}
#            DevACVals[sl.DispName] = ACVals
#            fig.canvas.draw()
#            plt.show()
#
#        FETAna.InterpolatePSD(DevACVals)
#        FETAna.FitACNoise(DevACVals, Fmin=150, Fmax=5e3, IsOkFilt=False)
#        pltPSD = FETplt.PyFETPlot()
#        pltPSD.AddAxes(('PSD', 'NoA', 'NoB'))
#        pltPSD.AddLegend()
#        pltPSD.PlotDataCh(DevACVals)
#
#        pickle.dump(DevACVals, open('test.pkl', 'wb'))
#
#        return DevACVals



def PlotPSD_SNR(Signals, TimeSig=None, TimeNoise=None, nFFT=2**17, FMin=None, Ax=None,
            scaling='density', Units=None, **LineKwargs):

    if Ax is None:
        Fig, Ax = plt.subplots()

    PSD = {}        
    for sl in Signals:
        if not hasattr(sl, 'GetSignal'):
            continue
        sig = sl.GetSignal(TimeSig, Units=Units)
        noise = sl.GetSignal(TimeNoise, Units=Units)

        if FMin is not None:
            nFFTsig = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/FMin))+1))
            nFFTnoise = int(2**(np.around(np.log2(noise.sampling_rate.magnitude/FMin))+1))

        ff, psdsig = signal.welch(x=sig, fs=sig.sampling_rate, axis=0,
                               window='hanning',
                               nperseg=nFFTsig,
                               scaling=scaling)
        
        ff, psdnoise = signal.welch(x=noise, fs=noise.sampling_rate, axis=0,
                               window='hanning',
                               nperseg=nFFTnoise,
                               scaling=scaling)
        
        SNR=psdsig/psdnoise
#        if scaling == 'density':
#            units = sig.units**2/pq.Hz
#        elif scaling == 'spectrum':
#            units = sig.units**2

        slN = sl.name
        PSD[slN] = {}
        PSD[slN]['SNR'] = SNR 
        PSD[slN]['ff'] = ff

        if hasattr(sl, 'LineKwargs'):
            lkwargs = sl.LineKwargs.copy()
            lkwargs.update(LineKwargs)
        else:
            lkwargs = LineKwargs

        if 'label' not in lkwargs:
            lkwargs['label'] = slN
        
        Ax.loglog(ff, SNR, **lkwargs)

    Ax.set_xlabel('Frequency [Hz]')
    Ax.set_ylabel('SNR')

    handles, labels = Ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))

#    nLines = len(by_label)
#    nlc = 4
#    if nLines > nlc:
#        ncol = (nLines / nlc) + ((nLines % nlc) > 0)
#    else:
#        ncol = 1
    Ax.legend(by_label.values(), by_label.keys(),
              loc='best',
#              ncol=1,
              fontsize='x-small')

    return PSD