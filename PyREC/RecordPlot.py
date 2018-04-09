#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:11:53 2017

@author: aguimera
"""

import neo
import numpy as np
import matplotlib.pyplot as plt
import quantities as pq
from scipy import signal
from scipy import interpolate
import pickle
from fractions import Fraction
from matplotlib import cm
import matplotlib.colors as colors


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

        time = signal.times
        events = time[cutout][take]

    if RelaxTime:
        outevents = []
        told = 0*pq.s
        for te in events:
            if (te-told) > RelaxTime:
                outevents.append(te)
                told = te
    else:
        outevents = events

    return outevents


class Filter():

    FTypes = {'lp': 'lowpass',
              'hp': 'highpass',
              'bp': 'bandpass',
              'bs': 'bandstop'}

    def __init__(self, Type, Order, Freq1, Freq2=None):
        self.Type = Type
        self.Freq1 = Freq1
        self.Freq2 = Freq2
        self.Order = Order

    def ApplyFilter(self, sig):
        st = np.array(sig)
        for nf, typ in enumerate(self.Type):
            if typ == 'lp' or typ == 'hp':
                FType = self.FTypes[typ]
                Freqs = self.Freq1[nf]/(0.5*sig.sampling_rate)
            elif typ == 'bp' or typ == 'bs':
                FType = self.FTypes[typ]
                Freqs = [self.Freq1[nf]/(0.5*sig.sampling_rate),
                         self.Freq2[nf]/(0.5*sig.sampling_rate)]
            else:
                print 'Filter Type error ', typ
                continue

#            print nf, self.Order[nf], Freqs, FType
            b, a = signal.butter(self.Order[nf], Freqs, FType)
            st = signal.filtfilt(b, a, st, axis=0)

        return neo.AnalogSignal(st,
                                units=sig.units,
                                t_start=sig.t_start,
                                sampling_rate=sig.sampling_rate,
                                name=sig.name)


class PltSlot():
    TrtNameStr = 'T'

    def __init__(self):
        self.SpecFmin = 0.5
        self.SpecFmax = 15
        self.SpecTimeRes = 0.05
        self.SpecMinPSD = -3
        self.SpecCmap = 'jet'

        self.Position = None
        self.DispName = None
        self.SigName = None
        self.FileName = None
        self.ResamplePoints = 10000
        self.ResampleFs = None
        self.PlotType = 'Wave'
        self.Ax = None
        self.RecN = 1
        self.TStart = None #  0*pq.s
        self.AutoScale = True
        self.FiltType = ('', )
        self.FiltOrder = ('', )
        self.FiltF1 = ('', )
        self.FiltF2 = ('', )
        self.Filter = None
        self.OutType = 'I'  # 'I' 'V' 'Vg'
        self.IVGain = 1e4
        self.Gm = -1e-4
        self.GmSignal = ''
        self.rec = None

        self.LSB = None

        self.Color = 'k'
        self.Line = '-'
        self.Alpha=1

        self.Ymax = 0
        self.Ymin = 0
        
        self.Fig = None

    def SetTstart(self):
        if self.TStart is not None:
            self.rec.SetTstart(self.SigName, self.TStart)

    def Resample(self, sig):
        if self.ResampleFs:
#            print 'Resamp freq', self.ResampleFs
            f = self.ResampleFs/sig.sampling_rate
            fact = Fraction(float(f)).limit_denominator()
            dowrate = fact.denominator
            uprate = fact.numerator
        else:
#            print 'Resamp points', self.ResamplePoints
            dowrate = sig.times.shape[0]/self.ResamplePoints
            if dowrate > 0:
                f = float(1/float(dowrate))
                uprate = 1

        if dowrate > 0:
            print sig.sampling_rate*f, f, uprate, dowrate
            rs = signal.resample_poly(sig, uprate, dowrate)
            sig = neo.AnalogSignal(rs,
                                   units=sig.units,
                                   t_start=sig.t_start,
                                   sampling_rate=sig.sampling_rate*f,
                                   name= sig.name)

        if self.LSB:
            rs = np.round(sig/self.LSB)
            rs = rs * self.LSB
            sig = neo.AnalogSignal(rs,
                                   units=sig.units,
                                   t_start=sig.t_start,
                                   sampling_rate=sig.sampling_rate,
                                   name=sig.name)

        return sig

    def Signal(self):
        return self.rec.Signal(self.SigName)

    def CheckTime(self, Time):
        if Time is None:
            return (self.Signal().t_start, self.Signal().t_stop)

        if Time[0] < self.Signal().t_start:
            Tstart = self.Signal().t_start
        else:
            Tstart = Time[0]

        if Time[1] > self.Signal().t_stop:
            Tstop = self.Signal().t_stop
        else:
            Tstop = Time[1]

        return (Tstart, Tstop)

    def ApplyGain(self, sig):
        if self.OutType == 'I':
            return sig/self.IVGain
        elif self.OutType == 'Vg':
            if self.GmSignal:
                gm = self.GmSignal[0].Signal(self.GmSignal[1])
                ti = gm.times
                gma = np.array(gm)

                if sig.t_start < gm.t_start:
                    ti = np.hstack((sig.t_start, ti))
                    gma = np.hstack((gma[0, None],
                                     gma.transpose())).transpose()

                if sig.t_stop > gm.t_stop:
                    ti = np.hstack((ti, sig.t_stop))
                    gma = np.hstack((gma.transpose(),
                                     gma[-1, None])).transpose()

                Gm = interpolate.interp1d(ti, gma, axis=0)(sig.times)
                return sig/(self.IVGain*Gm)

            return sig/(self.IVGain*self.Gm)
        else:
            return sig

    def GetSignal(self, Time, Resamp=True):
        sig = self.rec.GetSignal(self.SigName,
                                 self.CheckTime(Time))
        sig = self.ApplyGain(sig)

        if self.Filter:
            sig = self.Filter.ApplyFilter(sig)

        if Resamp:
            return self.Resample(sig)
        else:
            return sig

    def PlotSignal(self, Time, Resamp=True):
        if self.Ax:
            if self.PlotType == 'Spectrogram':
                sig = self.GetSignal(Time, Resamp=False)

                nFFT = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/self.SpecFmin))+1))
                Ts = sig.sampling_period.magnitude
                noverlap = int((Ts*nFFT - self.SpecTimeRes)/Ts)

                f, t, Sxx = signal.spectrogram(sig, sig.sampling_rate,
                                               window='hanning',
                                               nperseg=nFFT,
                                               noverlap=noverlap,
                                               scaling='density',
                                               axis=0)
                finds = np.where(f < self.SpecFmax)[0][1:]
                r, g, c = Sxx.shape
                S = Sxx.reshape((r, c))[finds][:]
                cax = self.Ax.pcolor(t*pq.s+sig.t_start, f[finds], S,
                                     cmap=self.SpecCmap,
                                     norm=colors.LogNorm(self.SpecMinPSD,
                                                         self.SpecMaxPSD))
                fig = plt.figure()
                fig.colorbar(cax)
#                self.Fig.colorbar(cax, orientation='horizontal')

#                self.Ax.set_yscale('log')
                
            elif self.PlotType == 'Wave':
                sig = self.GetSignal(Time, Resamp)
                self.Ax.plot(sig.times,
                             sig,
                             self.Line,
                             color=self.Color,
                             label=self.DispName,alpha=self.Alpha)

                if not(self.Ymin == 0 and self.Ymax == 0):
                    ylim = (self.Ymin, self.Ymax)
                    self.Ax.set_ylim(ylim)

    def PlotEvent(self, Time, color=None, alpha=0.5):
        if color is None:
            color = 'r--'
        self.Ax.plot((Time, Time), (1, -1), color, alpha=alpha)


class PlotRecord():
    FigFFT = None
    AxFFT = None
    Axs = None  # (Axhandler, Autoscale)
    LegNlabCol = 4  # Number of labels per col in legend
    LegFontSize = 'x-small'
    ShowAxis = True

    def ClearAxes(self):
        for sl in self.Slots:
            while sl.Ax.lines:
                sl.Ax.lines[0].remove()

    def FormatFigure(self):
        if self.ShowAxis:            
            for (Ax, _) in self.Axs:
    #            Ax.get_yaxis().set_visible(False)
                Ax.get_xaxis().set_visible(False)
                Ax.spines['top'].set_visible(False)
                Ax.spines['right'].set_visible(False)
    #            Ax.spines['left'].set_visible(False)
                Ax.spines['bottom'].set_visible(False)
                Ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        
            TimeAx = self.Axs[-1][0]
            TimeAx.set_xlabel('Time [s]')
            TimeAx.get_xaxis().set_visible(True)
        else:
            for (Ax, _) in self.Axs:
                Ax.get_yaxis().set_visible(False)
                Ax.get_xaxis().set_visible(False)
                Ax.spines['top'].set_visible(False)
                Ax.spines['right'].set_visible(False)
                Ax.spines['left'].set_visible(False)
                Ax.spines['bottom'].set_visible(False)
            
        self.Fig.tight_layout()
        self.Fig.subplots_adjust(hspace=0)

    def CreateFig(self, Slots, ShowLegend=False, figsize=None):
        self.ShowLegend = ShowLegend
        self.Slots = Slots

        Pos = []
        for sl in self.Slots:
            Pos.append(sl.Position)

        self.Fig, A = plt.subplots(max(Pos) + 1, 1, sharex=True, figsize=figsize)
        if type(A).__module__ == np.__name__:
            self.Axs = []
            for a in A:
                self.Axs.append([a, True])
        else:
            self.Axs = []
            self.Axs.append([A, True])

        for sl in self.Slots:  ## init slots axis
            sl.SetTstart()
            setattr(sl, 'Ax', self.Axs[sl.Position][0])

            if not ShowLegend:  # Add labels to axis
                lb = sl.Ax.get_ylabel()
                sl.Ax.set_ylabel(lb + ' ' + sl.DispName + '\n',
                                 rotation='horizontal',
                                 ha='center')
                sl.Ax.yaxis.set_label_coords(-0.1, 0.3) # TODO Check this label

            if not(sl.Ymin == 0 and sl.Ymax == 0):  # Check autoscale
                self.Axs[sl.Position][1] = False

#            if sl.PlotType == 'Spectrogram':
#                sl.Ax.set_yscale('log')
            sl.Fig = self.Fig
            if not sl.FiltType[0] == '':
                print sl.FiltType, sl.FiltF1, sl.FiltF2, sl.FiltOrder
                sl.Filter = Filter(Type=sl.FiltType,
                                   Freq1=sl.FiltF1,
                                   Freq2=sl.FiltF2,
                                   Order=sl.FiltOrder)
        self.FormatFigure()

    def PlotHist(self, Time, Resamp=False, Binds=250):
        fig, Ax = plt.subplots()

        for sl in self.Slots:
            Ax.hist(sl.GetSignal(Time, Resamp=Resamp),
                    Binds,
                    alpha=0.5)
            Ax.set_yscale('log')

    def PlotPSD(self, Time, nFFT=2**17, FMin=None, Resamp=False):

        if not self.FigFFT or not plt.fignum_exists(self.FigFFT.number):
            self.FigFFT, self.AxFFT = plt.subplots()

        PSD = {}
        for sl in self.Slots:
            sig = sl.GetSignal(Time, Resamp=Resamp)
            if FMin:
                nFFT = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/FMin))+1)) 

            ff, psd = signal.welch(x=sig, fs=sig.sampling_rate,
                                   window='hanning',
                                   nperseg=nFFT, scaling='density', axis=0)
            slN = sl.SigName
            PSD[slN] = {}
            PSD[slN]['psd'] = psd
            PSD[slN]['ff'] = ff
            self.AxFFT.loglog(ff, psd, label=sl.DispName)

        self.AxFFT.set_xlabel('Frequency [Hz]')
        self.AxFFT.set_ylabel('PSD [V^2/Hz]')
        self.AxFFT.legend()
        return PSD

    def PlotEventAvg(self, (EventRec, EventName), Time, TimeWindow,
                     OverLap=True, Std=False, Spect=False, Resamp=False):

        ft, Axt = plt.subplots()

        for sl in self.Slots:
            avg = np.array([])
            if Spect:
                ft, (Ax, AxS) = plt.subplots(2, 1, sharex=True)
            else:
                ft, Ax = plt.subplots()

            if Resamp:
                Fs = sl.ResampleFs.magnitude
            else:
                Fs = sl.Signal().sampling_rate.magnitude

            Ts = 1/Fs
            nSamps = int((TimeWindow[1]-TimeWindow[0])/Ts)
            t = np.arange(nSamps)*Ts*pq.s + TimeWindow[0]

            etimes = EventRec.GetEventTimes(EventName, Time)
            for et in etimes:
                start = et+TimeWindow[0]
                stop = et+TimeWindow[1]

                if sl.Signal().t_start < start and sl.Signal().t_stop > stop:
                    st = sl.GetSignal((start, stop), Resamp=Resamp)[:nSamps]
                    try:
                        avg = np.hstack([avg, st]) if avg.size else st
                        if OverLap:
                            Ax.plot(t, st, 'k-', alpha=0.1)
                    except:
                        print 'Error', nSamps, et, avg.shape, st.shape

            print EventName, 'Counts', len(etimes)

            MeanT = np.mean(avg, axis=1)
            Ax.plot(t, MeanT, 'r-')
            if Std:
                StdT = np.std(avg, axis=1)
                Ax.fill_between(t, MeanT+StdT, MeanT-StdT,
                                facecolor='r', alpha=0.5)

            ylim = Ax.get_ylim()
            Ax.plot((0, 0), (1, -1), 'g-.', alpha=0.7)
            Ax.set_ylim(ylim)
            Ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

            plt.title(sl.DispName)
            if Spect:
                nFFT = int(2**(np.around(np.log2(Fs/sl.SpecFmin))+1))
                noverlap = int((Ts*nFFT - sl.SpecTimeRes)/Ts)
                TWindOff = (nFFT * Ts)/8

                f, tsp, Sxx = signal.spectrogram(MeanT, Fs,
                                                 window='hanning',
                                                 nperseg=nFFT,
                                                 noverlap=noverlap,
                                                 scaling='density',
                                                 axis=0)

                finds = np.where(f < sl.SpecFmax)[0][1:]
                print Sxx.shape
                r, c = Sxx.shape
                S = Sxx.reshape((r, c))[finds][:]
                pcol = AxS.pcolormesh(tsp + TimeWindow[0].magnitude + TWindOff,
                                      f[finds],
                                      np.log10(S),
                                      vmin=np.log10(np.max(S))+sl.SpecMinPSD,
                                      vmax=np.log10(np.max(S)),
                                      cmap=sl.SpecCmap)
                f, a = plt.subplots(1, 1)
                f.colorbar(pcol)

            Axt.plot(t, np.mean(avg, axis=1), label=sl.DispName)
            ft.canvas.draw()

        Axt.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        Axt.legend()
        ft.canvas.draw()
        plt.show()

    def PlotPSDSNR(self, (EventRec, EventName), TDelay, TEval, DevDCVals,
                   Time=None, nFFT=2**17):

        etimes = EventRec.GetEventTimes(EventName, Time)

        DevACVals = {}

        for sl in self.Slots:
            fig, (ax1, ax2)  = plt.subplots(1,2)
            psd = np.array([])
            DCVals = DevDCVals[sl.SigName]
            for ne, et in enumerate(etimes):
                start = et+TDelay
                stop = et+TDelay+TEval

                sig = sl.GetSignal((start, stop), Resamp=False)
                fpsd, npsd = signal.welch(x=sig, fs=sig.sampling_rate,
                                          window='hanning',
                                          nperseg=nFFT, scaling='density', axis=0)                
                psd = np.hstack([psd, npsd]) if psd.size else npsd

#                Flin = fpsd[1:].magnitude
#                Flog = np.logspace(np.log10(Flin[0]),
#                                   np.log10(Flin[-1]), 100)
#                Flog = np.round(Flog, 9)
#                Flin = np.round(Flin, 9)
#                intpsd = interpolate.interp1d(Flin, npsd[1:].transpose())(Flog)
#                a, b, _ = FETAna.noise.FitNoise(Flog,
#                                                intpsd, Fmin=150, Fmax=5e3)

                ax1.loglog(fpsd, npsd)
#                ax2.loglog(Flog, intpsd/FETAna.noise(Flog, a, b))

            psdD = {'Vd0': psd.transpose()}
            ACVals = {'PSD': psdD,
                      'gm': None,
                      'Vgs': DCVals['Vgs'],
                      'Vds': DCVals['Vds'],
                      'Fpsd': fpsd.magnitude,
                      'Fgm': None,
                      'ChName': sl.SigName,
                      'Name': sl.DispName,
                      'GMPoly': DCVals['GMPoly'],
                      'IdsPoly': DCVals['IdsPoly'],
                      'Ud0': DCVals['Ud0'],
                      'IsOK': DCVals['IsOK'],
                      'DateTime': DCVals['DateTime']}
            DevACVals[sl.DispName] = ACVals
            fig.canvas.draw()
            plt.show()

        FETAna.InterpolatePSD(DevACVals)
        FETAna.FitACNoise(DevACVals, Fmin=150, Fmax=5e3, IsOkFilt=False)
        pltPSD = FETplt.PyFETPlot()
        pltPSD.AddAxes(('PSD', 'NoA', 'NoB'))
        pltPSD.AddLegend()
        pltPSD.PlotDataCh(DevACVals)

        pickle.dump(DevACVals, open('test.pkl', 'wb'))

        return DevACVals

    def AddLegend(self, Ax):
        nLines = len(Ax.lines)
        nlc = self.LegNlabCol
        if nLines > nlc:
            ncol = (nLines / nlc) + ((nLines % nlc) > 0)
        else:
            ncol = 1
        Ax.legend(loc='best',
                  fancybox=True,
                  shadow=True,
                  ncol=ncol,
                  fontsize=self.LegFontSize)

    def PlotChannels(self, Time, Resamp=True,
                     ResampPoints=None, ResampFs=None):

        if not self.Fig:
            return

        for sl in self.Slots:

            if ResampPoints:
                sl.ResamplePoints = ResampPoints
            if ResampFs:
                sl.ResampFs = ResampFs

            sl.PlotSignal(Time, Resamp=Resamp)

        if Time is not None:
            sl.Ax.set_xlim(Time[0], Time[1])

        if self.ShowLegend:
            for (Ax, AutoScale) in self.Axs:
                self.AddLegend(Ax)
                if AutoScale:
                    Ax.autoscale(enable=True, axis='y', tight=True)

    def PlotEvents(self, Time, (EventRec, EventName), color=None, alpha=0.5):

        Te = EventRec.GetEventTimes(EventName, Time)
        for te in Te:
            for sl in self.Slots:
                sl.PlotEvent(te, color=color, alpha=alpha)
