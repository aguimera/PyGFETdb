#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:27:23 2018

@author: aguimera
"""
import numpy as np
from scipy import signal
from fractions import Fraction
from PhyREC.NeoInterface import NeoSegment, NeoSignal, NeoTrain
import PhyREC.SignalAnalysis as Ran
import quantities as pq
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import elephant
import scipy.stats as stats
from scipy import interpolate
import sys
from scipy.stats.mstats import zscore
from . import DbgFplt
from scipy.interpolate import UnivariateSpline
from scipy.signal import medfilt


def Spectrogram(sig, Fres=2*pq.Hz, TimeRes=0.01*pq.s,
                Fmin=1*pq.Hz, Fmax=200*pq.Hz, Zscored=True, NormTime=None,
                dtype=np.float,
                **specKwarg):

    nFFT = int(2**(np.around(np.log2(sig.sampling_rate/Fres))+1))
    Ts = sig.sampling_period
    noverlap = int((Ts*nFFT - TimeRes)/Ts)

    f, t, Sxx = signal.spectrogram(sig,
                                   fs=sig.sampling_rate,
                                   nperseg=nFFT,
                                   noverlap=noverlap,
                                   axis=0,
                                   **specKwarg)

    finds = np.where((Fmin < f) & (f < Fmax))[0][1:]
    r, g, c = Sxx.shape
    data = Sxx.reshape((r, c))[finds][:]

    if Zscored and (NormTime is None):
        data = zscore(data, axis=1)

    s = sig.duplicate_with_new_data(data.astype(dtype).transpose())
    s.annotate(Freq=f[finds])
    s.annotate(spec=True)
    s.annotate(nFFT=nFFT)
    s.annotate(WindowTime=nFFT*Ts)
    s.sampling_period = np.mean(t[1:]-t[:-1])*pq.s
    s.t_start = s.t_start + (nFFT*Ts)/2

    if Zscored and (NormTime is None):
        return s

    if Zscored:
        NormSig = s.time_slice(NormTime[0], NormTime[1])
        mean = np.mean(NormSig, axis=0)
        std = np.std(NormSig, axis=0)
        return ((s - mean) / std)
    else:
        if NormTime is None:
            return s
        else:
            NormSig = s.time_slice(NormTime[0], NormTime[1])
            mean = np.mean(NormSig, axis=0)           
            return (s / mean)


def AvgSpectrogram(sig, TimesEvent, TimeAvg, SpecArgs,
                   AvgSpectNorm='Zscore', AvgSpectNormTime=None,
                   TrialProcessChain=None, **kwargs):

    Acc = np.array([])
    Trials = 0
    for et in TimesEvent:
        if et+TimeAvg[0] < sig.t_start:
            print(et, ' not valid')
            continue
        elif et+TimeAvg[1] > sig.t_stop:
            print(et, ' not valid')
            continue
        
        Trials += 1
        s = sig.time_slice(et+TimeAvg[0], et+TimeAvg[1])
        if TrialProcessChain is not None:
            st = ApplyProcessChain(s, TrialProcessChain)
        else:
            st = s  
      
        spect = Spectrogram(st, **SpecArgs)
        Acc = Acc + np.array(spect) if Acc.size else np.array(spect)

    AvgSpect = spect.duplicate_with_new_data(Acc / Trials)
    AvgSpect.t_start = TimeAvg[0] + AvgSpect.annotations['WindowTime']/2

    if AvgSpectNorm is None:
        return AvgSpect

    if AvgSpectNormTime is None:
        NormTime = (AvgSpect.t_start, -0*pq.s)
    else:
        NormTime = AvgSpectNormTime
        if NormTime[0] is None:
            NormTime[0] = AvgSpect.t_start
        if NormTime[1] is None:
            NormTime[1] = AvgSpect.t_stop

    NormSig = AvgSpect.time_slice(NormTime[0], NormTime[1])
    mean = np.mean(NormSig, axis=0)
    std = np.std(NormSig, axis=0)
    
    if AvgSpectNorm == 'Zscore':
        AvgNormSpect = (AvgSpect - mean) / std
    else:
        AvgNormSpect = AvgSpect / mean
    
    return AvgNormSpect


def TrigAveraging(sig, TimesEvent, TimeAvg, TrialProcessChain=None):
    Ts = sig.sampling_period
    nSamps = int((TimeAvg[1] - TimeAvg[0])/Ts)
    acc = None
    for et in TimesEvent:
        if et+TimeAvg[0] < sig.t_start:
            print(et, ' not valid')
            continue
        elif et+TimeAvg[1] > sig.t_stop:
            print(et, ' not valid')
            continue
    
        Samp1 = sig.time_index(et+TimeAvg[0])
        s = sig[Samp1:Samp1+nSamps, :]    
    
        if TrialProcessChain is not None:
            st = ApplyProcessChain(s, TrialProcessChain)
        else:
            st = s
    
        st.t_start = TimeAvg[0]
        st.array_annotations = {}
        st.annotations = {}
        if acc is None:
            acc = st
        else:
            acc = acc.merge(st)

    avg = acc.duplicate_with_new_data(np.mean(acc, axis=1))
    std = acc.duplicate_with_new_data(np.std(acc, axis=1))
    avg.annotate(std=std)
    avg.annotate(acc=acc)
    return avg

def Derivative(sig):
    derivative_sig = NeoSignal(
            np.diff(sig.as_quantity(), axis=0) / sig.sampling_period,
            t_start=sig.t_start+sig.sampling_period/2,
            sampling_period=sig.sampling_period,
            name=sig.name)

    return derivative_sig


def DownSampling(sig, Fact, zero_phase=True):
    print(sig.sampling_rate, sig.sampling_rate/Fact)
    rs = signal.decimate(np.array(sig),
                         q=Fact,
                         zero_phase=zero_phase,
                         ftype='iir',
                         axis=0)
    sig = sig.duplicate_with_new_data(signal=rs*sig.units)
    sig.sampling_rate = sig.sampling_rate/Fact
    return sig


def RemoveDC(sig, Type='constant'):
    st = np.array(sig)
    st = signal.detrend(st, type=Type, axis=0)
    return sig.duplicate_with_new_data(signal=st*sig.units)


def SetZero(sig, TWind=None):
    if TWind is None:
        TWind=(sig.t_start, sig.t_start+30*pq.s)
    st = np.array(sig)
    # offset = np.mean(sig.GetSignal(TWind))
    offset = np.mean(sig.time_slice(TWind[0], TWind[1]))
    print(sig.name, offset)
    st_corrected = st-offset.magnitude
    return sig.duplicate_with_new_data(signal=st_corrected*sig.units)


def Gain(sig, Gain):
    return sig*Gain


def Resample(sig, Fs=None, MaxPoints=None):
    if MaxPoints is None:
        f = Fs/sig.sampling_rate
        fact = Fraction(float(f)).limit_denominator()
        dowrate = fact.denominator
        uprate = fact.numerator
    else:
        dowrate = int(sig.times.shape[0]/MaxPoints)
        if dowrate > 0:
            f = float(1/float(dowrate))
            uprate = 1

    if dowrate > 0:
        print(sig.sampling_rate*f, f, uprate, dowrate)
        rs = signal.resample_poly(np.array(sig), uprate, dowrate)
        sig = sig.duplicate_with_new_data(signal=rs*sig.units)
        sig.sampling_rate = sig.sampling_rate*f
        return sig
    else:
        return sig


def Abs(sig):
    st = np.array(sig)
    st = np.abs(st)

    return sig.duplicate_with_new_data(signal=st*sig.units)


def power(sig):  # to solve units
    st = np.array(sig)**2
#    st = st**2
    return sig.duplicate_with_new_data(signal=st*sig.units)


def Filter(sig, Type, Order, Freqs):
    st = np.array(sig)
    Fs = sig.sampling_rate
    freqs = Freqs/(0.5 * Fs)

    # b, a = signal.butter(Order, freqs, Type)
    # st = signal.filtfilt(b, a, st, axis=0)
    
    sos = signal.butter(Order, freqs, Type, output='sos')
    st = signal.sosfiltfilt(sos, st, axis=0)
    
    # DbgFplt.PlotResponse(a, b, Fs)
    DbgFplt.PlotResponse(sos, Fs)
    
    return sig.duplicate_with_new_data(signal=st*sig.units)


def rms(x, axis=None):
    return np.sqrt(np.mean(x**2, axis=axis))


def power_sliding(x, axis=None):
    return np.mean(x**2, axis=axis)



def sliding_window(sig, timewidth, func=None, steptime=None,**kwargs):

    if steptime is None:
        steptime = timewidth/10

    if func is None:
        func = rms

    window_size = int(timewidth.rescale('s') / sig.sampling_period.rescale('s'))
    timewidth = sig.sampling_period.rescale('s')*window_size
    step_size = int(steptime.rescale('s') / sig.sampling_period.rescale('s'))
    steptime = sig.sampling_period.rescale('s')*step_size

    axis = 0
    shape = list(sig.shape)
    shape[axis] = np.floor(sig.shape[axis] / step_size - window_size / step_size + 1).astype(int)
    shape.append(window_size)

    strides = list(sig.strides)
    strides[axis] *= step_size
    strides.append(sig.strides[axis])

    strided = np.lib.stride_tricks.as_strided(sig, shape=shape, strides=strides)

    st = func(strided, axis=-1,**kwargs)
 
    return NeoSignal(signal=st,
                     units=sig.units,
                     t_start=sig.t_start + timewidth/2,
                     name=sig.name,
                     sampling_rate=1/steptime)


def ThresholdTrianGen(sig, RelaxTime=0.4*pq.s, threshold=None):

    if threshold is None:
        threshold = np.mean(sig) + np.std(sig)
    inttimes = Ran.threshold_detection(signal=sig,
                                       threshold=threshold,
                                       RelaxTime=RelaxTime)
    inttimes = np.array(inttimes)

    return NeoTrain(times=inttimes,
                    units='s',
                    t_start=sig.t_start,
                    t_stop=sig.t_stop,
                     )


def ThresholdInstantRate(sig, RelaxTime=0.1*pq.s, threshold=None,
                         OutSampling=0.01*pq.s,):

    

    return elephant.statistics.instantaneous_rate(ThresholdTrianGen(sig,
                                                                    RelaxTime,
                                                                    threshold),
                                                  sampling_period=OutSampling)


def HilbertInstantFreq(sig, MaxFreq=20, MinFreq=0):
    SigH = elephant.signal_processing.hilbert(sig)
    insfreq = np.diff(np.angle(SigH)[:, 0]) / np.diff(SigH.times)

    return NeoSignal(signal=np.clip(insfreq.magnitude, 0, 20),
                     units='Hz',
                     name=sig.name,
                     sampling_rate=SigH.sampling_rate,
                     t_start=SigH.t_start)

def HilbertAngle(sig):
    SigH = elephant.signal_processing.hilbert(sig)   
    return sig.duplicate_with_new_data(signal=np.angle(SigH),
                                       units=pq.radians)


def HilbertAmp(sig):
    SigH = elephant.signal_processing.hilbert(sig)
    
    return sig.duplicate_with_new_data(signal=np.array(np.abs(SigH))*sig.units)


def SplineSmooth(sig, sFact=2, **kwargs):
    s = sig.shape[0]/sFact
    spl = UnivariateSpline(sig.times, sig, s=s)    
    return sig.duplicate_with_new_data(signal=spl(sig.times)*sig.units)

def MedianFilt(sig, window_size=None, **kwargs): #window_size in pq.s 
    if window_size==None:
        kernel_size=int(sig.shape[0]/10)
    else:
        kernel_size=int(sig.sampling_rate*window_size)
    if bool(kernel_size & 1): ##kernel_size has to be an odd number
        kernel_size=kernel_size
    else:
        kernel_size=kernel_size-1
    SigFilt=medfilt(np.array(sig).reshape(len(sig),), kernel_size=kernel_size)    # TODO change by axis
    return sig.duplicate_with_new_data(signal=np.array(SigFilt)*sig.units)

   
def Strid_signal(sig, timewidth, steptime=None):
    if steptime is None:
        steptime = timewidth/10

    window_size = int(timewidth.rescale('s') / sig.sampling_period.rescale('s'))
    timewidth = sig.sampling_period.rescale('s')*window_size
    step_size = int(steptime.rescale('s') / sig.sampling_period.rescale('s'))
    steptime = sig.sampling_period.rescale('s')*step_size

    axis = 0
    shape = list(sig.shape)
    shape[axis] = np.floor(sig.shape[axis] / step_size - window_size / step_size + 1).astype(int)
    shape.append(window_size)

    strides = list(sig.strides)
    strides[axis] *= step_size
    strides.append(sig.strides[axis])

    strided = np.lib.stride_tricks.as_strided(sig, shape=shape, strides=strides)

    return strided

def CrossCorr(x1, x2):
    max_cross_corr=[]
    for row in range(x1.shape[0]):   
#        print(x1.shape)
#        print (x1[row,:,:].shape)
#        print(x2.shape)
#        print (x2[row,:,:].shape)
#        print(row)
        res=xcorr(
                        np.array(x1[row,:,:]).reshape((x1[row,:,:].shape[-1],)),
                      np.array(x2[row,:,:]).reshape((x2[row,:,:].shape[-1],)),maxlags=None)
        cross_correlation=max(res[1])
        max_cross_corr.append(cross_correlation)
        
    return max_cross_corr


def PearsonCorr(x1, x2):
    pearson_corr=[]
    for row in range(x1.shape[0]):   
        res=stats.pearsonr(
                        np.array(x1[row,:,:]).reshape((x1[row,:,:].shape[-1],)),
                      np.array(x2[row,:,:]).reshape((x2[row,:,:].shape[-1],)))
        pearson_corr.append(res[0])
        
    return pearson_corr


def sliding_window_2sigs(sig1, sig2, timewidth, func=None, steptime=None):

    if func is None:
        func = CrossCorr

    strided1 = Strid_signal(sig1, timewidth, steptime)
    strided2 = Strid_signal(sig2, timewidth, steptime)

    st = func(strided1, strided2)

    return NeoSignal(signal=st,
                     units=sig1.units,
                     t_start=sig1.t_start + timewidth/2,
                     name='CrossCorr'+sig1.name+sig1.name,
                     sampling_rate=1/steptime)    
    
    
def xcorr(x, y, normed=True, detrend=mlab.detrend_none,
              usevlines=True, maxlags=10, **kwargs):
        r"""
        Plot the cross correlation between *x* and *y*.

        The correlation with lag k is defined as
        :math:`\sum_n x[n+k] \cdot y^*[n]`, where :math:`y^*` is the complex
        conjugate of :math:`y`.

        Parameters
        ----------
        x : array-like of length n

        y : array-like of length n

        detrend : callable, optional, default: `mlab.detrend_none`
            *x* and *y* are detrended by the *detrend* callable. This must be a
            function ``x = detrend(x)`` accepting and returning an
            `numpy.array`. Default is no normalization.

        normed : bool, optional, default: True
            If ``True``, input vectors are normalised to unit length.

        usevlines : bool, optional, default: True
            Determines the plot style.

            If ``True``, vertical lines are plotted from 0 to the xcorr value
            using `Axes.vlines`. Additionally, a horizontal line is plotted
            at y=0 using `Axes.axhline`.

            If ``False``, markers are plotted at the xcorr values using
            `Axes.plot`.

        maxlags : int, optional, default: 10
            Number of lags to show. If None, will return all ``2 * len(x) - 1``
            lags.

        Returns
        -------
        lags : array (length ``2*maxlags+1``)
            The lag vector.
        c : array  (length ``2*maxlags+1``)
            The auto correlation vector.
        line : `.LineCollection` or `.Line2D`
            `.Artist` added to the axes of the correlation:

            - `.LineCollection` if *usevlines* is True.
            - `.Line2D` if *usevlines* is False.
        b : `.Line2D` or None
            Horizontal line at 0 if *usevlines* is True
            None *usevlines* is False.

        Other Parameters
        ----------------
        linestyle : `.Line2D` property, optional
            The linestyle for plotting the data points.
            Only used if *usevlines* is ``False``.

        marker : str, optional, default: 'o'
            The marker for plotting the data points.
            Only used if *usevlines* is ``False``.

        Notes
        -----
        The cross correlation is performed with :func:`numpy.correlate` with
        ``mode = "full"``.
        """
        Nx = len(x)
        if Nx != len(y):
            raise ValueError('x and y must be equal length')

        x = detrend(np.asarray(x))
        y = detrend(np.asarray(y))

        correls = np.correlate(x, y, mode="full")

        if normed:
            correls /= np.sqrt(np.dot(x, x) * np.dot(y, y))

        if maxlags is None:
            maxlags = Nx - 1

        if maxlags >= Nx or maxlags < 1:
            raise ValueError('maxlags must be None or strictly '
                             'positive < %d' % Nx)

        lags = np.arange(-maxlags, maxlags + 1)
        correls = correls[Nx - 1 - maxlags:Nx + maxlags]

#        if usevlines:
#            a = self.vlines(lags, [0], correls, **kwargs)
#            # Make label empty so only vertical lines get a legend entry
#            kwargs.pop('label', '')
#            b = self.axhline(**kwargs)
#        else:
#            kwargs.setdefault('marker', 'o')
#            kwargs.setdefault('linestyle', 'None')
#            a, = self.plot(lags, correls, **kwargs)
#            b = None
        return lags, correls


def CalcVgeff(Sig, Tchar, VgsExp=None, Regim='hole', CalType='interp'):
    Vgs = Tchar.GetVgs()
    vgs = np.linspace(np.min(Vgs), np.max(Vgs), 10000)

    if Regim == 'hole':
        Inds = np.where(vgs < Tchar.GetUd0())[1]
    else:
        Inds = np.where(vgs > Tchar.GetUd0())[1]

    Ids = Tchar.GetIds(Vgs=vgs[Inds]) * pq.A
    GM = Tchar.GetGM(Vgs=VgsExp) * pq.S


    IdsExp = Tchar.GetIds(Vgs=VgsExp) * pq.A
    IdsOff = np.mean(Sig)-IdsExp
    IdsBias = np.mean(Sig)

    Calibrated = np.array((True,))
    try:
        if CalType == 'interp':
            fgm = interpolate.interp1d(Ids[:, 0], vgs[Inds])
            st = fgm(np.clip(Sig, np.min(Ids), np.max(Ids)))*pq.V
        elif CalType=='linear':
            st = Sig/GM
        else:
            print('Calibration Not defined')
    except:
        print(Sig.name, "Calibration error:", sys.exc_info()[0])
        st = np.zeros(Sig.shape)
        Calibrated = np.array((False,))

    print(str(Sig.name), '-> ',
          'IdsBias', IdsBias,
          'IdsOff', IdsOff,
          'Vgs', np.mean(st),
          Tchar.IsOK)

    annotations = {'Calibrated': Calibrated,
                   'Working': Calibrated,
                   'IdsOff': IdsOff.flatten()[0],
                   'VgsCal': np.mean(st),
                   'IdsBias': IdsBias.flatten()[0],
                   'IsOK': Tchar.IsOK,
                   'Iname': Sig.name,
                   'GM': GM,
                   }

    CalSig = NeoSignal(st,
                       units='V',
                       t_start=Sig.t_start,
                       sampling_rate=Sig.sampling_rate,
                       name=str(Sig.name),
                       file_origin=Sig.file_origin)

#    CalSig.annotate(**annotations)
    CalSig.array_annotate(**annotations)

    return CalSig


def CalcVgeff2(Sig, Tchar, VgsExp, Regim='hole', CalType='interp'):
    Vgs = Tchar.GetVgs()
    vgs = np.linspace(np.min(Vgs), np.max(Vgs), 10000)

    if Regim == 'hole':
        Inds = np.where(vgs < Tchar.GetUd0())[1]
    else:
        Inds = np.where(vgs > Tchar.GetUd0())[1]

    IdsBias = np.array((np.nan,))*pq.A
    IdsOff = np.array((np.nan,))*pq.A
    GM = np.nan*pq.S
    try:    
        Ids = Tchar.GetIds(Vgs=vgs[Inds]).flatten() 
        GM = Tchar.GetGM(Vgs=VgsExp).flatten() 


        IdsExp = Tchar.GetIds(Vgs=VgsExp).flatten() 
        IdsOff = np.mean(Sig)-IdsExp
        IdsBias = np.mean(Sig)
    
        Calibrated = np.array((True,))
    
        if CalType == 'interp':
            fgm = interpolate.interp1d(Ids, vgs[Inds])
            st = fgm(np.clip(Sig, np.min(Ids), np.max(Ids)))*pq.V
        elif CalType=='linear':
            st = Sig/GM
        else:
            print('Calibration Not defined')
    except:
        print(Sig.name, "Calibration error:", sys.exc_info()[0])
        st = np.zeros(Sig.shape)*pq.V
        Calibrated = np.array((False,))

    print(str(Sig.name), '-> ',
          'IdsBias', IdsBias,
          'IdsOff', IdsOff,
          'Vgs', np.mean(st),
          'GM', GM,
          Tchar.IsOK)

    annotations = {'Calibrated': Calibrated[0],
                   'Working': Calibrated[0],
                   'IdsOff': float(IdsOff.flatten()[0].magnitude),
                   'VgsCal': float(np.mean(st).magnitude),
                   'IdsBias': float(IdsBias.flatten()[0].magnitude),
                   'IsOK': Tchar.IsOK,
                   'Iname': Sig.name,
                   'GM': float(GM.magnitude),
                   }

    CalSig = NeoSignal(st,
                       units='V',
                       t_start=Sig.t_start,
                       sampling_rate=Sig.sampling_rate,
                       name=str(Sig.name),
                       file_origin=Sig.file_origin)

    CalSig.annotate(**annotations)
    # CalSig.array_aºnnotate(**annotations)

    return CalSig

def CalcVgeffNoInterp(Sig, Tchar, VgsExp=None, Regim='hole'):    
    gm=Tchar.GetGM(Vgs=VgsExp)


    Calibrated = np.array((True,))
    try:
        st=Sig.magnitude/gm
    except:
        print(Sig.name, "Calibration error:", sys.exc_info()[0])
        st = np.zeros(Sig.shape)
        Calibrated = np.array((False,))
        
        
    print(str(Sig.name), '-> ', 'GM', gm, 'Vgs', np.mean(st), Tchar.IsOK)
    annotations = {'Calibrated': Calibrated,
                   'Working':Calibrated,
                   # 'IdsOff': IdsOff.flatten(),
                   'VgsCal': np.array((np.mean(st), )),
                   'IsOK': np.array((Tchar.IsOK, )),
                   'Iname': np.array((Sig.name, )),                   
                   }
    
    CalSig = NeoSignal(st*pq.V,
                       units='V',
                       t_start=Sig.t_start,
                       sampling_rate=Sig.sampling_rate,
                       name=str(Sig.name),
                       file_origin=Sig.file_origin)

#    CalSig.annotate(**annotations)    
    CalSig.array_annotate(**annotations)
        
    return CalSig



def ApplyProcessChain(sig, ProcessChain):
    sl = sig.copy()  
    for Proc in ProcessChain:
        sl = Proc['function'](sl, **Proc['args'])
    
    return sl



