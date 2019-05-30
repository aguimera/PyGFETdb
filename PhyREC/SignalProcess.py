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
import elephant

def DownSampling(sig, Fact, zero_phase=True):
    print (sig.sampling_rate, sig.sampling_rate/Fact)
    rs = signal.decimate(np.array(sig),
                         q=Fact,
                         zero_phase=zero_phase,
                         ftype='fir',
                         axis=0)
    sig = sig.duplicate_with_new_array(signal=rs*sig.units)
    sig.sampling_rate = sig.sampling_rate/Fact
    return sig


def RemoveDC(sig, Type='constant'):
    st = np.array(sig)
    st = signal.detrend(st, type=Type, axis=0)
    return sig.duplicate_with_new_array(signal=st*sig.units)


def SetZero(sig, TWind):
    st = np.array(sig)
    offset = np.mean(sig.GetSignal(TWind))
    print (sig.name, offset)
    st_corrected = st-offset.magnitude
    return sig.duplicate_with_new_array(signal=st_corrected*sig.units)


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
        print (sig.sampling_rate*f, f, uprate, dowrate)
        rs = signal.resample_poly(np.array(sig), uprate, dowrate)
        sig = sig.duplicate_with_new_array(signal=rs*sig.units)
        sig.sampling_rate = sig.sampling_rate*f
        return sig
    else:
        return sig


def Abs(sig):
    st = np.array(sig)
    st = np.abs(st)

    return sig.duplicate_with_new_array(signal=st*sig.units)


def Filter(sig, Type, Order, Freqs):
    st = np.array(sig)
    freqs = Freqs/(0.5*sig.sampling_rate)

    b, a = signal.butter(Order, freqs, Type)
    st = signal.filtfilt(b, a, st, axis=0)

    return sig.duplicate_with_new_array(signal=st*sig.units)


def rms(x, axis=None):
    return np.sqrt(np.mean(x**2, axis=axis))


def sliding_window(sig, timewidth, func=None, steptime=None):

    if steptime is None:
        steptime = timewidth/10

    if func is None:
        func = rms

    window_size = int(timewidth.rescale('s') / sig.sampling_period.rescale('s'))
    step_size = int(steptime.rescale('s') / sig.sampling_period.rescale('s'))

    axis = 0
    shape = list(sig.shape)
    shape[axis] = np.floor(sig.shape[axis] / step_size - window_size / step_size + 1).astype(int)
    shape.append(window_size)

    strides = list(sig.strides)
    strides[axis] *= step_size
    strides.append(sig.strides[axis])

    strided = np.lib.stride_tricks.as_strided(sig, shape=shape, strides=strides)

    st = func(strided, axis=-1)
 
    return NeoSignal(signal=st,
                     units=sig.units,
                     t_start=sig.t_start,
                     name=sig.name,
                     sampling_rate=1/steptime)


def HilbertInstantFreq(sig, MaxFreq=20, MinFreq=0):
    SigH = elephant.signal_processing.hilbert(sig)
    insfreq = np.diff(np.angle(SigH)[:, 0]) / np.diff(SigH.times)

    return NeoSignal(signal=np.clip(insfreq.magnitude, 0, 20),
                     units='Hz',
                     name=sig.name,
                     sampling_rate=SigH.sampling_rate,
                     t_start=SigH.t_start)
    
    
    
    
    
    
    
    


