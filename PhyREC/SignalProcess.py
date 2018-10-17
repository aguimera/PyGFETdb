#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:27:23 2018

@author: aguimera
"""
import numpy as np
from scipy import signal
from fractions import Fraction


def DownSampling(sig, Fact, zero_phase=True):
    print sig.sampling_rate, sig.sampling_rate/Fact
    rs = signal.decimate(np.array(sig),
                         q=Fact,
                         zero_phase=zero_phase,
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
    print sig.name, offset
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
        print sig.sampling_rate*f, f, uprate, dowrate
        rs = signal.resample_poly(np.array(sig), uprate, dowrate)
        sig = sig.duplicate_with_new_array(signal=rs*sig.units)
        sig.sampling_rate = sig.sampling_rate*f
        return sig
    else:
        return sig


def Filter(sig, Type, Order, Freqs):
    st = np.array(sig)
    freqs = Freqs/(0.5*sig.sampling_rate)

    b, a = signal.butter(Order, freqs, Type)
    st = signal.filtfilt(b, a, st, axis=0)

    return sig.duplicate_with_new_array(signal=st*sig.units)

def sliding_window(sig,func,timewidth,step):
    #func can be average/std or others to be added
    window_size = int(timewidth*sig.sampling_rate.item())
    stride = int(step*sig.sampling_rate.item())
    if func=='std':
        window_res = [ np.std(sig[i:i+window_size]) for i in range(0, len(sig), stride)
                       if i+window_size <= len(sig) ]
    elif func=='avg':    
         window_res = [ np.mean(sig[i:i+window_size]) for i in range(0, len(sig), stride)
                       if i+window_size <= len(sig) ]
         
    return window_res,stride
