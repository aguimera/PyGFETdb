#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:27:23 2018

@author: aguimera
"""
import numpy as np
from scipy import signal
from fractions import Fraction


def RemoveDC(sig):
    st = np.array(sig)
    st = signal.detrend(st, axis=0)
    return sig.duplicate_with_new_array(signal=st*sig.units)


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


class SignalProcess():
    def __init__(self, Signal, ProcessChain):
        self.Signal = Signal
        self.signal = Signal.signal
        self.ProcessChain = ProcessChain
        self.Name = Signal.Name

    def GetSignal(self, Time, Units=None):
        sig = self.Signal.GetSignal(Time, Units)
        for Proc in self.ProcessChain:
            sig = Proc['function'](sig, **Proc['args'])
        return sig
