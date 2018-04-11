#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 16:37:02 2018

@author: aguimera
"""

from PyREC.NeoInterface import NeoSegment, NeoSignal
import PyREC.Plots as PltRec
import PyREC.SignalProcess as RPro
import quantities as pq
import matplotlib.pyplot as plt
import neo 
from scipy import signal


plt.close('all')

Rec = NeoSegment('Exp1_full.h5')

TWind = (500*pq.s, 2500*pq.s)

ProcessChain = [
                {'function': RPro.RemoveDC, 'args': {}},
                {'function': RPro.Filter, 'args': {'Type':'highpass',
                                                   'Order':2,
                                                   'Freqs':(1,)}},
                {'function': RPro.Filter, 'args': {'Type':'bandstop',
                                                   'Order':2,
                                                   'Freqs':(48, 52)}},
#                {'function': Resample, 'args': {'Fs':50}}
                ]

MySigs = []
for sn in Rec.SigNames.keys():
    MySig = RPro.SignalProcess(Signal=NeoSignal(Rec, sn),
                               ProcessChain=ProcessChain)    
#    MySig=NeoSignal(Rec, sn)
    MySigs.append(MySig)

Slots = []
for MySig in MySigs:
    if MySig.Name not in ('Au-Ch01','Pt-Ch06','Vge-Ch02'):
        continue
    sl = PltRec.WaveSlot(MySig)
    if MySig.Name.startswith('Au'):
        sl.Color = 'b'
        sl.Position = 0
    if MySig.Name.startswith('Pt'):
        sl.Color = 'k'
        sl.Position = 1
    if MySig.Name.startswith('Vge'):
        sl.Color = 'r'
        sl.Position = 2
#    sl.Alpha = 0.2
#    sl.Ymax = 1
#    sl.Ymin = -1
#    sl.AutoScale = False
    Slots.append(sl)

PlotRecs = PltRec.PlotSlots(Slots)
PlotRecs.PlotChannels(TWind, Units='mV')

#sig = Rec.GetSignal('Au-Ch01')
#plt.plot(sig.times, sig)
#st = signal.detrend(sig, axis=0)
#plt.plot(sig.times, st)

