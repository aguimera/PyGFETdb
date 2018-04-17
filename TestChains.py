#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 16:37:02 2018

@author: aguimera
"""

from PyREC.NeoInterface import NeoSegment, NeoSignal
import PyREC.PlotWaves as Rplt
import PyREC.SignalProcess as RPro
import PyREC.SignalAnalysis as Ran
import quantities as pq
import matplotlib.pyplot as plt
import neo 
from scipy import signal
import gc

plt.close('all')
print 'Collect', gc.collect()

Rec = NeoSegment('Exp1_full.h5')

TWind = (1500*pq.s, 2500*pq.s)

ProcessChain = [
                {'function': RPro.RemoveDC, 'args': {}},
                {'function': RPro.Filter, 'args': {'Type':'highpass',
                                                   'Order':2,
                                                   'Freqs':(20,)}},
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
#    if MySig.Name not in ('Au-Ch01','Pt-Ch06','Vge-Ch02'):
#        continue
    sl = Rplt.WaveSlot(MySig)
    if MySig.Name.startswith('Au'):
        sl.Color = None
#        sl.Color = 'b'
#        sl.DispName = 'Au'
        sl.Position = 0        
    if MySig.Name.startswith('Pt'):
        sl.Color = 'k'
        sl.DispName = 'Pt'
        sl.Position = 1        
    if MySig.Name.startswith('Vge'):
        sl.Color = 'r'
        sl.DispName = 'SGFETs'
        sl.Position = 2
#        sls = Rplt.SpecSlot(MySig)
#        sls.Position = 3
#        sls.Fmin = 5
#        sls.Fmax = 100
#        Slots.append(sls)
            
    sl.Alpha = 0.2
#    sl.Ymax = 1
#    sl.Ymin = -1
#    sl.AutoScale = False
    Slots.append(sl)

PlotRecs = Rplt.PlotSlots(Slots, ShowNameOn='Legend', ShowAxis=True)
PlotRecs.PlotChannels(TWind, Units='V')

PSD = Ran.PlotPSD(Slots, Units='uV', FMin=0.1, scaling='spectrum')

#sig = Rec.GetSignal('Au-Ch01')
#plt.plot(sig.times, sig)
#st = signal.detrend(sig, axis=0)
#plt.plot(sig.times, st)

print 'Collect', gc.collect()