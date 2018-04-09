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


class NeoRecord():
    def __init__(self, RecordFile=None, UnitGain=1e6, Seg=None):

        self.UnitGain = UnitGain

        if not RecordFile:
            self.Seg = Seg
            if self.Seg:
                self.UpdateEventDict()
                self.UpdateSignalsDict()
            else:
                self.Seg = neo.Segment('New Seg')
            return

        ftype = RecordFile.split('.')[-1]
        if ftype == 'h5':
            self.RecFile = neo.io.NixIO(filename=RecordFile, mode='ro')
            Block = self.RecFile.read_block()
        elif ftype == 'smr':
            self.RecFile = neo.io.Spike2IO(filename=RecordFile)
            Block = self.RecFile.read()[0]

        self.Seg = Block.segments[0]

        self.UpdateSignalsDict()
        self.UpdateEventDict()

    def SaveRecord(self, FileName):
        out_f = neo.io.NixIO(filename=FileName)
        out_bl = neo.Block(name='NewBlock')

        out_bl.segments.append(self.Seg)
        out_f.write_block(out_bl)
        out_f.close()

    def UpdateSignalsDict(self):
        self.SigNames = {}
        for i, sig in enumerate(self.Seg.analogsignals):
            if sig.name is None:
                name = str(i)
                if sig.annotations is not None:
                    if 'nix_name' in sig.annotations.keys():
                        name = sig.annotations['nix_name']
            else:
                name = str(sig.name)

            self.SigNames.update({name: i})

    def UpdateEventDict(self):
        self.EventNames = {}
        for i, eve in enumerate(self.Seg.events):
            if eve.name is None:
                try:                    
                    name = eve.annotations['title']
                except:
                    print 'Event found no name ', i
                    name = str(i)
                self.EventNames.update({name: i})
            else:
                self.EventNames.update({eve.name: i})

    def GetEventTimes(self, EventName, Time=None):
        eve = self.Seg.events[self.EventNames[EventName]].times
        if Time:
            events = eve[np.where((eve > Time[0]) & (eve < Time[1]))]
        else:
            events = eve
        return events

    def GetTstart(self, ChName):
        return self.Seg.analogsignals[self.SigNames[ChName]].t_start

    def SetTstart(self, ChName, Tstart):
        self.Seg.analogsignals[self.SigNames[ChName]].t_start = Tstart

    def SetSignal(self, ChName, Sig):
        self.Seg.analogsignals[self.SigNames[ChName]] = Sig

    def AppendSignal(self, ChName, Vect):
        sig = self.Signal(ChName,
                          Scale=False)
        S_old = sig.copy()
        v_old = np.array(sig)
        v_new = np.vstack((v_old, Vect))
        sig_new = neo.AnalogSignal(v_new,
                                   units=S_old.units,
                                   sampling_rate=S_old.sampling_rate,
                                   t_start=S_old.t_start,
                                   name=S_old.name)
        self.SetSignal(ChName, sig_new)

    def Signal(self, ChName, Gain=None, Scale=True):
        sig = self.Seg.analogsignals[self.SigNames[ChName]]

        if Scale:
            if Gain:
                return sig/Gain
            else:
                return sig/self.UnitGain
        else:
            return sig

    def GetSignal(self, ChName, Time, Gain=None):
        sig = self.Seg.analogsignals[self.SigNames[ChName]]
        sl = sig.time_slice(Time[0], Time[1])

        if Gain:
            return sl/Gain
        else:
            return sl/self.UnitGain

    def AddEvent(self, Times, Name):
        eve = neo.Event(times=Times,
                        units=pq.s,
                        name=Name)

        self.Seg.events.append(eve)
        self.UpdateEventDict()

    def AddSignal(self, Sig):
        self.Seg.analogsignals.append(Sig)
        self.UpdateSignalsDict()

#    def AddSignal(self, ChName, Vect, Var): #  TODO check var???
#        sig = self.Signal(ChName,
#                          Scale=False)
#        S_old = sig.copy()
#        sig_new = neo.AnalogSignal(Vect,
#                                   units=S_old.units,
#                                   sampling_rate=S_old.sampling_rate,
#                                   t_start=S_old.t_start,
#                                   name=S_old.name+Var)
#        self.Seg.analogsignals.append(sig_new)
#        self.UpdateSignalsDict()


