#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:11:53 2017

@author: aguimera
"""

import neo
import numpy as np
import quantities as pq


class NeoSegment():
    def __init__(self, RecordFile=None, Seg=None):
        self.SigNames = {}
        self.EventNames = {}
        if RecordFile is None:
            if Seg is None:
                self.Seg = neo.Segment('New Seg')
            else:
                self.Seg = Seg
                self.UpdateEventDict()
                self.UpdateSignalsDict()
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
        if FileName.endswith('.h5'):
            out_f = neo.io.NixIO(filename=FileName)
        elif FileName.endswith('.mat'):
            out_f = neo.io.NeoMatlabIO(filename=FileName)
        else:
            return

        out_bl = neo.Block(name='NewBlock')
        out_bl.segments.append(self.Seg)
        out_f.write_block(out_bl)
        if FileName.endswith('.h5'):
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

    def GetSignal(self, ChName):
        return self.Seg.analogsignals[self.SigNames[ChName]]

    def AddEvent(self, Times, Name):
        eve = neo.Event(times=Times,
                        units=pq.s,
                        name=Name)

        self.Seg.events.append(eve)
        self.UpdateEventDict()

    def AddSignal(self, Sig):
        self.Seg.analogsignals.append(Sig)
        self.UpdateSignalsDict()

    def AppendSignal(self, ChName, Vect):
        sig = self.Seg.analogsignals[self.SigNames[ChName]]

        v_old = np.array(sig)
        v_new = np.vstack((v_old, Vect))

        sig_new = sig.duplicate_with_new_array(signal=v_new*sig.units)

        self.SetSignal(ChName, sig_new)


class NeoSignal():
    def __init__(self, NeoSeg=None, SigName=None, Signal=None):
        if Signal is None:
            self.Signal = NeoSeg.GetSignal(SigName)
            self.signal = NeoSeg.GetSignal(SigName)
        else:
            self.Signal = Signal
            self.signal = Signal
        self.Name = self.Signal.name

    def GetSignal(self, Time, Units=None):
        time = self.CheckTime(Time)

        sl = self.Signal.time_slice(time[0], time[1])

        if Units is None:
            return sl
        else:
            return sl.rescale(Units)

    def CheckTime(self, Time):
        if Time is None:
            return (self.Signal.t_start, self.Signal.t_stop)

        if Time[0] < self.Signal.t_start:
            Tstart = self.Signal.t_start
        else:
            Tstart = Time[0]

        if Time[1] > self.Signal.t_stop:
            Tstop = self.Signal.t_stop
        else:
            Tstop = Time[1]

        return (Tstart, Tstop)

