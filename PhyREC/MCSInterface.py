#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:41:42 2020

@author: aguimera
"""
import csv
import numpy as np
import PyGFETdb.AnalyzeData as FETana
import os
import datetime
from McsPy import McsData
import neo
import quantities as pq

def GetMCSChar(FileInput, Vds, ExcludeCh=(None, )):
    Fin = open(FileInput, encoding='latin')

    Time = datetime.datetime.fromtimestamp(os.path.getmtime(FileInput))

    reader = csv.reader(Fin, delimiter='\t')

    Vgs = []
    Ids = []
    for ir, r in enumerate(reader):
        if ir == 0:
            continue
        Vgs.append(np.float(r[0].replace(',', '.')))
        ids = []
        for i in r[1:]:
            ids.append(np.float(i.replace(',', '.')))
        Ids.append(ids)

    Ids = np.array(Ids)*1e-6
    Vgs = np.array(Vgs)
    Vds = np.array([Vds, ])

    DevDCVals = {}
    for i, ids in enumerate(Ids.transpose()):
        ChName = 'Ch{0:02d}'.format(i+1)
        if ChName in ExcludeCh:
            continue
        Vgs = Vgs
        DCVals = {'Ids': ids[:, None],
                  'Vds': Vds,
                  'Vgs': Vgs,
                  'ChName': ChName,
                  'Name': ChName,
                  'DateTime': Time}
        DevDCVals[ChName] = DCVals

    FETana.CheckIsOK(DevDCVals, RdsRange=[400, 40e3])
    FETana.CalcGM(DevDCVals)

    return DevDCVals

def LoadMCSFile(FileInput):
    Dat = McsData.RawData(FileInput)
    Rec = Dat.recordings[0]
    NSamps = Rec.duration

    Seg = neo.Segment()

    Sigs = None
    for AnaStrn, AnaStr in Rec.analog_streams.items():
        if len(AnaStr.channel_infos) == 1:
            continue 
        
        for Chn, Chinfo in AnaStr.channel_infos.items():
            print('Analog Stream ', Chinfo.label, Chinfo.sampling_frequency)
            ChName = str(Chinfo.label)
            print(ChName)
 
            Fs = Chinfo.sampling_frequency
            Var, Unit = AnaStr.get_channel_in_range(Chn, 0, NSamps)
            sig = neo.AnalogSignal(pq.Quantity(Var, pq.V),
                                   t_start=0*pq.s,
                                   sampling_rate=Fs.magnitude*pq.Hz,
                                   name=ChName)
 
            if Sigs is None:
                Sigs = sig
            else:
                Sigs = Sigs.merge(sig)

    return Sigs



