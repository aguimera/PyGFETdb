#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:41:42 2020

@author: aguimera
"""

import numpy as np
import PhyREC.PyGFETdb.AnalyzeData as FETana
import os
import datetime
from McsPy import McsData
from neo.core import AnalogSignal
import quantities as pq
import pandas as pd
import PhyREC.SignalProcess as Spro
from scipy import interpolate
import sys


def GetMCSChar(FileInput, Vds):
    df = pd.read_csv(FileInput, encoding_errors='ignore', sep='\t')
    Time = datetime.datetime.fromtimestamp(os.path.getmtime(FileInput))

    Ids = {}
    for col in df.columns:
        if col.startswith('Vgs'):
            Vgs = df[col].values
        elif col.startswith('Vds'):
            Vds = df[col].values[0]
        elif col.startswith('Result'):
            chn = col.split(' ')[1]
            Ids[chn] = df[col].values*1e-6

    DevDCVals = {}
    for k, v in Ids.items():
        DCVals = {'Ids': v[:, None],
                  'Vds': [Vds, ],
                  'Vgs': Vgs,
                  'ChName': k,
                  'Name': k,
                  'DateTime': Time}
        DevDCVals[k] = DCVals

    FETana.CheckIsOK(DevDCVals, RdsRange=[400, 40e3])
    FETana.CalcGM(DevDCVals)

    return DevDCVals


def LoadMCSFile(FileInput, DownFactor=1, AnaStrns=(1, )):

    Dat = McsData.RawData(FileInput)
    Rec = Dat.recordings[0]
    NSamps = Rec.duration

    Sigs = []
    for AnaStrn, AnaStr in Rec.analog_streams.items():
        # if AnaStrn not in AnaStrns:
        #     continue
        # if len(AnaStr.channel_infos) == 1:
        #     continue 

        for Chn, Chinfo in AnaStr.channel_infos.items():
            print('Analog Stream ', Chinfo.label, Chinfo.sampling_frequency)
            ChName = str(Chinfo.label).split(' ')[-1]
            print(ChName)

            Fs = Chinfo.sampling_frequency
            Var, Unit = AnaStr.get_channel_in_range(Chn, 0, NSamps)

            sig = AnalogSignal(pq.Quantity(Var, pq.V),
                               t_start=0*pq.s,
                               sampling_rate=Fs.magnitude*pq.Hz,
                               name=ChName)
            sig = Spro.DownSampling(sig, Fact=DownFactor)
            Sigs.append(sig)
    return Sigs


def CalcVgeff3(Sig, Tchar, VgsExp, Regim='hole', CalType='interp'):
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
        elif CalType == 'linear':
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

    annotations = {'Calibrated': Calibrated,
                   'Working': Calibrated,
                   # 'IdsOff': IdsOff.flatten().magnitude,
                   'VgsCal': np.mean(st).magnitude,
                   'IdsBias': IdsBias.flatten()[0].magnitude,
                   'IsOK': [Tchar.IsOK, ],
                   'Iname': Sig.name,
                   'GM': GM.magnitude,
                   }

    CalSig = AnalogSignal(st,
                          units='V',
                          t_start=Sig.t_start,
                          sampling_rate=Sig.sampling_rate,
                          name=str(Sig.name),
                          file_origin=Sig.file_origin)

    CalSig.annotate(**annotations)

    return CalSig



