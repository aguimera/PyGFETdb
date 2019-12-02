#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 12:09:07 2019

@author: aguimera
"""

import scipy.io as spio
from PhyREC.NeoInterface import NeoSignal
import PhyREC.SignalProcess as Spro
import quantities as pq
import numpy as np
import datetime
import deepdish as dd
import os

import PyGFETdb.PlotDataClass as FETplt
import PyGFETdb.DataStructures as FETdata
import PyGFETdb.DataClass as FETcl
import PyGFETdb.AnalyzeData as FETana


def LoadMatFile(FileName, InChannels=None, DownFact=None):
    """
    Load single file
    """
    out = spio.loadmat(FileName)
    Fs = out['SR\x00'][0, 0]*pq.Hz

    Time = out['y'][0, :].transpose()[:, 0]

    if Fs != int(1/(Time[1]-Time[0])):
        print('Warning FS', Fs, 1/(Time[1]-Time[0]))
    print('Detected sampling rate ', Fs)

    Data = out['y'][1:, 0, :].transpose()
    (Samps, nChs) = Data.shape

    fbChs = int(nChs/2)

    if InChannels is None:
        InChannels = range(fbChs)


    DCData = Spro.Filter(NeoSignal(Data[:, :fbChs][:, InChannels],
                                   units='A',
                                   t_start=Time[0]*pq.s,
                                   sampling_rate=Fs,
                                   name='DCChannels',
                                   file_origin=FileName.split('/')[-1],),
                         Type='lowpass',
                         Order=2,
                         Freqs=(1,)
                         )

    ACData = Spro.RemoveDC(NeoSignal(Data[:, fbChs:][:, InChannels],
                                     units='A',
                                     t_start=Time[0]*pq.s,
                                     sampling_rate=Fs,
                                     name='ACChannels',
                                     file_origin=FileName.split('/')[-1],),
                           Type='linear',
                           )

    if DownFact is not None:
        ACData = Spro.DownSampling(ACData,
                                   Fact=DownFact)
        DCData = Spro.DownSampling(DCData,
                                   Fact=DownFact)

    FBData = ACData + DCData
    FBData.name = 'FullBandChannels'

    return FBData, ACData, DCData


def AppendData(sig, data):
    v_old = np.array(sig)
    v_new = np.vstack((v_old, np.array(data)))
    return sig.duplicate_with_new_array(signal=v_new*sig.units)


def CheckFilesTime(FilesIn, TimeWind):
    if TimeWind is None:
        return FilesIn

    TimeWind = list(TimeWind)
    if TimeWind[0] is None:
        TimeWind[0] = 0 * pq.s

    if TimeWind[1] is None:
        TimeWind[1] = np.inf * pq.s

    Tstart = 0 * pq.s
    FileTimes = []
    for FileName in FilesIn:
        out = spio.loadmat(FileName)
        Data = out['y']
        Fs = out['SR\x00'][0, 0]*pq.Hz

        Tstop = Tstart + Data[0, 0, -1] * pq.s
        FileTimes.append((Tstart.copy(), Tstop.copy()))
        Tstart = Tstop + 1/Fs

    FileTimes = np.array(FileTimes)
    FileInds = np.where((FileTimes[:, 1] > TimeWind[0]) &
                        (FileTimes[:, 0] < TimeWind[1]))[0]

    FileInds = list(FileInds.flatten())
    FilesToRead = [FilesIn[i] for i in FileInds]

    return FilesToRead


def LoadMatFiles(FilesIn, InChannels=None, DownFact=None, TimeWind=None):
    FbData = None
    FilesToRead = CheckFilesTime(FilesIn, TimeWind)
    for ifi, FileIn in enumerate(FilesToRead):
        print(ifi, '-', len(FilesToRead), FileIn, )
        Data, _, _ = LoadMatFile(FileName=FileIn,
                                 InChannels=InChannels,
                                 DownFact=DownFact)

        if FbData is None:
            FbData = Data
        else:
            FbData = AppendData(FbData, Data)
    return FbData


def GetGtechChar(FileInput, ExcludeCh=(None,)):
    MatCon = spio.loadmat(FileInput)
    Data = np.array(MatCon['DCcharacterization'])

    Time = datetime.datetime.fromtimestamp(os.path.getmtime(FileInput))

    DevDCVals = {}
    for i, dat in enumerate(Data[1:, :]):
        ChName = 'Ch{0:02d}'.format(i+1)
        if ChName in ExcludeCh:
            continue
        Vgs = Data[0, :]
        Ids = np.array([dat, ]).transpose()
        Vds = np.array((0.05, ))

        DCVals = {'Ids': Ids,
                  'Vds': Vds,
                  'Vgs': Vgs,
                  'ChName': ChName,
                  'Name': ChName,
                  'DateTime': Time}
        DevDCVals[ChName] = DCVals

    FETana.CheckIsOK(DevDCVals, RdsRange=[400, 40e3])
    FETana.CalcGM(DevDCVals)
    return DevDCVals


def Calibrate(CalFile, FbData, VgsExp, CalTime=(60*pq.s, None),
              Regim='hole', PlotChar=True):

    DCChar, _ = FETdata.LoadDataFromFile(CalFile)

    if PlotChar:
        pltDC = FETplt.PyFETPlot()
        pltDC.AddAxes(('Ids', 'Gm', 'Rds'))
        pltDC.PlotDataCh(DCChar)
        pltDC.AddLegend()

    Vsigs = None
    for ic, sig in enumerate(FbData.transpose()):
        Csig = FbData.duplicate_with_new_array(sig)
        Csig.name = 'Ch{0:02d}'.format(ic+1)

        if Csig.name not in DCChar.keys():
            continue

        Tchar = FETcl.DataCharDC(DCChar[Csig.name])
        Vsig = Spro.CalcVgeff(Csig.GetSignal(CalTime),
                              Tchar=Tchar,
                              VgsExp=VgsExp,
                              Regim=Regim)

        if Vsigs is None:
            Vsigs = Vsig
        else:
            Vsigs = Vsigs.merge(Vsig)

    Vsigs.name = 'Calibrated Channels'
    return Vsigs
