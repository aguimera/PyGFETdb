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

import PhyREC.PyGFETdb.DataStructures as FETdata
import PhyREC.PyGFETdb.AnalyzeData as FETana
from PhyREC.PyGFETdb.DataClass import PyFETPlotDataClass as FETplt
from PhyREC.PyGFETdb.DataClass import DataCharAC as FETcl

from multiprocessing import Pool

def RemoveGlitch(dat, th=4e-8):
    dv = np.diff(dat, axis=0)
    oinp = np.where(dv>th) 
    oinn = np.where(-dv>th) 
    for i,(s, c) in enumerate(zip(oinp[0], oinp[1])):
        s = s + 1
        dat[s, c] = np.mean((dat[s-1, c], dat[s+1, c]))
        
    for i,(s, c) in enumerate(zip(oinn[0], oinn[1])):
        dat[s, c] = np.mean((dat[s-1, c], dat[s+1, c]))

    return len(oinn[0]) + len(oinp[0])


def LoadMatFile(FileName, FBChannels=None, InChannels=None, DownFact=None, GlitchTh=None):
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

    if FBChannels is None:
        fbChs = int(nChs/2)
    else:
        fbChs = FBChannels

    if InChannels is None:
        InChannels = range(fbChs)

    if GlitchTh is not None:
        Glitchs = RemoveGlitch(Data, GlitchTh)
        Atempts = 1
        while Glitchs > 0:
            Atempts += 1
            print('Glitch removal', Atempts, ' -->> ', Glitchs)
            Glitchs = RemoveGlitch(Data, GlitchTh)
            if Atempts > 50:
                break

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
    return sig.duplicate_with_new_data(signal=v_new*sig.units)


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

    return FilesToRead, FileTimes[FileInds[0], 0] * pq.s


def LoadMatFiles(FilesIn, FBChannels=None, InChannels=None, DownFact=None, TimeWind=None,
                 Multiproces=False, GlitchTh=None):

    FbData = None
    FilesToRead, Tstart = CheckFilesTime(FilesIn, TimeWind)
    
    if Multiproces:
        pass
        # Args = []
        # for ifi, FileIn in enumerate(FilesToRead):
        #     Args.append((FileIn, InChannels, DownFact))
        # Procs = Pool(len(FilesToRead))
        # Res = Procs.starmap(LoadMatFile, Args)
        
        # for ifi, FileIn in enumerate(FilesToRead):
        #     Data = Res[ifi, 0]
        #     if FbData is None:
        #         FbData = Data
        #     else:
        #         FbData = AppendData(FbData, Data)
    else:               
        for ifi, FileIn in enumerate(FilesToRead):
            print('Tstart -->> ', Tstart)
            print(ifi, '-', len(FilesToRead), FileIn, )
            Data, _, _ = LoadMatFile(FileName=FileIn,
                                     FBChannels=FBChannels,
                                     InChannels=InChannels,
                                     DownFact=DownFact,
                                     GlitchTh=GlitchTh)
        
            if FbData is None:
                FbData = Data
                FbData.t_start = Tstart
            else:
                FbData = AppendData(FbData, Data)
    
    return FbData


def GetGtechChar(FileInput, ExcludeCh=(None, )):
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


def Calibrate(DCChar, FbData, VgsExp, CalTime=(60*pq.s, None),
              Regim='hole', PlotChar=True, CalType='interp'):

    # DCChar, _ = FETdata.LoadDataFromFile(CalFile)
    DataTrts = {}
    for tn, datadic in DCChar.items():
        DataTrts[tn] = (FETcl(datadic),)

    if PlotChar:
        pltDC = FETplt()
        pltDC.AddAxes(('Ids', 'GM', 'Rds'))
        Trts = list(DCChar.keys())
        pltDC.PlotDataSet(DataTrts, Trts, PltIsOK=True, ColorOn='Trt')
        pltDC.AddLegend()

    Vsigs = None
    for ic, sig in enumerate(FbData.transpose()):
        Csig = FbData.duplicate_with_new_data(sig)
        Csig.name = 'Ch{0:02d}'.format(ic+1)

        if Csig.name not in DCChar.keys():
            continue

        Tchar = DataTrts[Csig.name][0]
        Vsig = Spro.CalcVgeff(Csig.time_slice(*CalTime),
                              Tchar=Tchar,
                              VgsExp=VgsExp,
                              Regim=Regim,
                              CalType=CalType)

        if Vsigs is None:
            Vsigs = Vsig
        else:
            Vsigs = Vsigs.merge(Vsig)

    Vsigs.name = 'Calibrated Channels'
    return Vsigs


