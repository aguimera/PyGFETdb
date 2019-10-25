#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 12:09:07 2019

@author: aguimera
"""

import scipy.io as spio
from PhyREC.NeoInterface import NeoSegment, NeoSignal
import quantities as pq
 

def LoadGetchData(FileName, ExcludeList=()):
    out = spio.loadmat(FileName)
    Data = out['y']
    Fs = out['SR\x00'][0, 0]*pq.Hz

    Time = Data[0, :].transpose()[:, 0]

    if Fs != int(1/(Time[1]-Time[0])):
        print('Warning FS', Fs, 1/(Time[1]-Time[0]))

    print('Detected sampling rate ', Fs)

    ##################################
    ### Load from file in Tmp seg
    TmpSeg = NeoSegment()
    ChCount = int(Data[1:, :, :].shape[0]/2)
    print('Channels ', ChCount)
    for i, dat in enumerate(Data[1:, :, :]):
        if i < ChCount:
            ChN = 'Ch{0:02d}'.format(i+1)
            ChName = 'DC-' + ChN
        elif i >= ChCount:
            ChN = 'Ch{0:02d}'.format(i-ChCount+1)
            ChName = 'AC-' + ChN

        if ChN in ExcludeList:
            continue

        print('Reading channel ', i , ' As ', ChName)

        val = dat.transpose()[:, 0]
        sig = NeoSignal(val,
                        units='A',
                        t_start=Time[0]*pq.s,
                        sampling_rate=Fs,
                        name=ChName,
                        file_origin=FileName.split('/')[-1])
        TmpSeg.AddSignal(sig)
    return TmpSeg
