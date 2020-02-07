#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 13:35:35 2019

@author: aguimera
"""

from PhyREC.GtechInterface import LoadMatFiles, Calibrate
import glob
import os

import scipy.io as spio
from PhyREC.NeoInterface import NeoSegment, NeoSignal
import quantities as pq
import PhyREC.SignalProcess as Spro
import matplotlib.pyplot as plt
import numpy as np

from neo.core import ChannelIndex, Block, Segment, Event
from neo.io import NixIO

#%% Configuration

ChannelCount = 16

CalFile = 'B12142W46-S1-currentDCcharacterization_20191218_131755.mat_py3.h5'
FilesIn = sorted(glob.glob('../RawData/*Rec1*.mat'))
FileOut = '../CalData/Rec1_PhyREC4.h5'
FileOutCur = None
OverWrite = True

DownFactor = 2
GlitchTh = None
CalibrationTime = (60*pq.s, None)
Regim = 'hole'
VgsExp = 0

#%% Load Files and Calibrate

FbData = LoadMatFiles(FilesIn,
                      FBChannels=ChannelCount,
                      DownFact=DownFactor,
                      GlitchTh=GlitchTh,
                      TimeWind=CalibrationTime)


Vsigs = Calibrate(CalFile, FbData,
                  VgsExp=VgsExp,
                  Regim=Regim,
                  CalTime=CalibrationTime,
                  CalType='interp',
                  )


#%% Define general channel index
#%% Define epicortical channels

Blk = Block()
Seg = Segment()

EpiPitch = 400*pq.um
IndexOffSet = 0
EpiSorting = {
                'E00': (12, (0, 0, 0)),
                'E01': (11, (0, 1, 0)),
                'E02': (3, (0, 2, 0)),
                'E03': (4, (0, 3, 0)),
                'E10': (14, (1, 0, 0)),
                'E11': (9, (1, 1, 0)),
                'E12': (1, (1, 2, 0)),
                'E13': (6, (1, 3, 0)),
                'E20': (8, (2, 0, 0)),
                'E21': (15, (2, 1, 0)),
                'E22': (7, (2, 2, 0)),
                'E23': (0, (2, 3, 0)),
                'E30': (10, (3, 0, 0)),
                'E31': (13, (3, 1, 0)),
                'E32': (5, (3, 2, 0)),
                'E33': (2, (3, 3, 0)),
                }

Index = [i+IndexOffSet for k, (i, p) in EpiSorting.items()]
Names = [k for k, (i, p) in EpiSorting.items()]
PltX = [p[0] for k, (i, p) in EpiSorting.items()]
PltY = [p[1] for k, (i, p) in EpiSorting.items()]
PltZ = [p[2] for k, (i, p) in EpiSorting.items()]
Coords = [p for k, (i, p) in EpiSorting.items()] * EpiPitch

ChxEpi = ChannelIndex(index=Index,
                      name='EpiCortical',
                      channel_names=Names,
                      coordinates=Coords)

ChxEpi.annotate(**{'Pitch': EpiPitch,
                   'Shape': (4, 4),
                   'PltX': PltX,
                   'PltY': PltY,
                   'PltZ': PltZ,
                   'Names': Names,
                   })

ChxEpi.analogsignals.append(Vsigs)
Vsigs.channel_index = ChxEpi

Seg.analogsignals.append(Vsigs)

Blk.segments.append(Seg)
Blk.channel_indexes.append(ChxEpi)


#%% Events definition

# for evName, evData in Events.items():
#     times = np.array([e[0] for e in evData['Events']]) * pq.s
#     labels = np.array([e[1] for e in evData['Events']])

#     ExpEvents = Event(times=times,
#                       labels=labels,
#                       name=evName,
#                       )

#     if len(evData['DataAnnotated'])>2:           
#         array_annotations = {}
#         for iann, annoLabel in enumerate(evData['DataAnnotated'][2:]):
#             data = np.array([e[iann+2] for e in evData['Events']])
#             array_annotations[annoLabel] = data    
#         ExpEvents.annotate(**array_annotations)
    
#     Seg.events.append(ExpEvents)

#%% Save file

if os.path.isfile(FileOut):
    if OverWrite:
        os.remove(FileOut)
        print ('File removed ', FileOut)
    else:
        print ('Warning File Exsist')

Outf = NixIO(FileOut)
Outf.write(Blk)
Outf.close()








