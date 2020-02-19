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
