#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:01:03 2020

@author: aguimera
"""

import copy
import PyGFETdb.DBSearch as DBs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib import colors
import matplotlib.gridspec as gridspec
import seaborn as sns
import math

BioSensMap = {'Shape': (6, 8),
              'Ch15OE': (0, 0),
              'Ch12OE': (1, 0),
              'Ch05OE': (2, 0),
              'Ch13OE': (3, 0),
              'Ch04OE': (4, 0),
              'Ch07OE': (5, 0),

              'Ch07NS': (0, 1),
              'Ch16OE': (1, 1),
              'Ch02OE': (2, 1),
              'Ch10OE': (3, 1),
              'Ch08OE': (4, 1),
              'Ch31NS': (5, 1),

              'Ch04NS': (0, 2),
              'Ch08NS': (1, 2),
              'Ch06OE': (2, 2),
              'Ch14OE': (3, 2),
              'Ch32NS': (4, 2),
              'Ch28NS': (5, 2),

              'Ch13NS': (0, 3),
              'Ch10NS': (1, 3),
              'Ch14NS': (2, 3),
              'Ch22NS': (3, 3),
              'Ch18NS': (4, 3),
              'Ch21NS': (5, 3),

              'Ch05NS': (0, 4),
              'Ch02NS': (1, 4),
              'Ch06NS': (2, 4),
              'Ch30NS': (3, 4),
              'Ch26NS': (4, 4),
              'Ch29NS': (5, 4),

              'Ch12NS': (0, 5),
              'Ch16NS': (1, 5),
              'Ch30OE': (2, 5),
              'Ch22OE': (3, 5),
              'Ch24NS': (4, 5),
              'Ch20NS': (5, 5),

              'Ch15NS': (0, 6),
              'Ch24OE': (1, 6),
              'Ch26OE': (2, 6),
              'Ch18OE': (3, 6),
              'Ch32OE': (4, 6),
              'Ch23NS': (5, 6),

              'Ch23OE': (0, 7),
              'Ch20OE': (1, 7),
              'Ch29OE': (2, 7),
              'Ch21OE': (3, 7),
              'Ch28OE': (4, 7),
              'Ch31OE': (5, 7),

               }


def GetBiosensSeriesDataSet(DevicesList, FuncStepsSort,
                            CharTable='DCcharacts',
                            AnalyteStep = ('Tromb', 'BSA'),
                            ValuesField = 'CharTable.AnalyteCon',
                            ):
    
    pdSeries = []
    for dev in DevicesList:
        devs = (dev+'NS', dev+'OE')
        Conditions = {'Devices.name = ': devs,
                      'CharTable.IsOK >': (0, ),
                      }
        GrBase = {'Conditions': Conditions,
                  'Table': CharTable,
                  'Last': True,
                  'GetGate': True,
                  }
        GTrts = DBs.GenGroups(GrBase,
                              GroupBy='Trts.name',
                              LongName=False)
    
        for Trtn, Gtrt in GTrts.items():
            Grs = DBs.GenGroups(Gtrt,
                                GroupBy='CharTable.FuncStep',
                                LongName=False)        
            for fn, Gr in Grs.items():
                s = Trtn.split('-')
                s[2] = s[2] + s[1][-2:]
                TnameBio = '-'.join(s)
    
                pdser = {}
                pdser['ChBioName'] = TnameBio.split('-')[-1]
                pdser['TrtName'] = Trtn
                pdser['DevName'] = dev
                # pdser['Batch'] = DevicesSortBatch[dev]
                # pdser['Linker'] = DevicesSortLinker[dev]
                pdser['FuncStep'] = FuncStepsSort[fn]
                pdser['MAPx'] = BioSensMap[TnameBio.split('-')[-1]][0]
                pdser['MAPy'] = BioSensMap[TnameBio.split('-')[-1]][1]
                
                if FuncStepsSort[fn] in AnalyteStep:
                    Gras = DBs.GenGroups(Gr,
                                          GroupBy=ValuesField,
                                          LongName=False)
                    for ian, (fna, Gra) in enumerate(Gras.items()):
                        pdsera = copy.deepcopy(pdser)
                        d, t = DBs.GetFromDB(**Gra)
                        dat = d[t[0]][0]
                        pdsera['Analyte'] = np.log10(dat.GetAnalyteCon().flatten())
                        pdsera['CharCl'] = dat
                        pdSeries.append(pdsera)
                else:
                    d, t = DBs.GetFromDB(**Gr)
                    dat = d[t[0]][0]
                    pdser['CharCl'] = dat
                    pdSeries.append(pdser)   
    return pdSeries


def CalcNormalization(df, RefStep, ScalarValues, Experiment=None,
                      StepField='FuncStep', SortTrtField='TrtName', **kwargs):
    pdSeries = []
    gTrts = df.groupby(SortTrtField, observed=True)
    for nTrt in gTrts.groups:
        Trt = gTrts.get_group(nTrt)
        gFuncs = Trt.groupby(StepField, observed=True)
        if RefStep in gFuncs.groups:
            refVal = gFuncs.get_group(RefStep).iloc[0,:]
        else:
            print(nTrt, 'Not Ref val')
            continue
        for nFunc in gFuncs.groups:
            Func = gFuncs.get_group(nFunc)
            ic, ir = Func.shape
            for icc in range(ic):
                gds = Func.iloc[icc,:]
                for par in ScalarValues: 
                    gds[par + 'N' + str(RefStep)] = refVal[par] - gds[par]
                if Experiment is not None:
                    gds['Experiment'] = Experiment
                pdSeries.append(gds)
    
    return pd.concat(pdSeries, axis=1).transpose()
    

def GenPDFLinePlots(df, GroupBy1, GroupBy2, vPars, PDF, Vgs, VgsNorm, cmap='jet'):
    plt.ioff()
    
    dg1 = df.groupby(GroupBy1, observed=True)

    for g1n in dg1.groups:
        fig, axs = plt.subplots(1, len(vPars))
        g1 = dg1.get_group(g1n)
        dg2 = g1.groupby(GroupBy2, observed=True)
        
        # create color map
        cm = ScalarMappable(norm=colors.Normalize(vmin=0, vmax=len(dg2)),
                            cmap=cmap)    
        G2ColorDict = {}
        for ic, k in enumerate(dg2.groups):
            G2ColorDict[k] = cm.to_rgba(ic)

        for g2n in dg2.groups:
            g2 = dg2.get_group(g2n)
            for ip, par in enumerate(vPars):
                y = np.array([])
                for index, row in g2.iterrows():
                    v = row[par].flatten()
                    y = np.vstack((y, v)) if y.size else v
                axs[ip].plot(Vgs, y.transpose(), color=G2ColorDict[g2n], label=g2n)
                y = np.array([])
                for index, row in g2.iterrows():
                    v = row[par+'Norm'].flatten()
                    y = np.vstack((y, v)) if y.size else v
                axs[ip].plot(VgsNorm, y.transpose(), '-.',
                             color=G2ColorDict[g2n],
                             alpha=0.5)
        axs[ip].legend()
        plt.title(g1n)
        PDF.savefig(fig)
        plt.close(fig)
    plt.ion()


def GenHeatMaps(df, GroupBy, value, nRows=2, figsize=(12, 7), **kwargs):

    dg = df.groupby(GroupBy)    
    nCols = math.ceil((len(dg.groups)-1)/nRows)

    Widths = np.ones(shape=(nCols+1, 1)).flatten()
    Widths[-1] = 0.1

    Fig = plt.figure(figsize=figsize, constrained_layout=True)
    gspec = gridspec.GridSpec(nrows=nRows,
                              ncols=nCols+1,
                              figure=Fig,
                              width_ratios=Widths)
    Axs = []
    Ax = None
    for ir in range(nRows):
        for ic in range(nCols):
            Ax = Fig.add_subplot(gspec[ir, ic],
                                 # sharey=Ax,
                                 # sharex=Ax,
                                 )
            Axs.append(Ax) 
    Cax = Fig.add_subplot(gspec[:,-1])
    
    ic = -1
    for ffg in dg.groups:
        if ffg in ('Tromb', 'BSA'):
            continue
        ic += 1
        d = dg.get_group(ffg)
        m = d.pivot('MAPx', 'MAPy', value).astype(np.float)
        Axs[ic].set_title(ffg)
        sns.heatmap(m,
                    ax=Axs[ic],
                    cbar_ax=Cax,
                    **kwargs)    
    
    Cax.set_ylabel(value)

# def GetStepData(Data, Parkwargs, RefSteps=None):
#     Data2Norm = copy.deepcopy(Data)
    
#     for devn, devD in Data2Norm.items():        
#         for trtn, trtD in devD.items():            
#             for fn, fnD in trtD.items():
#                 if type(fnD)==list:
#                     trtD[fn] = fnD
#                     continue
#                 dat = fnD.Get(**Parkwargs)
#                 trtD[fn] = dat.flatten()[0]
    
#     if RefSteps is not None:
#         for devn, devD in Data2Norm.items():
#             ToRemove = []
#             for trtn, trtD in devD.items():
#                 Ref = None
#                 for fn, fnD in trtD.items():
#                     if fn in RefSteps:
#                         Ref = fnD                    
#                 if Ref is None:
#                     print('Ref not found ->>', trtn, list(trtD.keys()))
#                     ToRemove.append(trtn)
#                     continue
#                 for fn, fnD in trtD.items():
#                     if type(fnD)==list:
#                         trtD[fn] = fnD
#                         continue
#                     trtD[fn] = fnD - Ref
            
#             if len(ToRemove):
#                 print('\n ************ REMOVING !!!!!!!!')
#                 print( ToRemove, '\n')
#                 for rm in ToRemove:
#                     devD.pop(rm)

#     return Data2Norm


# def GetDevicesData(DevicesList, AnalyteStep, ValuesField,
#                    CharTable='DCcharacts' ):

#     Data = {}
#     for dev in DevicesList:
#         devs = (dev+'NS', dev+'OE')
#         Conditions = {'Devices.name = ': devs,
#                       'CharTable.IsOK >': (0, ),
#                       }
#         GrBase = {'Conditions': Conditions,
#                   'Table': CharTable,
#                   'Last': True,
#                   'GetGate': True,
#                   }
#         GTrts = DBs.GenGroups(GrBase,
#                               GroupBy='Trts.name',
#                               LongName=False)
    
#         datdev = {}
#         for Trtn, Gtrt in GTrts.items():
#             Grs = DBs.GenGroups(Gtrt,
#                                 GroupBy='CharTable.FuncStep',
#                                 LongName=False)
#             dattrt = {}
#             for fn, Gr in Grs.items():
#                 if fn == AnalyteStep:
#                     Gras = DBs.GenGroups(Gr,
#                                          GroupBy=ValuesField,
#                                          LongName=False)
#                     anVals = []
#                     for ian, (fna, Gra) in enumerate(Gras.items()):
#                         anVals.append(fna)
#                         d, t = DBs.GetFromDB(**Gra)
#                         dattrt[fn+str(ian)] = d[t[0]][0]
#                     dattrt[fn] = anVals
#                 else:
#                     d, t = DBs.GetFromDB(**Gr)
#                     dattrt[fn] = d[t[0]][0]
            
#             s = Trtn.split('-')
#             s[2] = s[2] + s[1][-2:]
#             TnameBio = '-'.join(s)
#             datdev[TnameBio] = dattrt
#         Data[dev] = datdev

#     return Data

