#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 15:08:05 2017

@author: aguimera
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpcolors
import matplotlib.cm as cmx
import sys
from itertools import cycle
import statsmodels.api as sm
import xlsxwriter as xlsw
from PyGFETdb.DBSearch import GetFromDB, FindCommonValues


def CreateCycleColors(Vals):
    ncolors = len(Vals)
    cmap = cmx.ScalarMappable(mpcolors.Normalize(vmin=0, vmax=ncolors),
                              cmx.jet)
    colors = []
    for i in range(ncolors):
        colors.append(cmap.to_rgba(i))

    return cycle(colors)


def PlotMeanStd(Data, Xvar, Yvar, Vgs=None, Vds=None, Ax=None, Ud0Norm=True,
                Color='r', PlotOverlap=False, PlotOverlapMean=False,
                PlotStd=True,
                label=None, ScaleFactor=1, **kwargs):

    fontsize = 'medium'
    labelsize = 5
    scilimits = (-2, 2)
    PointsInRange = 100

    if Ax is None:
        fig, Ax = plt.subplots()

    if PlotOverlap:
        for Trtn, Datas in Data.items():
            for Dat in Datas:
                if Dat.IsOK:
                    funcX = Dat.__getattribute__('Get' + Xvar)
                    funcY = Dat.__getattribute__('Get' + Yvar)
                    Valy = funcY(Vds=Vds, Ud0Norm=Ud0Norm,
                                 **kwargs) * ScaleFactor
                    Valx = funcX(Vds=Vds, Ud0Norm=Ud0Norm, **kwargs)
                    if Valy is not None:
                        Ax.plot(Valx, Valy, color=Color, alpha=0.2)

    # Search Vgs Vals
    VxMin = []
    VxMax = []
    for Trtn, Datas in Data.items():
        for Dat in Datas:
            if Dat.IsOK:
                funcX = Dat.__getattribute__('Get' + Xvar)
                Valx = funcX(Vds=Vds, Ud0Norm=Ud0Norm, **kwargs)
                if Valx is not None:
                    VxMax.append(np.max(Valx))
                    VxMin.append(np.min(Valx))
    VxMax = np.min(VxMax)
    VxMin = np.max(VxMin)
    if Vgs is not None:
        PointsInRange = len(Vgs)
        if np.min(Vgs) > VxMin:
            VxMin = np.min(Vgs)
        if np.max(Vgs) < VxMax:
            VxMax = np.max(Vgs)
    ValX = np.linspace(VxMin, VxMax, PointsInRange)
    ValY = np.array([])

    if 'xlsSheet' in list(kwargs.keys()):
        xlscol = 0
        kwargs['xlsSheet'].write(0, xlscol, Xvar + ' - ' + Yvar)
        for ivr, vr in enumerate(ValX):
            kwargs['xlsSheet'].write(ivr+1, xlscol, vr)

    for Trtn, Datas in Data.items():
        for Dat in Datas:
            if Dat.IsOK:
                funcY = Dat.__getattribute__('Get' + Yvar)
                Valy = funcY(Vgs=ValX, Vds=Vds, Ud0Norm=Ud0Norm,
                             **kwargs) * ScaleFactor
                if Valy is not None:
                    ValY = np.hstack((ValY, Valy)) if ValY.size else Valy

                    if 'xlsSheet' in list(kwargs.keys()):
                        xlscol = xlscol + 1
                        kwargs['xlsSheet'].write(0, xlscol, Trtn)
                        for ivr, vr in enumerate(Valy):
                            kwargs['xlsSheet'].write(ivr+1, xlscol, vr)

                    if PlotOverlapMean:
                        Ax.plot(ValX, Valy, color=Color, alpha=0.2)

    if ValY.size:
        avg = np.mean(ValY, axis=1)
        std = np.std(ValY, axis=1)
        Ax.plot(ValX, avg, color=Color, label=label)
        if PlotStd: 
            Ax.fill_between(ValX, avg+std, avg-std,
                             color=Color,
                             linewidth=0.0,
                             alpha=0.3)

    if 'xlsSheet' in list(kwargs.keys()):
        xlscol = xlscol + 1
        kwargs['xlsSheet'].write(0, xlscol, 'Avg')
        for ivr, vr in enumerate(avg):
            kwargs['xlsSheet'].write(ivr+1, xlscol, vr)

    if 'xlsSheet' in list(kwargs.keys()):
        xlscol = xlscol + 1
        kwargs['xlsSheet'].write(0, xlscol, 'Std')
        for ivr, vr in enumerate(std):
            kwargs['xlsSheet'].write(ivr+1, xlscol, vr)

    if 'xscale' in list(kwargs.keys()):
        Ax.set_xscale(kwargs['xscale'])
    else:
        Ax.ticklabel_format(axis='x', style='sci', scilimits=scilimits)

    if 'yscale' in list(kwargs.keys()):
        Ax.set_yscale(kwargs['yscale'])
    else:
        Ax.ticklabel_format(axis='y', style='sci', scilimits=scilimits)

    Ax.set_ylabel(Yvar, fontsize=fontsize)
    Ax.set_xlabel(Xvar, fontsize=fontsize)
    Ax.tick_params(axis='both', which='Both', labelsize=labelsize)


def PlotXYVars(Data, Xvar, Yvar, Vgs, Vds, Ud0Norm=True, label=None,
               Ax=None, Color=None, **kwargs):

    fontsize = 'medium'
    labelsize = 5
    scilimits = (-2, 2)

    if Ax is None:
        fig, Ax = plt.subplots()

    for Trtn, Datas in Data.items():
        for Dat in Datas:
            if Dat.IsOK:
                funcX = Dat.__getattribute__('Get' + Xvar)
                funcY = Dat.__getattribute__('Get' + Yvar)

                try:
                    Valx = funcX(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm, **kwargs)
                    Valy = funcY(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm, **kwargs)
                    Ax.plot(Valx, Valy, '*', color=Color, label=label)
                except:  # catch *all* exceptions
                    print(Dat.Name, sys.exc_info()[0])

    if 'xscale' in list(kwargs.keys()):
        Ax.set_xscale(kwargs['xscale'])
    elif Xvar != 'DateTime':
        try:
            Ax.ticklabel_format(axis='x', style='sci', scilimits=scilimits)
        except:
            print('Formating error')

    if 'yscale' in list(kwargs.keys()):
        Ax.set_yscale(kwargs['yscale'])
    else:
        Ax.ticklabel_format(axis='y', style='sci', scilimits=scilimits)

    if 'ylim' in list(kwargs.keys()):
        Ax.set_ylim(kwargs['ylim'])

    Ax.set_ylabel(Yvar, fontsize=fontsize)
    Ax.set_xlabel(Xvar, fontsize=fontsize)
    Ax.tick_params(axis='both', which='Both', labelsize=labelsize)
    

def GetParam(Data, Param, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
    Vals = np.array([])
    for Trtn, Datas in Data.items():
        for Dat in Datas:
            func = Dat.__getattribute__('Get' + Param)

            try:
                Val = func(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm, **kwargs)
            except:  # catch *all* exceptions
                print(Dat.Name, sys.exc_info()[0])
                Val = None
    
            if Val is not None:
                Vals = np.hstack((Vals, Val)) if Vals.size else Val
    return Vals


def SearchAndGetParam(Groups, Plot=True, Boxplot=False, ParamUnits=None, **kwargs):
    if Plot:
        fig, Ax = plt.subplots()
        xLab = []
        xPos = []

    if 'XlsFile' in list(kwargs.keys()):
        xlswbook = xlsw.Workbook(kwargs['XlsFile'])
        xlssheet = xlswbook.add_worksheet('W1')

    Vals = {}
    for iGr, (Grn, Grc) in enumerate(sorted(Groups.items())):
        print('Getting data for ', Grn)
        Data, Trts = GetFromDB(**Grc)

        if len(Data) > 0:
            vals = GetParam(Data, **kwargs)
            if vals is None:
                continue
            if vals.size == 0:
                continue

            vals = vals[~np.isnan(vals)]
            Vals[Grn] = vals

            if 'XlsFile' in list(kwargs.keys()):
                xlssheet.write(0, iGr, Grn)
                for ivr, vr in enumerate(vals):
                    xlssheet.write(ivr+1, iGr, vr)

            if Plot:
                if Boxplot:
                    Ax.boxplot(vals.transpose(), positions=(iGr+1,))
                    xPos.append(iGr+1)
                else:
                    Ax.plot(np.ones(len(vals))*iGr, vals, '*')
                    xPos.append(iGr)
                xLab.append(Grn)
        else:
            print('Empty data for ', Grn)

    if Plot:
        plt.xticks(xPos, xLab, rotation=45)
        if ParamUnits is not None:
            Ax.set_ylabel(kwargs['Param']+ ParamUnits)
        else:
            Ax.set_ylabel(kwargs['Param'])
            
        Ax.grid()
        Ax.ticklabel_format(axis='y', style='sci', scilimits=(2, 2))
        if len(xPos) > 1:
            Ax.set_xlim(min(xPos)-0.5, max(xPos)+0.5)
        if 'Vgs' in kwargs and 'Vds' in kwargs:
            title = 'Vgs {} Vds {}'.format(kwargs['Vgs'], kwargs['Vds'])
            plt.title(title)
        plt.tight_layout()
        if 'xscale' in list(kwargs.keys()):
            Ax.set_xscale(kwargs['xscale'])
        if 'yscale' in list(kwargs.keys()):
            Ax.set_yscale(kwargs['yscale'])

    if 'XlsFile' in list(kwargs.keys()):
        xlswbook.close()

    return Vals


def SearchAndPlot(Groups, Func=PlotMeanStd, **kwargs):
    col = CreateCycleColors(Groups)

    if 'Ax' not in list(kwargs.keys()):
        fig, Ax = plt.subplots()
        kwargs['Ax'] = Ax

    if 'XlsFile' in list(kwargs.keys()):
        xlswbook = xlsw.Workbook(kwargs['XlsFile'])

    for Grn, Grc in sorted(Groups.items()):
        print('Getting data for ', Grn)
        Data, Trts = GetFromDB(**Grc)
        if len(Data) > 0:
            try:
                if 'XlsFile' in list(kwargs.keys()):
                    xlssheet = xlswbook.add_worksheet(str(Grn))
                    kwargs['xlsSheet'] = xlssheet
            except:
                print('Error in excel generation')

            try:
                Func(Data,
                     Color=next(col),
                     label=Grn,
                     **kwargs)
            except:
                print(Grn, 'ERROR --> ', sys.exc_info()[0])
        else:
            print('Empty data for ', Grn)

    Ax = kwargs['Ax']

    handles, labels = Ax.get_legend_handles_labels()
    hh = []
    ll = []
    for h, l in zip(handles, labels):
        if l not in ll:
            hh.append(h)
            ll.append(l)
    Ax.legend(hh, ll)

    if 'XlsFile' in list(kwargs.keys()):
        xlswbook.close()

    return plt.gcf(), Ax


def PlotGroupBy(GroupBase, GroupBy, **kwargs):

    GroupList = FindCommonValues(Table=GroupBase['Table'],
                                 Conditions=GroupBase['Conditions'],
                                 Parameter=GroupBy)

    Groups = {}
    for Item in sorted(GroupList):
        Cgr = GroupBase.copy()
        Cond = GroupBase['Conditions'].copy()
        Cond.update({'{}='.format(GroupBy): (Item,)})
        Cgr['Conditions'] = Cond
        Groups[Item] = Cgr

    return SearchAndPlot(Groups=Groups, **kwargs)

def PlotGroupBySearchAndGetParam(GroupBase, GroupBy, **kwargs):

    GroupList = FindCommonValues(Table=GroupBase['Table'],
                                 Conditions=GroupBase['Conditions'],
                                 Parameter=GroupBy)

    Groups = {}
    for Item in sorted(GroupList):
        Cgr = GroupBase.copy()
        Cond = GroupBase['Conditions'].copy()
        Cond.update({'{}='.format(GroupBy): (Item,)})
        Cgr['Conditions'] = Cond
        Groups[Item] = Cgr

    return SearchAndGetParam(Groups=Groups, **kwargs)


def CalcTLM(Groups, Vds=None, Ax=None, Color=None,
            DebugPlot=False, Label=None):
    if Ax is None:
        fig,  AxRs = plt.subplots()
        AxRc = AxRs.twinx()
        fig1,  AxLT = plt.subplots()
    else:
        AxRc = Ax[0]
        AxRs = Ax[1]
        AxLT = Ax[2]

    PointsInRange = 100
    DatV = []
    for Grn, Grc in sorted(Groups.items()):
        print('Getting data for ', Grn)
        Data, Trts = GetFromDB(**Grc)
        DatV.append(Data)

    VxMin = []
    VxMax = []
    for Data in DatV:
        if len(Data) > 0:
            for Trtn, Datas in Data.items():
                for Dat in Datas:
                    funcX = Dat.__getattribute__('GetVgs')
                    Valx = funcX(Vds=Vds, Ud0Norm=True)
                    if Valx is not None:
                        VxMax.append(np.max(Valx))
                        VxMin.append(np.min(Valx))

    VxMax = np.min(VxMax)
    VxMin = np.max(VxMin)
    VGS = np.linspace(VxMin, VxMax, PointsInRange)

    Rsheet = np.ones(VGS.size)*np.NaN
    RsheetMax = np.ones(VGS.size)*np.NaN
    RsheetMin = np.ones(VGS.size)*np.NaN
    Rc = np.ones(VGS.size)*np.NaN
    RcMax = np.ones(VGS.size)*np.NaN
    RcMin = np.ones(VGS.size)*np.NaN
    LT = np.ones(VGS.size)*np.NaN

    if DebugPlot:
        plt.figure()
    Width = None
    for ivg, vgs in enumerate(VGS):
        X = np.array([])
        Y = np.array([])
        for Data in DatV:
            if len(Data) > 0:
                for Trtn, Datas in Data.items():
                    for Dat in Datas:
                        rds = Dat.GetRds(Vgs=vgs, Vds=Vds, Ud0Norm=True)
                        Y = np.vstack((Y, rds)) if Y.size else rds
                        L = np.array((Dat.TrtTypes['Length']))
                        X = np.vstack((X, L)) if X.size else L
                        if Width is None:
                            Width = Dat.TrtTypes['Width']
                        else:
                            if not Width == Dat.TrtTypes['Width']:
                                print(Trtn, 'WARNING Bad width')
        if DebugPlot:
            plt.plot(X, Y, '*')

        X = sm.add_constant(X)
        res = sm.OLS(Y, X).fit()
        Rsheet[ivg] = res.params[1] * Width
        RsheetMax[ivg] = (res.bse[1]+res.params[1]) * Width
        RsheetMin[ivg] = (-res.bse[1]+res.params[1]) * Width
        Rc[ivg] = res.params[0]
        RcMax[ivg] = res.bse[0]+res.params[0]
        RcMin[ivg] = -res.bse[0]+res.params[0]
        LT[ivg] = (res.params[0]/res.params[1])/2

    AxRc.plot(VGS, Rc, color=Color, label=Label)
    AxRc.fill_between(VGS, RcMax, RcMin,
                      color=Color,
                      linewidth=0.0,
                      alpha=0.3)
    AxRs.plot(VGS, Rsheet, '--', color=Color)
    AxRs.fill_between(VGS, RsheetMax, RsheetMin,
                      color=Color,
                      linewidth=0.0,
                      alpha=0.3)

    AxLT.plot(VGS, LT, color=Color)

    AxRc.set_ylabel('Rc')
    AxRs.set_ylabel('Rsheet')
    AxRc.legend()

    ContactVals = {}
    ContactVals['Rsheet'] = Rsheet
    ContactVals['RsheetMax'] = RsheetMax
    ContactVals['RsheetMin'] = RsheetMin
    ContactVals['Rc'] = Rc
    ContactVals['RcMax'] = RcMax
    ContactVals['RcMin'] = RcMin
    ContactVals['VGS'] = VGS
    ContactVals['LT'] = LT

    return ContactVals


def CalcTLM2(Groups, Vds=None, Ax=None, Color=None,
            DebugPlot=False, Label=None, Lerror=0.4e-6, TrackResistance=None):
    if Ax is None:
        fig,  AxRs = plt.subplots()
        AxRc = AxRs.twinx()
        fig1,  AxLT = plt.subplots()
    else:
        AxRc = Ax[0]
        AxRs = Ax[1]
        AxLT = Ax[2]

    PointsInRange = 100
    DatV = []
    for Grn, Grc in sorted(Groups.items()):
        print('Getting data for ', Grn)
        Data, Trts = GetFromDB(**Grc)
        DatV.append(Data)

    # Find common Vgs Range
    VxMin = []
    VxMax = []
    for Data in DatV:
        if len(Data) > 0:
            for Trtn, Datas in Data.items():
                for Dat in Datas:
                    funcX = Dat.__getattribute__('GetVgs')
                    Valx = funcX(Vds=Vds, Ud0Norm=True)
                    if Valx is not None:
                        VxMax.append(np.max(Valx))
                        VxMin.append(np.min(Valx))

    VxMax = np.min(VxMax)
    VxMin = np.max(VxMin)
    VGS = np.linspace(VxMin, VxMax, PointsInRange)  # Generate VGS common range

    Rsheet = np.ones(VGS.size)*np.NaN
    RsheetMax = np.ones(VGS.size)*np.NaN
    RsheetMin = np.ones(VGS.size)*np.NaN
    Rc = np.ones(VGS.size)*np.NaN
    RcMax = np.ones(VGS.size)*np.NaN
    RcMin = np.ones(VGS.size)*np.NaN
    LT = np.ones(VGS.size)*np.NaN

    if DebugPlot:
        plt.figure()
    Width = None
    for ivg, vgs in enumerate(VGS):
        X = np.array([])
        Y = np.array([])
        for Data in DatV:
            if len(Data) > 0:
                for Trtn, Datas in Data.items():
                    for Dat in Datas:
                        rds = Dat.GetRds(Vgs=vgs, Vds=Vds, Ud0Norm=True)                        
                        if TrackResistance is not None:
                            rds = rds - TrackResistance[Dat.Name.split('-')[-1][1:]]
                        Y = np.vstack((Y, rds)) if Y.size else rds
                        L = np.array((Dat.TrtTypes['Length']))
                        L = L - Lerror
                        X = np.vstack((X, L)) if X.size else L
                        if Width is None:
                            Width = Dat.TrtTypes['Width']
                        else:
                            if not Width == Dat.TrtTypes['Width']:
                                print(Trtn, 'WARNING Bad width')
        if DebugPlot:
            plt.plot(X, Y, '*')

        X = sm.add_constant(X)
        res = sm.OLS(Y, X).fit()
        Rsheet[ivg] = res.params[1] * Width
        RsheetMax[ivg] = (res.bse[1]+res.params[1]) * Width
        RsheetMin[ivg] = (-res.bse[1]+res.params[1]) * Width
        Rc[ivg] = (res.params[0]/2) * (Width*1e6)
        RcError = (res.bse[0]/2) * (Width*1e6)
        RcMax[ivg] = RcError+Rc[ivg]
        RcMin[ivg] = -RcError+Rc[ivg]
        LT[ivg] = (res.params[0]/res.params[1])/2

    AxRc.plot(VGS, Rc, color=Color, label=Label)
    AxRc.fill_between(VGS, RcMax, RcMin,
                      color=Color,
                      linewidth=0.0,
                      alpha=0.3)
    AxRs.plot(VGS, Rsheet, '--', color=Color)
    AxRs.fill_between(VGS, RsheetMax, RsheetMin,
                      color=Color,
                      linewidth=0.0,
                      alpha=0.3)

    AxLT.plot(VGS, LT, color=Color)

    AxRc.set_ylabel('Rc')
    AxRs.set_ylabel('Rsheet')
    AxRc.legend()

    ContactVals = {}
    ContactVals['Rsheet'] = Rsheet
    ContactVals['RsheetMax'] = RsheetMax
    ContactVals['RsheetMin'] = RsheetMin
    ContactVals['Rc'] = Rc
    ContactVals['RcMax'] = RcMax
    ContactVals['RcMin'] = RcMin
    ContactVals['VGS'] = VGS
    ContactVals['LT'] = LT

    return ContactVals