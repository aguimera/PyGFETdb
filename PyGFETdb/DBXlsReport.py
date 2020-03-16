#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 12:51:22 2017
@author: aguimera
"""
import PyGFETdb.DBCore as PyFETdb
from PyGFETdb.DataClass import PyFETPlotDataClass as PyFETplt
import PyGFETdb.DBAnalyze as Dban
import PyGFETdb.DBSearch as DbSearch
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import cm
import matplotlib.colors as colors
import xlsxwriter
import numpy as np
import datetime
import tempfile
import shutil
import sys
from itertools import cycle
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import statsmodels.api as sm

Cortical30Map = {'Shape': (5, 6),
                 'Ch19': (0, 0),
                 'Ch10': (0, 1),     
                 'Ch26': (0, 2),  
                 'Ch02': (0, 3),     
                 'Ch18': (0, 4),   
                 'Ch03': (0, 5), 
                 'Ch09': (1, 0),
                 'Ch31': (1, 1),     
                 'Ch12': (1, 2),  
                 'Ch24': (1, 3),     
                 'Ch04': (1, 4),   
                 'Ch17': (1, 5),
                 'Ch25': (2, 0),
                 'Ch15': (2, 1),     
                 'Ch28': (2, 2),  
                 'Ch08': (2, 3),     
                 'Ch20': (2, 4),   
                 'Ch07': (2, 5), 
                 'Ch11': (3, 0),
                 'Ch29': (3, 1),     
                 'Ch14': (3, 2),  
                 'Ch22': (3, 3),     
                 'Ch05': (3, 4),   
                 'Ch23': (3, 5), 
                 'Ch27': (4, 0),
                 'Ch13': (4, 1),     
                 'Ch30': (4, 2),  
                 'Ch06': (4, 3),     
                 'Ch21': (4, 4),   
                 'Ch01': (4, 5),                       
                 
                 
                 
                 
                 }

Cortical16Map = {'Shape': (4, 4),
                  'Ch01': (2, 3),
                  'Ch02': (1, 2),
                  'Ch03': (3, 3),
                  'Ch04': (0, 2),
                  'Ch05': (0, 3),
                  'Ch06': (3, 2),
                  'Ch07': (1, 3),
                  'Ch08': (2, 2),
                  'Ch09': (2, 0),
                  'Ch10': (1, 1),
                  'Ch11': (3, 0),
                  'Ch12': (0, 1),
                  'Ch13': (0, 0),
                  'Ch14': (3, 1),
                  'Ch15': (1, 0),
                  'Ch16': (2, 1)}



def GetCycleColors(nColors, CMap=cm.jet):
    cmap = cm.ScalarMappable(colors.Normalize(vmin=0, vmax=nColors), CMap)
    col = []
    for i in range(nColors):
        col.append(cmap.to_rgba(i))
    return cycle(col)


def PlotXYLine(Data, Xvar, Yvar, Vgs, Vds, Ud0Norm=True, label=None,
               Ax=None, Color=None, **kwargs):

    fontsize = 'medium'
    labelsize = 5
    scilimits = (-2, 2)

    if Ax is None:
        fig, Ax = plt.subplots()

    for Trtn, Datas in Data.items():
        ValX = np.array([])
        ValY = np.array([])
        for Dat in Datas:
            if Dat.IsOK:
                funcX = Dat.__getattribute__('Get' + Xvar)
                funcY = Dat.__getattribute__('Get' + Yvar)
                try:
                    Valx = funcX(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
                    Valy = funcY(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
                    ValX = np.vstack((ValX, Valx)) if ValX.size else Valx
                    ValY = np.vstack((ValY, Valy)) if ValY.size else Valy
#                    Ax.plot(Valx, Valy, '*', color=Color, label=label)
                except:  # catch *all* exceptions
                    print (Dat.Name, sys.exc_info()[0])

        try:
            if ValX.size == 0:
                continue
            SortIndex = np.argsort(ValX, axis=0)
    #        Ax.plot(ValY[SortIndex][:,:,0], color=Color, label=label)
            Ax.plot(ValX[SortIndex][:, :, 0],
                    ValY[SortIndex][:, :, 0],
                    color=Color,
                    label=label)
        except:
            print (Trtn, sys.exc_info()[0])

    if 'xscale' in kwargs.keys():
        Ax.set_xscale(kwargs['xscale'])
    if 'yscale' in kwargs.keys():
        Ax.set_yscale(kwargs['yscale'])
    else:
        Ax.ticklabel_format(axis='y', style='sci', scilimits=scilimits)

    if 'ylim' in kwargs.keys():
        Ax.set_ylim(kwargs['ylim'])

    Ax.set_ylabel(Yvar, fontsize=fontsize)
    Ax.set_xlabel(Xvar, fontsize=fontsize)
    Ax.tick_params(axis='both', which='Both', labelsize=labelsize)


def CalcParMap(Data, ParMap, ParArgs, ProbeMap):
    Map = np.zeros(ProbeMap['Shape'])*np.NaN
    ret = False
    for Trtn, Dat in Data.items():
        if not Dat[0].IsOK:
            continue
        ch = Trtn.split('-')[-1]
        func = Dat[0].__getattribute__('Get' + ParMap)
        val = func(**ParArgs)
        if val is None:
            continue
        if val.size == 0:
            continue
        if ch not in ProbeMap:
            continue
        Map[ProbeMap[ch]] = val
        ret = True
    if ret:
        return Map
    else:
        return None


class XlsReportBase(object):
    """
    Generic class to define the most common fields and methods to generate XLS
    reports, the key methods and dictionaries are:
        - DBflieds dictionary structures:
            {'VTrts.DCMeas': ('DCMeas', 0, 0)}
            DBtable.Field: (HeaderName, XlsPosition, Not difined)
            - InfoTrtFields, InfoDevFields
    """
    InfoTrtFields = {'Trts.Name': ('Trt Name', 0, 0),
                     'VTrts.DCMeas': ('DCMeas', 1, 0),
                     'VTrts.ACMeas': ('ACMeas', 2, 0),
                     'VTrts.GMeas': ('GMeas', 3, 0),
                     'TrtTypes.Name': ('Trt Type', 4, 0),
                     'TrtTypes.Length': ('Lenght', 5, 0),
                     'TrtTypes.Width': ('Width', 6, 0),
                     'TrtTypes.Pass': ('Pass', 7, 0),
                     'TrtTypes.Area': ('Area', 8, 0),
                     'TrtTypes.Contact': ('Contact', 9, 0),
                     'Trts.Comments': ('T-Comments', 10, 0)}

    InfoDevFields = {'Devices.Name': ('Device', 0, 0),
                     'Devices.Comments': ('D-Comments', 1, 0),
                     'Devices.State': ('D-State', 2, 0),
                     'Devices.ExpOK': ('D-ExpOK', 3, 0),
                     'Wafers.Masks': ('W-Masks', 4, 0),
                     'Wafers.Graphene': ('W-Graphene', 5, 0),
                     'Wafers.Comments': ('W-Comments', 6, 0)}

    InfoDCMeasValues = {'Time': ('Meas Date', 0, {}),
                        'Vds': ('Vds', 1, {}),
                        'Ud0': ('Ud0', 2, {'Vgs': -0.1,
                                           'Vds': None,
                                           'Ud0Norm': True}),
                        'Rds': ('Rds', 3, {'Vgs': -0.1,
                                           'Vds': None,
                                           'Ud0Norm': True}),
                        'GMV': ('GMV', 4, {'Vgs': -0.1,
                                           'Vds': None,
                                           'Ud0Norm': True})}

    InfoACMeasValues = {'Time': ('Meas Date', 0, {}),
                        'Vrms': ('Vrms', 1, {'Vgs': -0.1,
                                             'Vds': None,
                                             'Ud0Norm': True})}

    InfoMeasValues = {'Name': ('Trt Name', 0, {}),
                      'Time': ('Meas Date', 1, {}),
                      'Vds': ('Vds', 2, {'Vgs': -0.1,
                                         'Vds': None,
                                         'Ud0Norm': True}),
                      'Ud0': ('Ud0', 3, {'Vgs': -0.1,
                                         'Vds': None,
                                         'Ud0Norm': True}),
                      'Rds': ('Rds', 4, {'Vgs': -0.1,
                                         'Vds': None,
                                         'Ud0Norm': True}),
                      'GMV': ('GMV', 5, {'Vgs': -0.1,
                                         'Vds': None,
                                         'Ud0Norm': True}),
                      'Vrms': ('Vrms', 6, {'Vgs': -0.1,
                                           'Vds': None,
                                           'Ud0Norm': True})}

    # ('IdTrt',0,0) (Header, position, CountOkDevices() Parameters)
    YeildOK = {'IsOK': ('Working', 0, {'Param' : None,
                                       'RefVal': None,
                                       'Lower': None,
                                       'ParArgs': None}),
               'GMV': ('GMV', 1, {'Param' : 'GMV',
                                  'RefVal': 5e-4,
                                  'Lower': False,
                                  'ParArgs': {'Vgs': -0.1,
                                              'Vds': None,
                                              'Ud0Norm': True}}),
               'Rds': ('Rds', 2, {'Param' : 'Rds',
                                  'RefVal': 10e3,
                                  'Lower': True,
                                  'ParArgs': {'Vgs': -0.1,
                                              'Vds': None,
                                              'Ud0Norm': True}}),
               'Vrms': ('Vrms', 3, {'Param' : 'Vrms',
                                    'RefVal': 100e-6,
                                    'Lower': True,
                                    'ParArgs': {'Vgs': -0.1,
                                                'Vds': None,
                                                'Ud0Norm': True}})}

    ProbeMap = Cortical16Map
    YeildMaps = {'Rds': ({'ParMap': 'Rds',
                          'ParArgs': {'Vgs': -0.1,
                                      'Vds': None,
                                      'Ud0Norm': True},
                          'ProbeMap': ProbeMap
                          },
                          colors.Normalize(200, 1e4),
                          '[Ohms]'),
                 'Vrms': ({'ParMap': 'Vrms',
                           'ParArgs': {'Vgs': -0.1,
                                       'Vds': None,
                                       'Ud0Norm': True},
                           'ProbeMap': ProbeMap
                           },
                           colors.LogNorm(1e-5, 1e-4),
                           '[Vrms]'),
                 'GMV': ({'ParMap': 'GMV',
                           'ParArgs': {'Vgs': -0.1,
                                       'Vds': None,
                                       'Ud0Norm': True},
                           'ProbeMap': ProbeMap
                           },
                           colors.Normalize(vmin=5e-3, vmax=0.5e-3),
                           '[Vrms]')}

    FigsDpi = 150  # Resolution for figures
    DtMax = np.timedelta64(1, 's')

    DictDC = None
    DictAC = None
    DataDC = None
    DataAC = None
    SortList = None

    DevGroups = {}

    def __init__(self, FileName):
        # Init WorkBook, formats
        self.WorkBook = xlsxwriter.Workbook(FileName)
        self.FYeild = self.WorkBook.add_format({'num_format': '0.00%',
                                                'font_color': 'blue',
                                                'bold': True})
        self.Fbold = self.WorkBook.add_format({'bold': True})
        self.FOK = self.WorkBook.add_format({'num_format': '###.00E+2',
                                             'font_color': 'black'})
        self.FNOK = self.WorkBook.add_format({'num_format': '###.00E+2',
                                              'font_color': 'red'})
        # Init temp folder
        self.TmpPath = tempfile.mkdtemp(suffix='PyFET')

        self.DevGroups = {}

        # Init Db connection TODO hide credentials
        self.Mydb = PyFETdb.PyFETdb()

    def GetSortData(self, TrtName):
        Conditions = {'Trts.Name=': (TrtName, )}
        CharTable = 'DCcharacts'
        DatDC, _ = Dban.GetFromDB(Conditions=Conditions,
                                  Table=CharTable,
                                  Last=False,
                                  GetGate=True)
        CharTable = 'ACcharacts'
        DatAC, _ = Dban.GetFromDB(Conditions=Conditions,
                                  Table=CharTable,
                                  Last=False,
                                  GetGate=True)

        DataDC = DatDC.values()[0]
        if len(DatAC.values()) == 0:
            DataAC = []
        else:
            DataAC = DatAC.values()[0]

        DCtimes = []
        for dat in DataDC:
            DCtimes.append(dat.GetTime()[0, 0])
        ACtimes = []
        for dat in DataAC:
            ACtimes.append(dat.GetTime()[0, 0])
        SortList = []
        for idct, DCt in enumerate(DCtimes):
            idacdc = None
            for iact, ACt in enumerate(ACtimes):
                if np.abs(DCt - ACt) < self.DtMax:
                    idacdc = iact
            SortList.append((DCt, idct, idacdc))
        SortList.sort()

        self.DictDC = DatDC
        self.DictAC = DatAC
        self.DataDC = DataDC
        self.DataAC = DataAC
        self.SortList = SortList

    def GetDeviceData(self, DeviceName):
        DeviceNames = (DeviceName, )

        CondBase = {}
        CondBase['Table'] = 'ACcharacts'
        CondBase['Last'] = True
        CondBase['GetGate'] = True
        CondBase['Conditions'] = {'Devices.Name=': DeviceNames,
                                  'CharTable.FuncStep=': ('Report', )}
        Data, _ = Dban.GetFromDB(**CondBase)
        if len(Data) > 0:
            print (DeviceName, 'Getting data from Report Flag')
            self.DevGroups[DeviceName] = CondBase
            return Data, CondBase

        CondBase['Conditions'] = {'Devices.Name=': DeviceNames}
        Data, _ = Dban.GetFromDB(**CondBase)
        if len(Data) > 0:
            print (DeviceName, 'Getting data from last ACcharacts')
            self.DevGroups[DeviceName] = CondBase
            return Data, CondBase

        CondBase['Table'] = 'DCcharacts'
        Data, _ = Dban.GetFromDB(**CondBase)
        if len(Data) > 0:
            print (DeviceName, 'Getting data from last DCcharacts')
            self.DevGroups[DeviceName] = CondBase
            return Data, CondBase

    def GetOKTrts(self, DeviceName=None, Data=None,
                  Param=None, ParArgs=None,
                  RefVal=None, Lower=True):
        Count = 0
        if Data is None:
            Data, _ = self.GetDeviceData(DeviceName)
        for Trtn, Dat in sorted(Data.items()):
            if not Dat[0].IsOK:
                continue

            if Param is 'IsOK' or Param is None:
                Count += 1
                continue

            func = Dat[0].__getattribute__('Get' + Param)
            val = func(**ParArgs)
            if val is None:
                continue
            if hasattr(val, '__iter__'):
                if val.size == 0:
                    continue
            
            # if len(val)>1: #test
            #     val=val[0]
                
            if Lower:
                try:
                    if val < RefVal:
                        Count += 1
                except:
                    continue
            else:
                try:
                    if val > RefVal:
                        Count += 1
                except:
                    continue

        return Count

    def WriteHeaders(self, Sheet, DictList, LocOff=(0, 0),
                     Vertical=True, WriteKeys=False):
        RowOff = LocOff[0]
        ColOff = LocOff[1]
        DictOff = 0
        for Dict in DictList:
            for k, val in Dict.items():
                if Vertical:
                    row = RowOff + val[1] + DictOff
                    col = ColOff
                else:
                    row = RowOff
                    col = ColOff + val[1] + DictOff
                if WriteKeys:
                    header = k
                else:
                    header = val[0]
                Sheet.write(row, col, header, self.Fbold)
            DictOff += len(Dict)

    def WriteDBValues(self, Sheet, DictList, DBSearch, LocOff=(0, 0),
                      Vertical=True, WriteHeader=True, Format=None):

        if WriteHeader:
            self.WriteHeaders(Sheet=Sheet,
                              DictList=DictList,
                              LocOff=LocOff,
                              Vertical=Vertical)
            if Vertical:
                RowOff = LocOff[0]
                ColOff = LocOff[1] + 1
            else:
                RowOff = LocOff[0] + 1
                ColOff = LocOff[1]
        else:
            if Vertical:
                RowOff = LocOff[0]
                ColOff = LocOff[1] + 1
            else:
                RowOff = LocOff[0] + 1
                ColOff = LocOff[1]

        DictOff = 0
        for Dict in DictList:
            DBRes = self.Mydb.GetCharactInfo(Table='DCcharacts',
                                             Conditions=DBSearch,
                                             Output=Dict.keys())
            for ir, Res in enumerate(DBRes):
                for k, val in Res.items():
                    if Vertical:
                        row = RowOff + Dict[k][1] + DictOff
                        col = ColOff + ir
                    else:
                        row = RowOff + ir
                        col = ColOff + Dict[k][1] + DictOff
                    Sheet.write(row, col, val, Format)
            DictOff += len(Dict)

    def WriteMeasValues(self, Sheet, DictList, DataList,
                        LocOff=(0, 0), Vertical=True, Format=None):

        RowOff = LocOff[0]
        ColOff = LocOff[1]
        DictOff = 0

        for idi, Dict in enumerate(DictList):
            if DataList[idi] is None:
                continue
            for par, v in Dict.items():
                if Vertical:
                    row = RowOff + v[1] + DictOff
                    col = ColOff
                else:
                    row = RowOff
                    col = ColOff + v[1] + DictOff

                func = DataList[idi].__getattribute__('Get' + par)
                val = func(**v[2])
                if val is None:
                    continue
                # if hasattr(val, '__iter__'):
                    # print(val, par) #test
                    # if val.size == 0:
                    #     continue
                if par == 'Time':
                    val = val[0, 0].astype(datetime.datetime).strftime('%x %X')
                    Sheet.set_column(col, col, width=20)

                try:
                    Sheet.write(row, col, val, Format)
                except:
                    print ('Error writing', par, val)
            DictOff += len(Dict)

    def WriteTrtMeasHist(self, TrtName, Sheet, Loc, Vertical=False):
        RowOff = Loc[0]
        ColOff = Loc[1]
        DictList = (self.InfoDCMeasValues,
                    self.InfoACMeasValues)

        self.WriteHeaders(Sheet,
                          DictList=DictList,
                          LocOff=Loc,
                          Vertical=Vertical)

        self.GetSortData(TrtName)

        for iMea, SortInd in enumerate(self.SortList):
            DCi = SortInd[1]
            ACi = SortInd[2]

            if self.DataDC[DCi].IsOK:
                Format = self.FOK
            else:
                Format = self.FNOK

            if ACi is None:
                DataList = (self.DataDC[DCi], None)
            else:
                DataList = (self.DataDC[DCi], self.DataAC[ACi])

            self.WriteMeasValues(Sheet=Sheet,
                                 DictList=DictList,
                                 DataList=DataList,
                                 LocOff=(RowOff + iMea + 1, ColOff),
                                 Vertical=Vertical,
                                 Format=Format)

    def WriteDevTrtsMeas(self, Sheet, Data, Loc, Vertical=False):
        #  TODO Check the vertical behivior
        RowOff = Loc[0]
        ColOff = Loc[1]

# Write Header
        self.WriteHeaders(Sheet=Sheet,
                          DictList=(self.InfoMeasValues,
                                    self.InfoTrtFields),
                          LocOff=Loc,
                          Vertical=Vertical)

# Iter for each Trt in Data
        for iTrt, (Trtn, Dat) in enumerate(sorted(Data.items())):
            if Dat[0].IsOK:
                Format = self.FOK
            else:
                Format = self.FNOK
# Fill Meas fields
            row = RowOff + iTrt + 1
            col = ColOff
            self.WriteMeasValues(Sheet=Sheet,
                                 DictList=(self.InfoMeasValues, ),
                                 DataList=(Dat[0], ),
                                 LocOff=(row, col),
                                 Vertical=Vertical,
                                 Format=Format)
            row = RowOff + iTrt
            col = ColOff + len(self.InfoMeasValues)
            self.WriteDBValues(Sheet=Sheet,
                               DictList=(self.InfoTrtFields, ),
                               DBSearch={'Trts.Name=': (Trtn, )},
                               LocOff=(row, col),
                               Vertical=Vertical,
                               WriteHeader=False,
                               Format=Format)

    def WriteOKcount(self, Sheet, Loc, Data, Vertical=True, WriteHeader=True):
        RowOff = Loc[0]
        ColOff = Loc[1]
        TotalCount = float(len(self.ProbeMap)-1)

        if WriteHeader:
            self.WriteHeaders(Sheet=Sheet,
                              DictList=(self.YeildOK,),
                              LocOff=Loc,
                              Vertical=Vertical)

        for v in self.YeildOK.values():
            count = self.GetOKTrts(Data=Data, **v[2])
            if Vertical:
                row = RowOff + v[1]
                col = ColOff + 1
            else:
                row = RowOff + 1
                col = ColOff + v[1]
            Sheet.write(row, col, count)
            if Vertical:
                row = RowOff + v[1]
                col = ColOff + 2
            else:
                row = RowOff + 1
                col = ColOff + v[1] + len(self.YeildOK)
            Sheet.write(row, col, count/TotalCount, self.FYeild)

    def InsertPyGFETplot(self, Sheet, PlotPars, DataDict,
                         ColorOn, Loc=(0, 0), Legend=False):
        Plot = PyFETplt()
        Plot.AddAxes(PlotPars)
        Plot.PlotDataSet(DataDict=DataDict,
                         ColorOn=ColorOn,
                         MarkOn=None)
        if Legend:
            Plot.AddLegend()
        self.InsertFigure(Sheet=Sheet, Loc=Loc)

    def InsertCharMap(self, Sheet, Data, Loc, CalcMapArgs, norm,
                      Units, figsize=(3, 3)):
        sfmt = ticker.ScalarFormatter(useMathText=True)
        sfmt.set_powerlimits((2, 2))

        Map = CalcParMap(Data=Data, **CalcMapArgs)
        if Map is None:
            return

        Fig, Ax = plt.subplots(figsize=figsize)
        Cax = Ax.imshow(Map, cmap=cm.afmhot, norm=norm)
        Ax.set_xlabel('column')
        Ax.set_ylabel('row')
        Ax.set_xticks(np.arange(CalcMapArgs['ProbeMap']['Shape'][1]))
        Ax.set_yticks(np.arange(CalcMapArgs['ProbeMap']['Shape'][0]))
        Ax.set_title(CalcMapArgs['ParMap'])

        cbar = Fig.colorbar(Cax, format=sfmt)
        cbar.set_label(Units, rotation=270, labelpad=10)

        self.InsertFigure(Sheet=Sheet, Loc=Loc, Fig=Fig)

    def InsertFigure(self, Sheet, Loc, Fig=None, Close=True):
        fname = tempfile.mktemp(suffix='.png', dir=self.TmpPath)
        if Fig is None:
            Fig = plt.gcf()
        Fig.savefig(fname, dpi=self.FigsDpi)
        Sheet.insert_image(Loc[0], Loc[1], fname)

        if Close:
            plt.close('all')

    def close(self):
        #  Close and save the XLS file and clear temp folder
        self.WorkBook.close()
        shutil.rmtree(self.TmpPath)


class GenXlsTrtsHistory(XlsReportBase):
    TimePlotParsDC = ('Ids', 'Rds', 'GM', 'Ig')
    TimePlotParsAC = ('Ids', 'GM', 'Vrms', 'Irms')

    TimeEvolPars = {'Ud0': ('DC', {'Vgs': None,
                                   'Vds': None}),
                    'Rds': ('DC', {'Vgs': -0.1,
                                   'Vds': None,
                                   'Ud0Norm': True,
                                   'ylim': (500, 10e3)}),
                    'Ids': ('DC', {'Vgs': -0.1,
                                   'Vds': None,
                                   'Ud0Norm': True}),
                    'GM': ('DC', {'Vgs': -0.1,
                                  'Vds': None,
                                  'Ud0Norm': True,
                                  'ylim': (-5e-4, 0)}),
                    'Vrms': ('AC', {'Vgs': -0.1,
                                    'Vds': None,
                                    'Ud0Norm': True,
                                    'yscale': 'log',
                                    'ylim': (1e-5, 2e-4)})}

    def __init__(self, FileName, Conditions):
        super(GenXlsTrtsHistory, self).__init__(FileName=FileName)

        self.dbConditions = Conditions
        GroupBy = 'Trts.Name'
        self.TrtsList = Dban.FindCommonValues(Table='DCcharacts',
                                              Parameter=GroupBy,
                                              Conditions=Conditions)
        GroupBy = 'Devices.Name'
        self.DeviceName = Dban.FindCommonValues(Table='DCcharacts',
                                                Parameter=GroupBy,
                                                Conditions=Conditions)

        self.WorkBook.add_worksheet('Summary')
        for TrtName in sorted(self.TrtsList):
            self.WorkBook.add_worksheet(TrtName)

    def GenFullReport(self):
        #  Init Summary plot
        self.SummaryFigCol = GetCycleColors(len(self.TrtsList))
        self.SummaryFig, self.SummaryAx = plt.subplots(len(self.TimeEvolPars),
                                                       1,
                                                       sharex=True,
                                                       figsize=(12, 12))

        for TrtName in self.TrtsList:
            self.GenTrtReport(TrtName)

#  Genertate Summary sheet
        Sheet = self.WorkBook.sheetnames['Summary']
        self.WriteDBValues(Sheet,
                           LocOff=(0, 0),
                           DictList=(self.InfoTrtFields, ),
                           DBSearch=self.dbConditions,
                           Vertical=False)
#  Insert summary plot generated updated at GenTrtReport
        self.SummaryFig.tight_layout()
        self.SummaryFig.subplots_adjust(hspace=0)
        self.InsertFigure(Sheet=Sheet, Fig=self.SummaryFig, Loc=(0, 15))

    def GenTrtReport(self, TrtName):
        Sheet = self.WorkBook.sheetnames[TrtName]
        self.WriteDBValues(Sheet,
                           LocOff=(0, 0),
                           DictList=(self.InfoTrtFields, ),
                           DBSearch={'Trts.Name=': (TrtName, )})

        self.WriteTrtMeasHist(TrtName=TrtName,
                              Sheet=Sheet,
                              Loc=(len(self.InfoTrtFields)+2, 0))

        Loc = (0, len(self.InfoDCMeasValues) + len(self.InfoACMeasValues) + 1)
        self.InsertPyGFETplot(Sheet=Sheet,
                              Loc=Loc,
                              PlotPars=self.TimePlotParsDC,
                              DataDict=self.DictDC,
                              ColorOn='Date')

        Loc = (30, len(self.InfoDCMeasValues) + len(self.InfoACMeasValues) + 1)
        self.InsertPyGFETplot(Sheet=Sheet,
                              Loc=Loc,
                              PlotPars=self.TimePlotParsAC,
                              DataDict=self.DictAC,
                              ColorOn='Date')

#TODO fix this part in a generic fucntion
# Generate time evolution plots
        Fig, Ax = plt.subplots(len(self.TimeEvolPars), 1, sharex=True)
        for i, (k, v) in enumerate(self.TimeEvolPars.items()):
            if v[0] == 'DC':
                dat = self.DictDC
            elif v[0] == 'AC':
                dat = self.DictAC
            
            Dban.PlotXYVars(Data=dat,
                            Xvar='DateTime',
                            Yvar=k,
                            Ax=Ax[i],
                            **v[1])
# Update summary plot
            PlotXYLine(Data=dat,
                       Xvar='Time',
                       Yvar=k,
                       Ax=self.SummaryAx[i],
                       Color=self.SummaryFigCol.next(),
                       **v[1])

        Fig.tight_layout()
        Fig.subplots_adjust(hspace=0)
# insert time evolution plots
        self.InsertFigure(Sheet=Sheet, Fig=Fig, Loc=(60, 10))


class GenXlsReport(XlsReportBase):
    
    CharPars = ('Ids', 'GM', 'Vrms', 'Rds')

    SummaryPlotSpacing = 25
    SummaryLinePlots = ({'Xvar': 'Vgs', 'Yvar': 'Ids',
                         'Vds': None, 'PlotOverlap': False, 'Ud0Norm': False},
                        {'Xvar': 'Vgs', 'Yvar': 'GM',
                        'Vds': None, 'PlotOverlap': False, 'Ud0Norm': False},
                        {'Xvar': 'Vgs', 'Yvar': 'Vrms',
                         'Vds': None, 'PlotOverlap': False, 'Ud0Norm': True,
                         'yscale': 'log'})

    SummaryBoxPlots = ({'Plot': True, 'Boxplot': False,
                        'Param': 'Vrms', 'Vgs': -0.1,
                        'Ud0Norm': True, 'Vds': None, 'yscale': 'log'},
                       {'Plot': True, 'Boxplot': True,
                        'Param': 'Ud0', 'Vgs': -0.1,
                        'Ud0Norm': True, 'Vds': None})

    def __init__(self, FileName, Conditions, PathHistory=None):
        super(GenXlsReport, self).__init__(FileName=FileName)
        
# Find devicelist that will be reported
        self.dbConditions = Conditions
        GroupBy = 'Devices.Name'        
        self.DevicesList = Dban.FindCommonValues(Table='DCcharacts',
                                                 Parameter=GroupBy,
                                                 Conditions=Conditions)
#  Add Worksheets
        self.WorkBook.add_worksheet('Summary')
        for DevName in sorted(self.DevicesList):
            self.WorkBook.add_worksheet(DevName)
#            if PathHistory is not None:
#                Cond = {'Devices.Name=': (DevName, )}
#                XlsHist = GenXlsDeviceHistory(PathHistory + DevName +'.xlsx', Cond)
#                XlsHist.GenFullReport()
#                XlsHist.close()

    def GenFullReport(self):
        Sheet = self.WorkBook.sheetnames['Summary']
        self.WriteHeaders(Sheet=Sheet,
                          LocOff=(0, 0),
                          DictList=(self.InfoDevFields,
                                    self.YeildOK,
                                    self.YeildOK),
                          Vertical=False)

        for idev, DevName in enumerate(sorted(self.DevicesList)):
            Data = self.GenDeviceReport(DevName)
            self.WriteOKcount(Sheet=Sheet,
                              Data=Data,
                              Loc=(idev, len(self.InfoDevFields)),
                              Vertical=False,
                              WriteHeader=False)
            self.WriteDBValues(Sheet,
                               LocOff=(idev, 0),
                               DictList=(self.InfoDevFields, ),
                               DBSearch={'Devices.Name=': (DevName, )},
                               WriteHeader=False,
                               Vertical=False)

        Sheet.set_column(0, len(self.InfoDevFields), width=12)

        for ipl, LiPlots in enumerate(self.SummaryLinePlots):
            Fig, _ = Dban.SearchAndPlot(Groups=self.DevGroups, **LiPlots)
            self.InsertFigure(Sheet=Sheet,
                              Loc=(idev + 2 + ipl*self.SummaryPlotSpacing, 7),
                              Fig=Fig)
        for ipl, BoxPlots in enumerate(self.SummaryBoxPlots):
            Dban.SearchAndGetParam(Groups=self.DevGroups, **BoxPlots)
            self.InsertFigure(Sheet=Sheet,
                              Loc=(idev + 2 + ipl*self.SummaryPlotSpacing, 0))

    def GenDeviceReport(self, DeviceName):
        Sheet = self.WorkBook.sheetnames[DeviceName]
        Sheet.set_column(0, 0, width=20)
        self.WriteDBValues(Sheet,
                           LocOff=(0, 0),
                           DictList=(self.InfoDevFields, ),
                           DBSearch={'Devices.Name=': (DeviceName, )})

        Data, _ = self.GetDeviceData(DeviceName)
        self.WriteOKcount(Sheet=Sheet,
                          Loc=(len(self.InfoDevFields), 0),
                          Data=Data)
        Loc = (15, 0)
        self.WriteDevTrtsMeas(Sheet, Data, Loc)

        Loc = (Loc[0] + len(Data) + 2,
               0)
        self.InsertPyGFETplot(Sheet=Sheet,
                              Loc=Loc,
                              PlotPars=self.CharPars,
                              DataDict=Data,
                              ColorOn='Trt',
                              Legend=True)

        for i, v in enumerate(self.YeildMaps.values()):
            self.InsertCharMap(Sheet=Sheet,
                               Data=Data,
                               Loc=(0, 5+i*5),
                               CalcMapArgs=v[0],
                               norm=v[1],
                               Units=v[2])
        return Data


class FittingReport(object):
    XVar = 'IonStrength'
    XVarLog = True
    YVar = 'Ud0'
    YVarLog = False
    GroupBase = None

    figsize=(4, 4)

    InfoMeasValues = {'Name': ('Trt Name', 0, {}),
                      'Time': ('Meas Date', 1, {}),
                      'Vds': ('Vds', 2, {'Vgs': -0.1,
                                         'Vds': None,
                                         'Ud0Norm': True}),
                      'Ud0': ('Ud0', 3, {'Vgs': -0.1,
                                         'Vds': None,
                                         'Ud0Norm': True}),
                      'IonStrength': ('IonStrength', 4, {}),
                      'Ph': ('Ph', 5, {}),
                      'FuncStep': ('FuncStep', 6, {}),
                      'Comments': ('Comments', 7, {})}

    FittingValues = {'R2': ('rsquared', 0, None),
                     'Slope': ('params', 1, 1),
                     'OffSet': ('params', 2, 0)}

    def WriteFittingResults(self, Sheet, Loc, WriteHeaders=True):
        if WriteHeaders:
            self.XlsRep.WriteHeaders(Sheet,
                                     DictList=(self.FittingValues, ),
                                     LocOff=Loc,
                                     Vertical=False,
                                     WriteKeys=True)
        row = Loc[0]+1
        for v in self.FittingValues.values():
            try:
                val = self.Res.__getattribute__(v[0])
                if v[2] is not None:
                    val = val[v[2]]
                col = Loc[1]+v[1]
                Sheet.write(row, col, val)
            except:
                print ('Error in fitting write')

    def __init__(self, XVar, XVarLog, YVar, YVarLog, GroupBase, XlsRep):
        self.XVar = XVar
        self.XVarLog = XVarLog
        self.YVar = YVar
        self.YVarLog = YVarLog
        self.GroupBase = GroupBase
        self.XlsRep = XlsRep

    def GenDebugPlot(self, Sheet, Loc):
        fig, Ax = plt.subplots(figsize=self.figsize)
        fig, _ = Dban.PlotGroupBy(GroupBase=self.GroupBase,
                                  GroupBy='CharTable.{}'.format(self.XVar),
                                  Xvar='Vgs',
                                  Yvar='Ids',
                                  PlotOverlap=True,
                                  Ud0Norm=False,
                                  Ax=Ax)
        self.XlsRep.InsertFigure(Sheet, Loc)

    def CalcLinearFitting(self, Sheet, Loc):
        self.XlsRep.WriteHeaders(Sheet,
                                 DictList=(self.InfoMeasValues, ),
                                 LocOff=Loc,
                                 Vertical=False)
        RowOff = Loc[0]
        ColOff = Loc[1]
        col = ColOff + len(self.InfoMeasValues)
        Sheet.write(RowOff, col, self.XVar, self.XlsRep.Fbold)
        Sheet.write(RowOff, col+1, self.YVar, self.XlsRep.Fbold)

        Grs = DbSearch.GenGroups(GroupBase=self.GroupBase,
                                 GroupBy='CharTable.{}'.format(self.XVar))
        DatFit = []
        for gr in Grs.values():
            Dat, Trts = DbSearch.GetFromDB(**gr)
            for dat in Dat[Trts[0]]:
                DatFit.append(dat)

        ValY = np.array([])
        ValX = np.array([])
        for ip, dat in enumerate(DatFit):
            row = RowOff + ip + 1
            self.XlsRep.WriteMeasValues(Sheet,
                                        DictList=(self.InfoMeasValues, ),
                                        DataList=(dat, ),
                                        LocOff=(row, ColOff),
                                        Vertical=False)
            funcx = dat.__getattribute__('Get' + self.XVar)
            funcy = dat.__getattribute__('Get' + self.YVar)
            valx = funcx()
            valy = funcy()
            ValY = np.vstack((ValY, valy)) if ValY.size else valy
            ValX = np.vstack((ValX, valx)) if ValX.size else valx
            Sheet.write(row, col, valx, self.XlsRep.Fbold)
            Sheet.write(row, col+1, valy, self.XlsRep.Fbold)

        si = np.argsort(ValX[:, 0])
        ValX = ValX[si, 0]
        ValY = ValY[si, 0]

        if self.XVarLog:
            ValX = np.log10(ValX)

        if self.YVarLog:
            ValY = np.log10(ValY)

        plt.figure(figsize=self.figsize)
        plt.plot(ValX, ValY, '*')

        X = sm.add_constant(ValX)
        Res = sm.OLS(ValY, X).fit()
        self.Res = Res
        prstd, iv_l, iv_u = wls_prediction_std(Res)
        plt.plot(ValX, Res.fittedvalues, 'k--')
        plt.fill_between(ValX, iv_u, iv_l,
                         color='b',
                         linewidth=0.0,
                         alpha=0.3)
        self.XlsRep.InsertFigure(Sheet, Loc=(RowOff+2, col+3))
        self.WriteFittingResults(Sheet, Loc=(RowOff, col+3))


class GenXlsFittingReport(XlsReportBase):
    TimePlotParsDC = ('Ids', 'Rds', 'GM', 'Ig')
    TimeEvolPars = {'Ud0': ('DC', {'Vgs': None,
                                   'Vds': None})}

    CalTypes = {'pHCal': {'XVar': 'Ph',
                          'XVarLog': False,
                          'YVar': 'Ud0',
                          'YVarLog': False},
                'IonCal': {'XVar': 'IonStrength',
                           'XVarLog': True,
                           'YVar': 'Ud0',
                           'YVarLog': False},
                'Tromb': {'XVar': 'AnalyteCon',
                          'XVarLog': True,
                          'YVar': 'Ud0',
                          'YVarLog': False}}

    DBCalField = 'CharTable.Comments'
    DBCalFlags = ('%Cal%', )

    CalMaps = None
    CalMapVars = None
    figsize = (4, 4)

    def __init__(self, FileName, GroupBase,
                 DBCalField='CharTable.Comments', DBCalFlags=('%Cal%', )):

        super(GenXlsFittingReport, self).__init__(FileName=FileName)

        self.InfoDCMeasValues.update({'IonStrength': ('IonStrength', 5, {})})
        self.InfoDCMeasValues.update({'Ph': ('pH', 6, {})})
        self.InfoDCMeasValues.update({'FuncStep': ('FuncStep', 7, {})})
        self.InfoDCMeasValues.update({'AnalyteCon': ('AnalyteCon', 8, {})})
        self.InfoDCMeasValues.update({'Comments': ('Comments', 9, {})})

        self.DBCalField = DBCalField
        self.DBCalFlags = DBCalFlags
        self.GroupBase = GroupBase

        GroupBy = 'Trts.Name'
        self.TrtsList = Dban.FindCommonValues(Parameter=GroupBy,
                                              **GroupBase)

        self.WorkBook.add_worksheet('Summary')
        for TrtName in sorted(self.TrtsList):
            self.WorkBook.add_worksheet(TrtName)

    def GenFullReport(self):
        SheetSummary = self.WorkBook.sheetnames['Summary']
        self.WriteHeaders(SheetSummary,
                          DictList=(self.InfoTrtFields, ),
                          Vertical=False)
        for it, TrtName in enumerate(sorted(self.TrtsList)):
            self.WriteDBValues(SheetSummary,
                               LocOff=(it, 0),
                               DictList=(self.InfoTrtFields, ),
                               Vertical=False,
                               DBSearch={'Trts.Name=': (TrtName, )},
                               WriteHeader=False)
            self.GenTrtReport(TrtName, it)

    def GenTrtReport(self, TrtName, itrt):
        Sheet = self.WorkBook.sheetnames[TrtName]
        SheetSummary = self.WorkBook.sheetnames['Summary']

        self.WriteDBValues(Sheet,
                           DictList=(self.InfoTrtFields,),
                           DBSearch={'Trts.Name=': (TrtName, )})

        GrBase = self.GroupBase.copy()
        GrBase['Conditions'].update({'Trts.Name=': (TrtName, )})
        CalFlag = {'{} like'.format(self.DBCalField): self.DBCalFlags}
        GrBase['Conditions'].update(CalFlag)
        CalList = DbSearch.FindCommonValues(Parameter=self.DBCalField,
                                            **GrBase)

        ic = 0
        row = len(self.InfoTrtFields) + 10
        for cal in CalList:
            caltype = cal.split(' ')[0]
            if caltype in self.CalTypes:
                GrBase = self.GroupBase.copy()
                GrBase['Conditions'].update({'Trts.Name=': (TrtName,)})
                CalFlag = {'{} like'.format(self.DBCalField): (cal, )}
                GrBase['Conditions'].update(CalFlag)
                FitRep = FittingReport(XlsRep=self,
                                       GroupBase=GrBase,
                                       **self.CalTypes[caltype])
                Loc = (row+2+25*ic, len(FitRep.InfoMeasValues)+10)
                FitRep.GenDebugPlot(Sheet, Loc=Loc)
                FitRep.CalcLinearFitting(Sheet, Loc=(row+25*ic, 0))
                srow = itrt + 1
                scol = ic*(len(FitRep.FittingValues)+1) + len(self.InfoTrtFields)
                SheetSummary.write(srow, scol, cal, self.Fbold)
                FitRep.WriteFittingResults(SheetSummary,
                                           (srow-1, scol+1),
                                           WriteHeaders=False)
                ic += 1
            else:
                print (cal, caltype, 'Not found')


        row = len(self.InfoTrtFields) + 12 + len(CalList)*25
        self.WriteTrtMeasHist(TrtName=TrtName,
                              Sheet=Sheet,
                              Loc=(row, 0))
        
        Loc = (row + 25,
               len(self.InfoDCMeasValues) + len(self.InfoACMeasValues) + 1)
        self.InsertPyGFETplot(Sheet=Sheet,
                              Loc=Loc,
                              PlotPars=self.TimePlotParsDC,
                              DataDict=self.DictDC,
                              ColorOn='Date')

#TODO fix this part in a generic fucntion
# Generate time evolution plots
        Loc = (row,
               len(self.InfoDCMeasValues) + len(self.InfoACMeasValues) + 1)
        Fig, Ax = plt.subplots(len(self.TimeEvolPars), 1, sharex=True)
        if len(self.TimeEvolPars) == 1:
            Ax = [Ax,]
        for i, (k, v) in enumerate(self.TimeEvolPars.items()):
            if v[0] == 'DC':
                dat = self.DictDC
            elif v[0] == 'AC':
                dat = self.DictAC
            
            Dban.PlotXYVars(Data=dat,
                            Xvar='FuncStep',
                            Yvar=k,
                            Ax=Ax[i],
                            **v[1])
        Fig.tight_layout()
        Fig.subplots_adjust(hspace=0)
# insert time evolution plots
        self.InsertFigure(Sheet=Sheet, Fig=Fig, Loc=Loc)

        
        if self.CalMaps is not None:
            for ic, calmap in enumerate(self.CalMaps):
                self.InsertCalMap(Sheet,
                                  (0, ic*6 + 6),
                                  TrtName,
                                  calmap)

    def InsertCalMap(self, Sheet, Loc, TrtName, CalMap):

        Gr = self.GroupBase.copy()
        Gr['Conditions'].update({'Trts.Name=': (TrtName,)})
        CalFlag = {'{} like'.format(self.DBCalField): CalMap}
        Gr['Conditions'].update(CalFlag)

        dat, _ = Dban.GetFromDB(**Gr)
        Dats = dat[TrtName]

        x = np.ones([len(Dats)])
        y = np.ones([len(Dats)])
        z = np.ones([len(Dats)])
        for i, d in enumerate(Dats):
            funcx = d.__getattribute__('Get' + self.CalMapVars['XVar'])
            funcy = d.__getattribute__('Get' + self.CalMapVars['YVar'])
            funcz = d.__getattribute__('Get' + self.CalMapVars['ZVar'])
            x[i] = funcx()
            y[i] = funcy()
            z[i] = funcz()

        if self.CalMapVars['YVarLog']:
            y = np.log10(y)
        if self.CalMapVars['XVarLog']:
            x = np.log10(x)

        plt.figure(figsize=self.figsize)
        plt.tricontourf(x, y, z, 100, cmap='seismic')
        plt.plot(x, y, 'ko')
        plt.colorbar()
        plt.title(CalMap[0]+CalMap[1])

        self.InsertFigure(Sheet, Loc)


class XlsGraphPadPrism(XlsReportBase):
    TimePlotParsDC = ('Ids', 'Rds', 'GM', 'Ig')
    TimeEvolPars = {'Ud0': ('DC', {'Vgs': None,
                                   'Vds': None})}

    CalTypes = {'pHCal': {'XVar': 'Ph',
                          'XVarLog': False,
                          'YVar': 'Ud0',
                          'YVarLog': False},
                'IonCal': {'XVar': 'IonStrength',
                           'XVarLog': True,
                           'YVar': 'Ud0',
                           'YVarLog': False},
                'Tromb': {'XVar': 'AnalyteCon',
                          'XVarLog': True,
                          'YVar': 'Ud0',
                          'YVarLog': False}}

    DBCalField = 'CharTable.Comments'
    DBCalFlags = ('%Cal%', )

    CalMaps = None
    CalMapVars = None
    figsize = (4, 4)

    def __init__(self, FileName, GroupBase,
                 DBCalField='CharTable.Comments', DBCalFlags=('%Cal%', )):

        super(XlsGraphPadPrism, self).__init__(FileName=FileName)

        self.InfoDCMeasValues.update({'IonStrength': ('IonStrength', 5, {})})
        self.InfoDCMeasValues.update({'Ph': ('pH', 6, {})})
        self.InfoDCMeasValues.update({'FuncStep': ('FuncStep', 7, {})})
        self.InfoDCMeasValues.update({'AnalyteCon': ('AnalyteCon', 8, {})})
        self.InfoDCMeasValues.update({'Comments': ('Comments', 9, {})})

        self.DBCalField = DBCalField
        self.DBCalFlags = DBCalFlags
        self.GroupBase = GroupBase

        GroupBy = 'Trts.Name'
        self.TrtsList = Dban.FindCommonValues(Parameter=GroupBy,
                                              **GroupBase)

        self.WorkBook.add_worksheet('GraphPadPrism')
        self.WorkBook.add_worksheet('Summary')
        for TrtName in sorted(self.TrtsList):
            self.WorkBook.add_worksheet(TrtName)

    def GenFullReport(self):
        self.GenGraphPadPrism(self.WorkBook.sheetnames['GraphPadPrism'])
        SheetSummary = self.WorkBook.sheetnames['Summary']
        self.WriteHeaders(SheetSummary,
                          DictList=(self.InfoTrtFields, ),
                          Vertical=False)
        for it, TrtName in enumerate(sorted(self.TrtsList)):
            self.WriteDBValues(SheetSummary,
                               LocOff=(it, 0),
                               DictList=(self.InfoTrtFields, ),
                               Vertical=False,
                               DBSearch={'Trts.Name=': (TrtName, )},
                               WriteHeader=False)
            self.GenTrtReport(TrtName, it)       

    def GenGraphPadPrism(self, Sheet, Loc=(0, 0)):        
        RowOff = Loc[0]
        ColOff = Loc[1]
        
        # Gen trt Headers
        Trts = DbSearch.FindCommonValues('Trts.Name', **self.GroupBase)        
        TrtLoc = {}
        for itrt, trt in enumerate(sorted(Trts)):
            col = itrt + ColOff + 2
            TrtLoc.update({trt: col })
            Sheet.write(RowOff, col, trt, self.Fbold)
    
        row = RowOff
        GrFunc = DbSearch.GenGroups(self.GroupBase.copy(),
                                    'CharTable.FuncStep',
                                    LongName=False)
        for GrfN, Grf in GrFunc.items():
            GrAnalyte = DbSearch.GenGroups(Grf.copy(),
                                          'CharTable.AnalyteCon',
                                          LongName=False)
            for GrAn, GrA in GrAnalyte.items():
                row += 1
                Sheet.write(row, 0, GrfN, self.Fbold)
                Sheet.write(row, 1, GrAn, self.Fbold)
                DatDict, _ = DbSearch.GetFromDB(**GrA)
                for Dat in DatDict.values():
                    func = Dat[0].__getattribute__('Get' + 'Ud0')
                    Sheet.write(row, TrtLoc[Dat[0].Name],
                                func())
                    

    def GenTrtReport(self, TrtName, itrt):
        Sheet = self.WorkBook.sheetnames[TrtName]
        SheetSummary = self.WorkBook.sheetnames['Summary']

        self.WriteDBValues(Sheet,
                           DictList=(self.InfoTrtFields,),
                           DBSearch={'Trts.Name=': (TrtName, )})

        GrBase = self.GroupBase.copy()
        GrBase['Conditions'].update({'Trts.Name=': (TrtName, )})
        CalFlag = {'{} like'.format(self.DBCalField): self.DBCalFlags}
        GrBase['Conditions'].update(CalFlag)
        CalList = DbSearch.FindCommonValues(Parameter=self.DBCalField,
                                            **GrBase)

        ic = 0
        row = len(self.InfoTrtFields) + 10
        for cal in CalList:
            caltype = cal.split(' ')[0]
            if caltype in self.CalTypes:
                GrBase = self.GroupBase.copy()
                GrBase['Conditions'].update({'Trts.Name=': (TrtName,)})
                CalFlag = {'{} like'.format(self.DBCalField): (cal, )}
                GrBase['Conditions'].update(CalFlag)
                FitRep = FittingReport(XlsRep=self,
                                       GroupBase=GrBase,
                                       **self.CalTypes[caltype])
                Loc = (row+2+25*ic, len(FitRep.InfoMeasValues)+10)
                FitRep.GenDebugPlot(Sheet, Loc=Loc)
                FitRep.CalcLinearFitting(Sheet, Loc=(row+25*ic, 0))
                srow = itrt + 1
                scol = ic*(len(FitRep.FittingValues)+1) + len(self.InfoTrtFields)
                SheetSummary.write(srow, scol, cal, self.Fbold)
                FitRep.WriteFittingResults(SheetSummary,
                                           (srow-1, scol+1),
                                           WriteHeaders=False)
                ic += 1
            else:
                print (cal, caltype, 'Not found')


        row = len(self.InfoTrtFields) + 12 + len(CalList)*25
        self.WriteTrtMeasHist(TrtName=TrtName,
                              Sheet=Sheet,
                              Loc=(row, 0))
        
        Loc = (row + 25,
               len(self.InfoDCMeasValues) + len(self.InfoACMeasValues) + 1)
        self.InsertPyGFETplot(Sheet=Sheet,
                              Loc=Loc,
                              PlotPars=self.TimePlotParsDC,
                              DataDict=self.DictDC,
                              ColorOn='Date')

#TODO fix this part in a generic fucntion
# Generate time evolution plots
        Loc = (row,
               len(self.InfoDCMeasValues) + len(self.InfoACMeasValues) + 1)
        Fig, Ax = plt.subplots(len(self.TimeEvolPars), 1, sharex=True)
        if len(self.TimeEvolPars) == 1:
            Ax = [Ax,]
        for i, (k, v) in enumerate(self.TimeEvolPars.items()):
            if v[0] == 'DC':
                dat = self.DictDC
            elif v[0] == 'AC':
                dat = self.DictAC
            
            Dban.PlotXYVars(Data=dat,
                            Xvar='FuncStep',
                            Yvar=k,
                            Ax=Ax[i],
                            **v[1])
        Fig.tight_layout()
        Fig.subplots_adjust(hspace=0)
# insert time evolution plots
        self.InsertFigure(Sheet=Sheet, Fig=Fig, Loc=Loc)

        
        if self.CalMaps is not None:
            for ic, calmap in enumerate(self.CalMaps):
                self.InsertCalMap(Sheet,
                                  (0, ic*6 + 6),
                                  TrtName,
                                  calmap)

    def InsertCalMap(self, Sheet, Loc, TrtName, CalMap):

        Gr = self.GroupBase.copy()
        Gr['Conditions'].update({'Trts.Name=': (TrtName,)})
        CalFlag = {'{} like'.format(self.DBCalField): CalMap}
        Gr['Conditions'].update(CalFlag)

        dat, _ = Dban.GetFromDB(**Gr)
        Dats = dat[TrtName]

        x = np.ones([len(Dats)])
        y = np.ones([len(Dats)])
        z = np.ones([len(Dats)])
        for i, d in enumerate(Dats):
            funcx = d.__getattribute__('Get' + self.CalMapVars['XVar'])
            funcy = d.__getattribute__('Get' + self.CalMapVars['YVar'])
            funcz = d.__getattribute__('Get' + self.CalMapVars['ZVar'])
            x[i] = funcx()
            y[i] = funcy()
            z[i] = funcz()

        if self.CalMapVars['YVarLog']:
            y = np.log10(y)
        if self.CalMapVars['XVarLog']:
            x = np.log10(x)

        plt.figure(figsize=self.figsize)
        plt.tricontourf(x, y, z, 100, cmap='seismic')
        plt.plot(x, y, 'ko')
        plt.colorbar()
        plt.title(CalMap[0]+CalMap[1])

        self.InsertFigure(Sheet, Loc)


class GenXlsFittingReport1(XlsReportBase):

    CalTypes = {'pHCal': {'XVar': 'Ph',
                          'XVarLog': False,
                          'YVar': 'Ud0',
                          'YVarLog': False},
                'IonCal': {'XVar': 'IonStrength',
                           'XVarLog': True,
                           'YVar': 'Ud0',
                           'YVarLog': False},
                'Tromb': {'XVar': 'AnalyteCon',
                          'XVarLog': True,
                          'YVar': 'Ud0',
                          'YVarLog': False}}

    DBCalField = 'CharTable.Comments'
    DBCalFlags = ('%Cal%', )

    CalMaps = None
    CalMapVars = None
    figsize = (4, 4)

    def __init__(self, FileName, GroupBase,
                 DBCalField='CharTable.Comments', DBCalFlags=('%Cal%', )):

        super(GenXlsFittingReport, self).__init__(FileName=FileName)

        self.InfoDCMeasValues.update({'IonStrength': ('IonStrength', 5, {})})
        self.InfoDCMeasValues.update({'Ph': ('pH', 6, {})})
        self.InfoDCMeasValues.update({'FuncStep': ('FuncStep', 7, {})})
        self.InfoDCMeasValues.update({'AnalyteCon': ('AnalyteCon', 8, {})})
        self.InfoDCMeasValues.update({'Comments': ('Comments', 9, {})})

        self.DBCalField = DBCalField
        self.DBCalFlags = DBCalFlags
        self.GroupBase = GroupBase

        GroupBy = 'Trts.Name'
        self.TrtsList = Dban.FindCommonValues(Parameter=GroupBy,
                                              **GroupBase)

        self.WorkBook.add_worksheet('Summary')
        for TrtName in sorted(self.TrtsList):
            self.WorkBook.add_worksheet(TrtName)

    def GenFullReport(self):
        SheetSummary = self.WorkBook.sheetnames['Summary']
        self.WriteHeaders(SheetSummary,
                          DictList=(self.InfoTrtFields, ),
                          Vertical=False)
        for it, TrtName in enumerate(sorted(self.TrtsList)):
            self.WriteDBValues(SheetSummary,
                               LocOff=(it, 0),
                               DictList=(self.InfoTrtFields, ),
                               Vertical=False,
                               DBSearch={'Trts.Name=': (TrtName, )},
                               WriteHeader=False)
            self.GenTrtReport(TrtName, it)

    def GenTrtReport(self, TrtName, itrt):
        Sheet = self.WorkBook.sheetnames[TrtName]
        SheetSummary = self.WorkBook.sheetnames['Summary']

        self.WriteDBValues(Sheet,
                           DictList=(self.InfoTrtFields,),
                           DBSearch={'Trts.Name=': (TrtName, )})

        GrBase = self.GroupBase.copy()
        GrBase['Conditions'].update({'Trts.Name=': (TrtName, )})
        CalFlag = {'{} like'.format(self.DBCalField): self.DBCalFlags}
        GrBase['Conditions'].update(CalFlag)
        CalList = DbSearch.FindCommonValues(Parameter=self.DBCalField,
                                            **GrBase)

        ic = 0
        row = len(self.InfoTrtFields) + 10
        for cal in CalList:
            caltype = cal.split(' ')[0]
            if caltype in self.CalTypes:
                GrBase = self.GroupBase.copy()
                GrBase['Conditions'].update({'Trts.Name=': (TrtName,)})
                CalFlag = {'{} like'.format(self.DBCalField): (cal, )}
                GrBase['Conditions'].update(CalFlag)
                FitRep = FittingReport(XlsRep=self,
                                       GroupBase=GrBase,
                                       **self.CalTypes[caltype])
                Loc = (row+2+25*ic, len(FitRep.InfoMeasValues)+10)
                FitRep.GenDebugPlot(Sheet, Loc=Loc)
                FitRep.CalcLinearFitting(Sheet, Loc=(row+25*ic, 0))
                srow = itrt + 1
                scol = ic*(len(FitRep.FittingValues)+1) + len(self.InfoTrtFields)
                SheetSummary.write(srow, scol, cal, self.Fbold)
                FitRep.WriteFittingResults(SheetSummary,
                                           (srow-1, scol+1),
                                           WriteHeaders=False)
                ic += 1
            else:
                print (cal, caltype, 'Not found')


        row = len(self.InfoTrtFields) + 12 + len(CalList)*25
        self.WriteTrtMeasHist(TrtName=TrtName,
                              Sheet=Sheet,
                              Loc=(row, 0))

        if self.CalMaps is not None:
            for ic, calmap in enumerate(self.CalMaps):
                self.InsertCalMap(Sheet,
                                  (0, ic*6 + 6),
                                  TrtName,
                                  calmap)

    def InsertCalMap(self, Sheet, Loc, TrtName, CalMap):

        Gr = self.GroupBase.copy()
        Gr['Conditions'].update({'Trts.Name=': (TrtName,)})
        CalFlag = {'{} like'.format(self.DBCalField): CalMap}
        Gr['Conditions'].update(CalFlag)

        dat, _ = Dban.GetFromDB(**Gr)
        Dats = dat[TrtName]

        x = np.ones([len(Dats)])
        y = np.ones([len(Dats)])
        z = np.ones([len(Dats)])
        for i, d in enumerate(Dats):
            funcx = d.__getattribute__('Get' + self.CalMapVars['XVar'])
            funcy = d.__getattribute__('Get' + self.CalMapVars['YVar'])
            funcz = d.__getattribute__('Get' + self.CalMapVars['ZVar'])
            x[i] = funcx()
            y[i] = funcy()
            z[i] = funcz()

        if self.CalMapVars['YVarLog']:
            y = np.log10(y)
        if self.CalMapVars['XVarLog']:
            x = np.log10(x)

        plt.figure(figsize=self.figsize)
        plt.tricontourf(x, y, z, 100, cmap='seismic')
        plt.plot(x, y, 'ko')
        plt.colorbar()
        plt.title(CalMap[0]+CalMap[1])

        self.InsertFigure(Sheet, Loc)

#if __name__ == "__main__":
#    plt.ioff()
#
#    plt.close('all')   
#    
#    CharTable = 'DCcharacts'
#    DeviceNames = ('B10803W17-Xip7N','B10803W17-Xip7S')
#    Conditions = {'Devices.Name=': DeviceNames,
#                  'CharTable.IsOK>': (0, ),
#                  'CharTable.AnalyteCon<': (1e-7, )}
#    
#    GroupBase = {}
#    GroupBase['Table'] = CharTable
#    GroupBase['Last'] = False
#    GroupBase['Conditions'] = Conditions
#    
#    GenFit = GenXlsFittingReport('../testfb.xls', GroupBase)
#    GenFit.DBCalField = 'CharTable.FuncStep'
#    GenFit.DBCalFlags = ('Tromb', )
#    GenFit.GenFullReport()
#    GenFit.close()

#    plt.close('all')   
#    
#    CharTable = 'DCcharacts'
#    DeviceNames = ('B10179W15-T1',)
#    Conditions = {'Devices.Name=': DeviceNames,
#                  'CharTable.IsOK>': (0, ),
#                  'CharTable.Comments like':('%Cal%', )}
#    
#    GroupBase = {}
#    GroupBase['Table'] = CharTable
#    GroupBase['Last'] = False
#    GroupBase['Conditions'] = Conditions
#    
#    GenFit = GenXlsFittingReport('../testf.xls', GroupBase)
#    GenFit.CalMaps = (('PhCal 1', 'IonCal 1'), ('PhCal 2', 'IonCal 2'))
#    GenFit.CalMapVars = {'XVar': 'Ph',
#                  'XVarLog': False,
#                  'YVar': 'IonStrength',
#                  'YVarLog': True,
#                  'ZVar': 'Ud0'}
#    GenFit.GenFullReport()
#    GenFit.close()

#Caltypes={'pHCal':('Ph', False),
#          'IonCal':('IonStrength', True)}

#CalList = DbSearch.FindCommonValues(Conditions=Conditions,
#                                    Parameter='CharTable.Comments',
#                                    Table=CharTable)
#
#for cal in CalList:
#    caltype = cal.split(' ')[0]
#    print 'Cal type', caltype
#    print 'Gen excel', cal
#
## Fixed Conditions
#    GroupBase = {}
#    GroupBase['Table'] = CharTable
#    GroupBase['Last'] = False
#    Conditions['DCcharacts.Comments like']=(cal,)
#    GroupBase['Conditions'] = Conditions.copy()
#    
#    XlsRep=GenXlsFittingReport(cal+'.xls', GroupBase)
#    XlsRep.XVar = Caltypes[caltype][0]
#    XlsRep.XVarLog = Caltypes[caltype][1]
#    XlsRep.GenFullReport()
#    XlsRep.close()



#    Name = 'B10179W15'
#    Conditions = {'Wafers.name=':(Name, )}
#    XlsHist = GenXlsReport('../' +Name+'.xlsx', Conditions)
#    XlsHist.GenFullReport()
#    XlsHist.close()

    
#    Name = 'B10179W15-B4'
#    Conditions = {'Devices.name=':(Name, )}
#    XlsHist = GenXlsTrtsHistory('../' +Name+'.xlsx', Conditions)
#    XlsHist.GenFullReport()
#    XlsHist.close()