#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:12:37 2018

@author: aguimera
"""

from PyGFETdb.DataClass import DataCharAC, DataCharDC
import PyGFETdb.DBCore as PyFETdb
import numpy as np
import pandas as pd
import quantities as pq
import pickle

# Array values on measured Vgs range
OutUnits = {'Ids': 'uA',
            'GM': 'mS',
            'GMV': 'mS/V',
            'Irms': 'uA',
            'Vrms': 'uV',
            'NoA': 'A**2',
            'NoC': 'A**2',
            }

########################################################################################################################
## Define default scalar parameters
########################################################################################################################
ScalarParamsDC = ['Ids', 'GM', 'GMV']
ScalarParamsAC = ['Irms', 'Vrms', 'NoA', 'NoB', 'NoC']
VgsScalar = -0.1*pq.V
ScalarQueriesDC = {'CNP': {'Param': 'Ud0',
                            'Units': 'mV'},
                   'IgMax': {'Param': 'IgMax',
                              'Units': 'nA'
                              },
                   'IdsV01': {'Param': 'Ids',
                              'Vgs': 0.1*pq.V,
                              'Ud0Norm': False,
                              'Units': 'uA'
                              },
                   'RdsCNP': {'Param': 'Rds',
                              'Vgs': 0*pq.V,
                              'Ud0Norm': True,
                              'Units': 'kOhm'
                              },
                   'IdsCNP': {'Param': 'Ids',
                              'Vgs': 0*pq.V,
                              'Ud0Norm': True,
                              'Units': 'uA'
                              } ,
                   }

for par in ScalarParamsDC:
    d = {'Param': par,
         'Vgs': VgsScalar,
         'Ud0Norm': True}
    if par in OutUnits:
        d['Units'] = OutUnits[par]
    ScalarQueriesDC[par+'01'] = d

ScalarQueries = ScalarQueriesDC.copy()
for par in ScalarParamsAC:
    d = {'Param': par,
         'Vgs': VgsScalar,
         'Ud0Norm': True}
    if par in OutUnits:
        d['Units'] = OutUnits[par]
    ScalarQueries[par+'01'] = d

########################################################################################################################
## Define default Array parameters
########################################################################################################################
Vgs = np.linspace(-0.1, 0.6, 100) * pq.V
VgsNorm = np.linspace(-0.4, 0.4, 100) * pq.V
ArrayParamsDC = ['Ids', 'GM', 'GMV', 'Ig']
ArrayParamsAC = ['Irms', 'Vrms', 'NoA', 'NoB', 'NoC']
ArrayQueriesDC = {}
for par in ArrayParamsDC:
    d = {'Param': par,
         'Vgs': Vgs,
         'Ud0Norm': False}
    if par in OutUnits:
        d['Units'] = OutUnits[par]
    ArrayQueriesDC[par] = d

    d = {'Param': par,
         'Vgs': VgsNorm,
         'Ud0Norm': True}
    if par in OutUnits:
        d['Units'] = OutUnits[par]
    ArrayQueriesDC[par+'Norm']=d

ArrayQueries = ArrayQueriesDC.copy()
for par in ArrayParamsAC:
    d = {'Param': par,
         'Vgs': Vgs,
         'Ud0Norm': False}
    if par in OutUnits:
        d['Units'] = OutUnits[par]
    ArrayQueries[par] = d

    d = {'Param': par,
         'Vgs': VgsNorm,
         'Ud0Norm': True}
    if par in OutUnits:
        d['Units'] = OutUnits[par]
    ArrayQueries[par+'Norm']=d

########################################################################################################################
## Define default Array parameters
########################################################################################################################

ClassQueries = ScalarQueries.copy()
ClassQueries.update(ArrayQueries)
pdAttr = {'Vgs': Vgs,
          'VgsNorm': VgsNorm,
          'ScalarCols': list(ScalarQueries.keys()),
          'ArrayCols': list(ArrayQueries.keys()),
          }

ClassQueriesDC = ScalarQueriesDC.copy()
ClassQueriesDC.update(ArrayQueriesDC)
pdAttrDC = {'Vgs': Vgs,
            'VgsNorm': VgsNorm,
            'ScalarCols': list(ScalarQueriesDC.keys()),
            'ArrayCols': list(ArrayQueriesDC.keys()),
            }


def CalcElectricalParams(dbRaw, ClsQueries, dfAttr, pErrors=False):
    PdSeries = []
    for index, row in dbRaw.iterrows():
        print("Calculating {} of {}".format(index, dbRaw.shape[0]))
        char = row['CharCl']
        Vds = char.GetVds()
        for vds in Vds:
            vals = {}
            vals['Vds'] = vds.flatten()
            for parn, park in ClsQueries.items():
                park['Vds'] = vds
                try:
                    val = char.Get(**park)
                except:
                    if pErrors:
                        print('Error', parn, char.Name)
                    continue
                if val is None:
                    continue
                if val.size > 1:
                    vals[parn] = val
                else:
                    vals[parn] = val.flatten()
            PdSeries.append(pd.concat((row, pd.Series(vals))))
    
    dfDat = pd.concat(PdSeries, axis=1).transpose()       
    
    DataTypes = dbRaw.dtypes
    DataTypes['Vds'] = np.float
    for col in dfAttr['ScalarCols']:
        if col in dfDat.columns:
            DataTypes[col] = np.float
    for col in dfAttr['ArrayCols']:
        if col in dfDat.columns:
            DataTypes[col] = object
    dfDat = dfDat.astype(DataTypes)
    
    ColUnits = {}
    for col, p in ClsQueries.items():
        if 'Units' in p:
            ColUnits[col] = p['Units']
        else:
            ColUnits[col] = ''
    
    dfAttr['ColUnits'] = ColUnits
    
    dfDat.attrs.update(dfAttr)

    return dfDat

def LoadPickleData(FileIn):
    DataIn = pickle.load(open(FileIn, 'rb'), encoding='latin')
    DevDCVals = DataIn['DevDC']
    DevACVals = DataIn['DevAC']

    DictClDC = {}
    if DevACVals is not None:
        DictClAC = {}
    else:
        DictClAC = None

    GateDict = None
    for chn, dat in DevDCVals.items():
        if chn == 'Gate':
            GateDict = dat
            continue
        DCclass = DataCharDC(dat)
        DictClDC[chn] = DCclass

        if DevACVals is None:
            continue

        acDat = DCclass.__dict__.copy()
        acDat.update(DevACVals[chn])
        acDat.pop('Ids')
        acDat['Vgs'] = acDat.pop('VgsAC')
        acDat['Vds'] = acDat.pop('VdsAC')
        DictClAC[chn] = DataCharAC(acDat)
    
    return DictClDC, DictClAC, GateDict


def CheckConditionsCharTable(Conditions, Table):
    for k in list(Conditions.keys()):
        if k.startswith('CharTable'):
            nk = k.replace('CharTable', Table)
            Conditions.update({nk: Conditions[k]})
            del(Conditions[k])
    return Conditions


def GetFromDB(Conditions, Table='ACcharacts', Last=True, GetGate=True):
    """
    Get data from data base

    This function returns data which meets with "Conditions" dictionary for sql
    selct query constructor.

    Parameters
    ----------
    Conditions : dictionary, conditions to construct the sql select query.
        The dictionary should follow this structure:\n
        {'Table.Field <sql operator>' : iterable type of values}
        \nExample:\n
        {'Wafers.Name = ':(B10803W17, B10803W11),
        'CharTable.IsOK > ':(0,)}
    Table : string, optional. Posible values 'ACcharacts' or 'DCcharacts'.
        The default value is 'ACcharacts'. Characterization table to get data
        \n
        The characterization table of Conditions dictionary can be indicated
        as 'CharTable'. In that case 'CharTable' will be replaced by Table
        value. 
    Last : bolean, optional. If True (default value) just the last measured
        data for each transistor is returned. If False, all measured data is
        returned
    Last : bolean, optional. If True (default value) the gate measured data
        is also obtained

    Returns
    -------
    Return : tupple of (Data, Trts)
    Data: Dictionary with the data arranged as follows:
        {'Transistor Name':list of PyGFET.DataClass.DataCharAC classes}

    Trts: List of transistors

    Examples
    --------

    """

    Conditions = CheckConditionsCharTable(Conditions, Table)

    MyDb = PyFETdb.PyFETdb()

    DataD, Trts = MyDb.GetData2(Conditions=Conditions,
                                Table=Table,
                                Last=Last,
                                GetGate=GetGate)

    del(MyDb)
    Trts = list(Trts)

    Data = {}
    for Trtn, Cys in DataD.items():
        Chars = []
        for Cyn, Dat in sorted(Cys.items()):
            try:
                Char = DataCharAC(Dat)
                Chars.append(Char)
            except:
                print('Failed loading', Trtn)
                print(Dat)
        Data[Trtn] = Chars

    print ('Trts Found ->', len(Trts))
    return Data, Trts


def Data2Pandas(DictData):
    pdSeries = []
    for tn, dats in DictData.items():
        for dd in dats:
            pdser = {}
            pdser['Name'] = tn
            pdser['CharCl'] = DataCharAC(dd)
            pdser['IsOk'] = dd.IsOK
            pdser['ChName'] = dd.ChName
            pdser['Date'] = dd.GetDateTime()
            for k, v in dd['Wafers'].items():
                if k == 'Name':
                    pdser['Wafer'] = v
                else:
                    pdser['Waf' + k] = v
            for k, v in dd['Devices'].items():
                if k == 'Name':
                    pdser['Device'] = v
                else:
                    pdser['Dev' + k] = v
            for k, v in dd['TrtTypes'].items():
                if k == 'Name':
                    pdser['TrtType'] = v
                else:
                    pdser[k] = v
            for k, v in dd['Info'].items():
                if k == 'Gate_id':
                    continue
                pdser[k] = v
            pdSeries.append(pd.Series(pdser))
        
    dfRaw = pd.concat(pdSeries, axis=1).transpose()
    DataTypes = {}
    for col in dfRaw.keys():    
        if col in ('Width', 'Length', 'Pass', 'Area', 'Ph', 'IonStrength', 'AnalyteCon'):
            DataTypes[col] = np.float
        else:
            DataTypes[col] = 'category'
    
    DataTypes['CharCl'] = object
    DataTypes['IsOk'] = bool
    DataTypes['Date'] = 'datetime64[ns]'
    
    return dfRaw.astype(DataTypes)


