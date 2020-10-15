#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:12:37 2018

@author: aguimera
"""

from PyGFETdb.DataClass import DataCharAC
import PyGFETdb.DBCore as PyFETdb
import numpy as np
import logging


def GenGroups(GroupBase, GroupBy, LongName=True):
    GroupList = FindCommonValues(Table=GroupBase['Table'],
                                 Conditions=GroupBase['Conditions'],
                                 Parameter=GroupBy)

    Groups = {}
    for Item in GroupList:
        if Item is None:
            continue
        Cgr = GroupBase.copy()
        Cond = GroupBase['Conditions'].copy()
        if Item is None:
            Cond.update({'{} is '.format(GroupBy): (Item,)})
        else:
            Cond.update({'{}='.format(GroupBy): (Item,)})
        Cgr['Conditions'] = Cond
        if LongName:
            GroupName = '{}-{}'.format(GroupBy.split('.')[-1], Item)
        else:
            GroupName = Item
        Groups[GroupName] = Cgr

    return Groups


def GenBiosensGroups(CondBase,
                     GroupBy='CharTable.FuncStep',
                     AnalyteStep='Tromb',
                     AnalyteGroupBy='CharTable.AnalyteCon'):

    Cond = CondBase.copy()
    Conditions = Cond['Conditions'].copy()
    FuncStepList = FindCommonValues(Table=Cond['Table'],
                                    Parameter=GroupBy,
                                    Conditions=Conditions)
    Conditions = Cond['Conditions'].copy()
    Conditions['{}='.format(GroupBy)] = (AnalyteStep,)
    AnalyteConList = FindCommonValues(Table=Cond['Table'],
                                      Parameter=AnalyteGroupBy,
                                      Conditions=Conditions)

    Groups = {}
    for FuncStep in FuncStepList:
        if FuncStep == AnalyteStep:
            for AnalyteCon in AnalyteConList:
                Cgr = CondBase.copy()

                Cond = CondBase['Conditions'].copy()
                Cgr['Conditions'] = Cond
                Cond.update({'{}='.format(AnalyteGroupBy): (AnalyteCon, )})
                Groups['{} {}'.format(FuncStep, AnalyteCon)] = Cgr
        else:
            Cgr = CondBase.copy()

            Cond = CondBase['Conditions'].copy()
            Cgr['Conditions'] = Cond
            Cond.update({'{}='.format(GroupBy): (FuncStep,)})
            Groups[FuncStep] = Cgr

    return Groups


def CheckConditionsCharTable(Conditions, Table):
    for k in list(Conditions.keys()):
        if k.startswith('CharTable'):
            nk = k.replace('CharTable', Table)
            Conditions.update({nk: Conditions[k]})
            del(Conditions[k])
    return Conditions


def FindCommonValues(Parameter, Conditions, Table='ACcharacts', **kwargs):
    Conditions = CheckConditionsCharTable(Conditions, Table)

    if Parameter.startswith('CharTable'):
        Parameter = Parameter.replace('CharTable', Table)

    MyDb = PyFETdb.PyFETdb()
#    MyDb = PyFETdb.PyFETdb()

    Output = (Parameter,)
    Res = MyDb.GetCharactInfo(Table=Table,
                              Conditions=Conditions,
                              Output=Output)

    del (MyDb)
    #  Generate a list of tupples with devices Names and comments
    Values = []
    for Re in Res:
        Values.append(Re[Parameter])

    return set(Values)


def GetFromDB(Conditions, Table='ACcharacts', Last=True, GetGate=True,
              OutilerFilter=None, DataSelectionConfig=None):
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
    OutilerFilter : dictionary, optional. (default 'None'),
        If defined, dictionary to perform a statistical pre-evaluation of the
        data. The data that are not between the p25 and p75 percentile are
        not returned. The dictionary should follow this structure:
        {'Param':Value, --> Characterization parameter, ie. 'Ids', 'Vrms'...
         'Vgs':Value,   --> Vgs evaluation point
         'Vds':Value,   --> Vds evaluationd point
         'Ud0Norm':Boolean} --> Indicates if Vgs is normalized to CNP

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
    Total = float(len(Trts))

    Data = {}
    for Trtn, Cys in DataD.items():
        Chars = []
        for Cyn, Dat in sorted(Cys.items()):
            try:
                Char = DataCharAC(Dat)
                Chars.append(Char)
            except:
                print('Failed loading', Trtn ,' \n')
                print(Dat)                
        Data[Trtn] = Chars

    logging.debug('Getting Data from %s', Conditions)
    print ('Trts Found ->', len(Trts))

    if OutilerFilter is not None:
        logging.debug('Look for Outliers %s', OutilerFilter)
        Data = RemoveOutilers(Data, OutilerFilter)
        Trts = Data.keys()
        logging.debug('Input Trts %d Output Trts d', Total, len(Trts))
        print ('Outlier filter Yield -> ', len(Trts)/Total)

    if DataSelectionConfig is not None:
        Trts = {}
        Trts['Total'] = Data.keys()
        for DataSel in DataSelectionConfig:
            logging.debug('Look for data in range %s', DataSel)
            if 'Name' not in DataSel:
                DataSel['Name'] = DataSel['Param']
            Data = DataSelection(Data, **DataSel)
            Trts[DataSel['Name']] = Data.keys()
            logging.debug('Input Trts %d Output Trts %d', Total, len(Trts))

        Trts['Final'] = Data.keys()
        for DataSel in DataSelectionConfig:
            name = DataSel['Name']
            v = float(len(Trts[name]))
            if Total>0:
                print (name, ' Yield -> ', v/Total)

    return Data, Trts


def RemoveOutilers(Data, OutilerFilter):
    Vals = np.array([])
    for Trtn, Datas in Data.items():
        for Dat in Datas:
            func = Dat.__getattribute__('Get' + OutilerFilter['Param'])
            Val = func(Vgs=OutilerFilter['Vgs'],
                       Vds=OutilerFilter['Vds'],
                       Ud0Norm=OutilerFilter['Ud0Norm'])
            if Val is not None:
                Vals = np.hstack((Vals, Val)) if Vals.size else Val

    p25 = np.percentile(Vals, 25)
    p75 = np.percentile(Vals, 75)
    lower = p25 - 1.5 * (p75 - p25)
    upper = p75 + 1.5 * (p75 - p25)

    DataFilt = {}
    for Trtn, Datas in Data.items():
        Chars = []
        for Cyn, Char in enumerate(Datas):
            func = Char.__getattribute__('Get' + OutilerFilter['Param'])
            Val = func(Vgs=OutilerFilter['Vgs'],
                       Vds=OutilerFilter['Vds'],
                       Ud0Norm=OutilerFilter['Ud0Norm'])

            if (Val <= lower or Val >= upper):
                print ('Outlier Removed ->', Val, Trtn, Cyn)
            else:
                Chars.append(Char)
        DataFilt[Trtn] = Chars

    return DataFilt


def DataSelection(Data, Param, Range, Function=None, InSide=True, Name=None,
                  ParArgs={'Vgs': None,
                           'Vds': None,
                           'Ud0Norm': False}):

    DataFilt = {}
    for Trtn, Datas in Data.items():
        DatTrt = []
        for Dat in Datas:
            func = Dat.__getattribute__('Get' + Param)
            Val = func(**ParArgs)
            if Function is not None:
                Val = Function(Val)

            if InSide:
                if Range[0] is None:
                    MinCond = False
                else:
                    MinCond = (Val < Range[0]).any()

                if Range[1] is None:
                    MaxCond = False
                else:
                    MaxCond = (Val > Range[1]).any()
                FinalCond = MinCond | MaxCond
            else:
                if Range[0] is None:
                    MinCond = True
                else:
                    MinCond = (Val > Range[0]).any()

                if Range[1] is None:
                    MaxCond = True
                else:
                    MaxCond = (Val < Range[1]).any()
                FinalCond = MinCond & MaxCond

            if FinalCond:
                logging.debug('Meas Out %s %s %f' , Trtn, Dat.GetTime(), Val)
                continue
            DatTrt.append(Dat)
        if len(DatTrt) > 0:
            DataFilt[Trtn] = DatTrt
    return DataFilt



def UpdateCharTableField(Conditions, Value,
                         Table='ACcharacts', Field='Comments'):

    Conditions = CheckConditionsCharTable(Conditions, Table)

    MyDb = PyFETdb.PyFETdb()

    out = '{}.id{}'.format(Table, Table)
    re = MyDb.GetCharactInfo(Table=Table,
                             Conditions=Conditions,
                             Output=(out, ))

    print (re)
    text = "Do you wan to update {} in {} for {} y/n ?".format(Field, Table, Value)
    inText = input(text)
    if inText =='y':
        print ('Updated')
        field = '{}.{}'.format(Table, Field)
        fields = {field: Value}
        for r in re:
            condition = ('{}.id{}='.format(Table, Table), list(r.values())[0])
            MyDb.UpdateRow(Table=Table, Fields=fields, Condition=condition)
    else:
        print ('Cancelled')

    MyDb.db.commit()
