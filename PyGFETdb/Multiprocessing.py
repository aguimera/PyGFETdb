# -*- coding: utf-8 -*-
"""

@author: dragc

"""

import itertools
import sys

import numpy as np
import quantities as pq

import PyGFETdb
import PyGFETdb.DBAnalyze as DbAn
import PyGFETdb.DBCore as PyFETdb
import PyGFETdb.DBSearch as DbSe
import PyGFETdb.DataClass as pData
from PyGFETdb import qty, multithrds, superthreading, Thread
from PyGFETdb.DataClass import DataCharAC

getParam = 'GetParam'
getFromDB = 'GetFromDB'

if not superthreading:
    getparamclass = DbAn
    getclass = DbSe
else:
    getparamclass = ''
    getclass = ''


########################################################################
#
#  NORMAL FUNCTIONS
#
########################################################################
def SearchDB(GrWfs, **kwargs):
    """

    :param Group: A group of conditions
    :return: The results of the search in the database
    """
    ResultsDB = {}
    for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
        ResultsDB[Wfn] = {}
        if type(Wfc) is dict and Wfc.get('Conditions') is None:
            for iGr, (Grn, Grc) in enumerate(Wfc.items()):
                kwargs.update(**Grc)
                Data, Trts = DbSe.GetFromDB(**kwargs)
                if Data is not None:
                    ResultsDB[Wfn][Grn] = dict(Data)
        else:
            kwargs.update(**Wfc)
            Data, Trts = DbSe.GetFromDB(**kwargs)
            if Data is not None:
                ResultsDB[Wfn] = dict(Data)

    return dict(ResultsDB)


def GetParams(ResultsDB, GrWfs, args: dict, **kwargs):
    """

    :param ResultsDB: The results of a search in the database
    :param Group: A group of conditions to analyse
    :param args: Arguments for getting the parameters
    :param Plot: Bool that activates plotting
    :param kwargs:
    :return: A dict of args of dicts of groupnames and parameter found in a previous search
    """
    Results = {}
    for karg, arg in args.items():
        Results[karg] = {}
        for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
            Results[karg][Wfn] = {}
            if type(Wfc) is dict and Wfc.get('Conditions') is None:
                for iGr, (Grn, Grc) in enumerate(sorted(Wfc.items())):
                    Data = ResultsDB.get(Grn)
                    if Data is not None:
                        ParamData = DbAn.GetParam(Data, **arg)
                        if ParamData is not None:
                            Results[karg][Wfn][Grn] = ParamData
            else:
                Data = ResultsDB.get(Wfn)
                if Data is not None:
                    ParamData = DbAn.GetParam(Data, **arg)
                    if ParamData is not None:
                        Results[karg][Wfn] = ParamData
    return processResults(Results, args)


def processResults(Results, args):
    """

    :param Results: results to classify
    :param args: arguments to search
    :return: A dict of arguments with each corresponding result
    """
    Ret = {}
    if args is not None:
        for karg, arg in args.items():
            Ret[karg] = {}
            #            for rk,rd in Results.items():
            if type(Results) is dict:
                tdict = Results.get(karg)
                if tdict is not None:
                        for iWf, (Wfn, Wfc) in enumerate(tdict.items()):
                            Ret[karg][Wfn] = {}
                            if type(Wfc) is dict and Wfc.get('Conditions') is None:
                                for iGr, (Grn, Grc) in enumerate(tdict.items()):
                                    Ret[karg][Wfn][Grn] = Grc
                            else:
                                Ret[karg][Wfn] = Wfc
    else:
        Ret = Results

    return Ret


########################################################################
#
#  MULTIPROCESSING
#
########################################################################
def SearchDB_MP(GrWfs, **kwargs):
    """

    :param Group: A group of conditions
    :return: The results of the search in the database
    """
    ResultsDB = {}

    if PyGFETdb.Multiprocessing.getclass == '':
        getclass = PyGFETdb.Multiprocessing
    else:
        getclass = PyGFETdb.Multiprocessing.getclass

    p = (len(GrWfs) * len(GrWfs.items()))
    n = int(10 / p + p / 100 / 5000 + 5)
    print()
    print("Searching in DB {} times... forking {} threads".format(p, n))
    thread = Thread.MultiProcess(getclass, n)
    for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
        if type(Wfc) is dict and Wfc.get('Conditions') is None:
            for iGr, (Grn, Grc) in enumerate(Wfc.items()):
                key = Wfn + ' ' + Grn
                thread.initcall(key, getclass)
                kwargs.update(**Grc)
                thread.call(key, getclass, getFromDB, {}, **kwargs)
        else:
            thread.initcall(Wfn, DbSe)
            kwargs.update(**Wfc)
            thread.call(Wfn, getclass, getFromDB, {}, **kwargs)

    for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
        ResultsDB[Wfn] = {}
        if type(Wfc) is dict and Wfc.get('Conditions') is None:
            for iGr, (Grn, Grc) in enumerate(Wfc.items()):
                ResultsDB[Grn] = {}
                key = str(Wfn + ' ' + Grn)
                ret = thread.getResults(key)[key]
                ResultsDB[Wfn][Grn] = processGetFromDB(ret)[0]
        else:
            ret = thread.getResults(Wfn)[Wfn]
            ResultsDB[Wfn] = processGetFromDB(ret)[0]
    thread.end()
    del thread
    print('DB Search... Done')
    return ResultsDB


def processGetFromDB(results):
    """
        **Gathers the results of different processes after a search in the database**

    :param results: Results obtained with multi-processing from the database
    :return: All the results from the search gathered in a tuple
    """
    ret = (None, None)
    if not multithrds:
        ret = results
    else:
        if type(results) is dict:
            Data = {}
            Trts = []
            try:
                for i, (data, trts) in results.items():
                    Data.update(data)
                    Trts.append(trts)
                ret = (Data, Trts)
            except ValueError:
                print(sys.exc_info())
        elif type(results) is list:
            Data = {}
            Trts = []
            try:
                for data, trts in results:
                    Data.update(data)
                    Trts.append(trts)
                ret = (Data, Trts)
            except ValueError:
                print(sys.exc_info())

    return ret


def GetParams_MP(ResultsDB, GrWfs, arguments, **kwargs):
    """

    :param ResultsDB: The results of a search in the database
    :param Group: A group of conditions to analyse
    :param args: Arguments for getting the parameters
    :param Plot: Bool that activates plotting
    :param kwargs:
    :return: A dict of args of dicts of groupnames and parameter found in a previous search
    """

    if PyGFETdb.Multiprocessing.getparamclass == '':
        getparamclass = PyGFETdb.Multiprocessing
    else:
        getparamclass = PyGFETdb.Multiprocessing.getparamclass

    p = (len(arguments) * len(GrWfs) * len(GrWfs.items()))
    n = int(100 / p + p / 5000 / 50000 + 5)
    print("Searching {} Parameters... forking {} threads".format(p, n))
    thread = Thread.MultiProcess(getparamclass, n)
    for karg, arg in arguments.items():
        print('Searching Parameter {}'.format(arg.get('Param')))
        for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
            if type(Wfc) is dict and Wfc.get('Conditions') is None:
                for iGr, (Grn, Grc) in enumerate(sorted(Wfc.items())):
                    # print('Searching Parameter {} for {}/{}'.format(arg.get('Param'), Wfn,Grn))
                    Data = ResultsDB.get(Grn)
                    if Data is not None:
                        key = str(karg + ' ' + Wfn + ' ' + Grn)
                        kargs = {'Data': Data, 'args': arg}
                        thread.initcall(key, getparamclass)
                        thread.call(key, getparamclass, getParam, kargs, **kargs)

            else:
                # print('Searching Parameter {} for {}'.format(arg.get('Param'),Wfn))
                Data = ResultsDB.get(Wfn)
                if Data is not None:
                    key = str(karg + ' ' + Wfn)
                    thread.initcall(key, getparamclass)
                    kargs = {'Data': Data, 'args': arg}
                    thread.call(key, getparamclass, getParam, kargs, **kargs)

    print('Retrieving Parameters...')
    # """""
    res = {}
    for karg, arg in arguments.items():
        for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
            if type(Wfc) is dict and Wfc.get('Conditions') is None:
                for iGr, (Grn, Grc) in enumerate(sorted(Wfc.items())):
                    print('Searching Parameter {} for {}/{}'.format(arg.get('Param'), Wfn, Grn))
                    key = str(karg + ' ' + Wfn + ' ' + Grn)
                    res[key] = thread.getResults(key)
            else:
                print('Retrieving Parameter {} for {}'.format(arg.get('Param'), Wfn))
                key = str(karg + ' ' + Wfn)
                res[key] = thread.getResults(key)
    thread.end()
    del thread
    print('Searching Parameters Done.')
    Results = processGetParams_MP(GrWfs, res, arguments)
    return Results


def processGetParams_MP(GrWfs, Results, args):
    """

    :param GrWfs: A group of conditions
    :param Results: The results of the parameter search
    :param args: a dict with the parameters to be searched
    :return: The results ordered by parameter
    """
    Ret = {}
    if args is not None:
        for karg, arg in args.items():
            Ret[karg] = {}
            if type(Results) is dict:
                processOneArg(GrWfs, Results, Ret, karg)
    else:
        Ret = Results

    return Ret


def processOneArg(GrWfs, Results, Ret, karg):
    """

    :param GrWfs: A group of conditions
    :param Results: Results of the parameter search
    :param Ret: A dict to update with Eesults
    :param karg: Key-name of the argument to process
    :return: None but modifies Ret dict
    """
    for iWf, (Wfn, Wfc) in enumerate(GrWfs.items()):
        Ret[karg][Wfn] = {}
        if type(Wfc) is dict and Wfc.get('Conditions') is None:
            for iGr, (Grn, Grc) in enumerate(Wfc.items()):
                for key, result in Results.items():
                    if key.startswith(karg + ' ' + Wfn + ' ' + Grn):
                        k = result.get(key)
                        if k is not None:
                            l = list(itertools.chain(k))
                            Ret[karg][Wfn][Grn] = l[0][0]
        else:
            for key, result in Results.items():
                if key.startswith(karg + ' ' + Wfn):
                    k = result.get(key)
                    if k is not None:
                        l = list(itertools.chain(k))
                        Ret[karg][Wfn] = l[0][0]


def GetParam(Data, Param, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
    """
    **Multiprocessing Version of GetParam**

    :param Data:
    :param Param:
    :param Vgs:
    :param Vds:
    :param Ud0Norm:
    :param kwargs:
    :return:
    """
    Vals = qty.createQuantityList()

    if Data is None:
        return Vals

    ret = []
    # print('Getting Parameter {}...'.format(Param))
    thread = Thread.MultiProcess(pData.DataCharAC)
    kwargs.update({'Param': Param, 'Vds': Vds, 'Vgs': Vgs, 'Ud0Norm': Ud0Norm})
    results = {}
    for Trtn, Datas in Data.items():
        results[Trtn] = {}
        key = thread.initcall(Thread.key(), pData.DataCharAC)
        for Dat in Datas:
            args = {'args': {'self': Dat}}
            args.update(kwargs)
            thread.call(key, pData.DataCharAC, 'Get' + Param, args, **args)
        results[Trtn] = thread.getResults(key)[key][0]
    for Trtn, Val in results.items():
        if Val is not None:
            if type(Val) is pq.Quantity \
                    or Param == 'PSD' or Param == 'Fpsd' or Param == 'NoA' or Param == 'NoB':
                Vals = qty.appendQuantity(Vals, Val)
            else:
                Vals = np.array(Vals)
                try:
                    Vals = np.hstack(((Vals), Val)) if Vals.size else Val
                except ValueError:
                    # print(sys.exc_info())
                    ret.append(Vals)
                    Vals = qty.createQuantityList()
                    Vals = qty.appendQuantity(Vals, Val)
                    # raise ArithmeticError # FIXME:
    thread.end()
    del thread

    # print("Getting Parameter {} Done.".format(Param))
    if len(ret) > 1:
        return ret, results
    else:
        return Vals, results


def GetFromDB(Conditions, Table='ACcharacts', Last=True, GetGate=True,
              OutilerFilter=None, DataSelectionConfig=None, remove50Hz=False):
    """

        **Get data from data base**

        This function returns data which meets with "Conditions" dictionary for sql
        select query constructor.

        :param Conditions: dictionary, conditions to construct the sql select query.

        The dictionary should follow this structure:

            {'Table.Field <sql operator>' : iterable type of values}

            - Example:

                {'Wafers.Name = ':(B10803W17, B10803W11),'CharTable.IsOK > ':(0,)}

        :param Table: string, optional.

        Posible values 'ACcharacts' or 'DCcharacts'.

        The default value is 'ACcharacts'. Characterization table to get data

        The characterization table of Conditions dictionary can be indicated
        as 'CharTable'. In that case 'CharTable' will be replaced by Table value.

        :param Last: boolean, optional.

        If True (default value) just the last measured data for each transistor is returned.

        If False, all measured data is returned

        :param Last: boolean, optional.

        If True (default value) the gate measured data is also obtained

        :param OutilerFilter: dictionary, optional.  (default 'None'),

        If defined, dictionary to perform a statistical pre-evaluation of the
        data.

        The data that are not between the p25 and p75 percentile are
        not returned.

            The dictionary should follow this structure:

                {'Param':Value, --> Characterization parameter, ie. 'Ids', 'Vrms'...

                'Vgs':Value,   --> Vgs evaluation point

                'Vds':Value,   --> Vds evaluationd point

                'Ud0Norm':Boolean} --> Indicates if Vgs is normalized to CNP

        :param remove50Hz: bool to activate the removal of frequency 50Hz


        :return: A Dictionary with the data arranged as follows:

            {'Transistor Name':list of PyGFET.DataClass.DataCharAC classes}

        :return: A List of transistors


    """
    selfclass = PyGFETdb.Multiprocessing

    # logging.basicConfig(filename=log, level=logging.DEBUG)

    Conditions = DbSe.CheckConditionsCharTable(Conditions, Table)

    MyDb = PyFETdb.PyFETdb()

    DataD, Trts = MyDb.GetData2(Conditions=Conditions,
                                Table=Table,
                                Last=Last,
                                GetGate=GetGate,
                                remove50Hz=remove50Hz)

    del (MyDb)
    Trts = list(Trts)
    Total = float(len(Trts))

    Data = {}
    for Trtn, Cys in DataD.items():
        Chars = []
        for Cyn, Dat in sorted(Cys.items()):
            Char = DataCharAC(Dat)
            Chars.append(Char)
        Data[Trtn] = Chars

    #    logging.debug('Getting Data from %s', Conditions)
    print('Trts Found ->', len(Trts))

    if OutilerFilter is not None:
        #       logging.debug('Look for Outliers %s', OutilerFilter)
        Data = DbSe.RemoveOutilers(Data, OutilerFilter)
        Trts = Data.keys()
        #      logging.debug('Input Trts %d Output Trts d', Total, len(Trts))
        print('Outlier filter Yield -> ', qty.Divide(len(Trts), Total))

    if DataSelectionConfig is not None:
        thread = Thread.MultiProcess(selfclass)
        key = thread.initcall(Thread.key(), selfclass)
        Trts = {}
        Trts['Total'] = Data.keys()
        # print('DataSelection...')
        for DataSel in DataSelectionConfig:
            #         logging.debug('Look for data in range %s', DataSel)
            if 'Name' not in DataSel:
                DataSel['Name'] = DataSel['Param']
            args = {'args': {'Data': Data}}
            args.update(DataSel)
            thread.call(key, selfclass, 'DataSelection', args, **args)

        DataFilt = {}
        ret = thread.getResults(key)[key]  # _DataSelection(Data, **DataSel)
        for item in ret:
            DataFilt.update(item)
        Trts = DataFilt.keys()
        Data = DataFilt
        #print('DataSelection Done.')
        thread.end()
        del thread
    return Data, Trts


def DataSelection(Data, Param, Range, Function=None, InSide=True, Name=None, Units=None,
                  ParArgs={'Vgs': None,
                            'Vds': None,
                            'Ud0Norm': False,
                            'Units': None}):
    rPAr, tPar = GetParam(Data, Param, **ParArgs)
    DataFilt = {}
    for Trtn, Datas in tPar.items():
        if Datas is not None and Datas.size > 0:
            if type(Datas) is list:
                for Val in Datas:
                    if checkRanges(Val, Units, Range, InSide):
                        DataFilt.update({Trtn: Data[Trtn]})
            else:
                if checkRanges(Datas, Units, Range, InSide):
                    DataFilt.update({Trtn: Data[Trtn]})
    return DataFilt


def checkRanges(Val, Units, Range, InSide):
    # Added extra checks to make scripts compatible
    # begin fix
    if Val is None:
        return
    if type(Val) is not tuple:
        Val = (Val)
    # Added Quantity Support for the Range
    if qty.isActive() and Units is not None:
        try:
            Range = pq.Quantity(Range, Units)
        except TypeError as e:
            print("Range Units Error: ", sys.exc_info())
    # end fix
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
        # logging.debug('Meas Out %s %s %f', Trtn, Dat.GetTime(), Val)
        return
    return True
