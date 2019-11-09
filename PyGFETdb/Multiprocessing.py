import sys

import PyGFETdb.DBAnalyze as DbAn
import PyGFETdb.DBSearch as DbSe
from PyGFETdb import multithrds, Thread


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
    return Results


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
                Data, Trts = processGetFromDB(Thread.call(DbSe, 'GetFromDB', **Grc))
                if Data is not None:
                    ResultsDB[Wfn][Grn] = dict(Data)
        else:
            Data, Trts = processGetFromDB(Thread.call(DbSe, 'GetFromDB', **Wfc))
            if Data is not None:
                ResultsDB[Wfn] = dict(Data)

    return dict(ResultsDB)


def processGetFromDB(results):
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


def processResults(Results, args):
    Ret = {}
    if args is not None:
        for karg, arg in args.items():
            Ret[karg] = {}
            for rd in Results:
                if type(rd) is dict:
                    tdict = rd.get(karg)
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


def SearchDB_MP(GrWfs, **kwargs):
    """

    :param Group: A group of conditions
    :return: The results of the search in the database
    """
    ResultsDB = {}
    thread = Thread.MultiProcess(DbSe)
    for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
        if type(Wfc) is dict and Wfc.get('Conditions') is None:
            for iGr, (Grn, Grc) in enumerate(Wfc.items()):
                key = Wfn + ' ' + Grn
                thread.initcall(key, DbSe)
                thread.call(key, DbSe, 'GetFromDB', **Wfc)
        else:
            thread.initcall(Wfn, DbSe)
            thread.call(Wfn, DbSe, 'GetFromDB', **Wfc)

    for iWf, (Wfn, Wfc) in enumerate(sorted(GrWfs.items())):
        ResultsDB[Wfn] = {}
        if type(Wfc) is dict and Wfc.get('Conditions') is None:
            for iGr, (Grn, Grc) in enumerate(Wfc.items()):
                ResultsDB[Grn] = {}
                key = Wfn + ' ' + Grn
                ResultsDB[Wfn][Grn] = processGetFromDB(thread.getResults(key))[0]
        else:
            ResultsDB[Wfn] = processGetFromDB(thread.getResults(Wfn))[0]

    return ResultsDB
