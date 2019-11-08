# import threading as th
import random
from multiprocessing import pool, Lock

# import numpy as np
import PyGFETdb
from PyGFETdb import multithrds

NUMTHREADS = 10000
rets = {}
lock = Lock()


class Thread:
    def __init__(self, package):  # , host=DBhost, user=DBuser, passwd=DBpasswd, db=DBdb, Update=True):
        self.pool = pool.Pool()
        self.parent = package
        rets = []

    def call(self, funcname, arguments):
        func = self.parent.__getattribute__(funcname)
        i = arguments.get('args')
        if i is not None:
            for karg, varg in i.items():
                args = []
                for kargument, vargument in arguments.items():
                    if kargument != "args":
                        args.append(vargument)
                args.append(dict({karg: varg}))
                ret = self.pool.apply_async(func, args, error_callback=errorlog, callback=addResult)
        else:
            args = []
            for kargument, vargument in arguments.items():
                args.append(vargument)
            ret = self.pool.apply_async(func, [], arguments, error_callback=errorlog, callback=addResult)
        return ret

    # def getPool(self):
    #    return self.pool

    def __del__(self):
        self.pool.close()
        self.pool.join()
        rets = []

    def getResults(self):
        self.pool.close()
        self.pool.join()
        self.pool.terminate()
        return rets


def addResult(result):
    PyGFETdb.Thread.lock.acquire()
    randomkey = random.randint(0, NUMTHREADS)
    rets[randomkey] = result
    PyGFETdb.Thread.lock.release()


def errorlog(e):
    print(e)


def call(klass, function, **kwargs):
    if multithrds:  # is not None:
        pool = Thread(klass)
        pool.call(function, kwargs)
        tdict = pool.getResults()
        del pool
        return processResults(tdict, kwargs.get('args'))
    else:
        func = klass.__getattribute__(function)
        res = func(**kwargs)
    return res


def processResults(ResultsDict, args):
    Results = {}
    for karg, arg in args.items():
        Results[karg] = {}
        for r, rd in ResultsDict.items():
            tdict = rd.get(karg)
            if tdict is not None:
                for iWf, (Wfn, Wfc) in enumerate(tdict.items()):
                    Results[karg][Wfn] = {}
                    if type(Wfc) is dict and Wfc.get('Conditions') is None:
                        for iGr, (Grn, Grc) in enumerate(tdict.items()):
                            Results[karg][Wfn][Grn] = Grc
                    else:
                        Results[karg][Wfn] = Wfc
    return Results
