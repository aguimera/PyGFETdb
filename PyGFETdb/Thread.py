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
        self.pool = pool.ThreadPool()
        self.parent = package
        rets = []

    def call(self, funcname, arguments):
        ret = None
        if self.parent is not None:
            func = getattr(self.parent, funcname)
            if not callable(func):
                func = self.parent.__getattribute__(self.parent, funcname)

            # elif hasattr(self.parent,funcname) and inspect.ismethod(getattr(self.parent,funcname)):
            #   inspect.getargs()

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

    def __del__(self):
        self.pool.close()
        self.pool.join()
        rets = {}


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
        return tdict
    else:
        func = klass.__getattribute__(function)
        res = func(**kwargs)
    return res
