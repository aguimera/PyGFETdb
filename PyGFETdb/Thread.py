# import threading as th
import random
from multiprocessing import pool, Lock

# import numpy as np
import PyGFETdb

NUMTHREADS = 10000
rets = {}
lock = Lock()


class PyFETdb:
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
