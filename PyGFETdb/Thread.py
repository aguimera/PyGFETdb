# import threading as th
from multiprocessing import pool, Lock

# import numpy as np
import PyGFETdb

rets = {}
lock = Lock()


class PyFETdb:
    def __init__(self, package):  # , host=DBhost, user=DBuser, passwd=DBpasswd, db=DBdb, Update=True):
        self.pool = pool.Pool()
        self.parent = package
        rets = []

    def call(self, funcname, args, **kwargs):
        func = self.parent.__getattribute__(funcname)
        ret = self.pool.apply_async(func, args, kwargs, error_callback=errorlog, callback=addResult)
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
    rets[result[0]] = result[2]
    PyGFETdb.Thread.lock.release()


def errorlog(e):
    print(e)
