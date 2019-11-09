# import threading as th
from multiprocessing import pool, Lock

# import numpy as np
from PyGFETdb import multithrds


class Thread(pool.ThreadPool):
    def __init__(self, package):  # , host=DBhost, user=DBuser, passwd=DBpasswd, db=DBdb, Update=True):
        pool.ThreadPool.__init__(self)
        self.pool = self
        self.parent = package
        self._NUMTHREADS = 10000000
        self._rets = []
        self.lock = Lock()
        self.args = None

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
                self.args = i
                for karg, varg in i.items():
                    args = []
                    for kargument, vargument in arguments.items():
                        if kargument != "args":
                            args.append(vargument)
                    args.append(dict({karg: varg}))
                    ret = self.pool.apply_async(func, args, error_callback=self.errorlog, callback=self.addResult)
            else:
                ret = self.pool.apply_async(func, [], arguments, error_callback=self.errorlog, callback=self.addResult)
        return ret

    def __del__(self):
        self.pool.close()
        self.pool.join()
        self.rets = {}


    def getResults(self):
        self.pool.close()
        self.pool.join()
        self.pool.terminate()
        return self._rets

    def addResult(self, result):
        self.lock.acquire()
        self._rets.append(result)
        self.lock.release()

    def errorlog(self, e):
        print(e)


##########################################################

lock = Lock()


class MultiProcess():
    def __init__(self, klass):
        self.pool = {}

    def initcall(self, key, klass):
        self.pool[key] = Thread(klass)

    def call(self, key, klass, function, **kwargs):
        res = None
        if multithrds:  # is not None:
            self.pool[key].call(function, kwargs)
        else:
            func = getattr(klass, function)
            if not callable(func):
                func = klass.__getattribute__(function)
                res = klass.func(**kwargs)
            else:
                res = func(**kwargs)
        return res

    def getResults(self, key):
        ret = {}
        res = self.pool[key].getResults()
        for item in res:
            ret.update({key: item})
        del self.pool[key]
        return ret


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