# -*- coding: utf-8 -*-
"""

@author: dragc

"""
import random
from multiprocessing import pool, Lock

from PyGFETdb import multithrds


class Thread(pool.ThreadPool):
    def __init__(self, package):
        """

        :param package: Module where the function to process is

        """

        pool.ThreadPool.__init__(self)
        self.pool = self
        self.parent = package
        self._NUMTHREADS = 10000000
        self._rets = []
        self.lock = Lock()
        self.args = None

    def call(self, funcname, arguments, **kwargs):
        """
            **Calls a function for the Thread to process**

        :param funcname: Name of the function to call
        :param arguments: Arguments of the call
        :param kwargs:
        :return: None

        """

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
                # """""
                if len(i) > 0:
                    # for karg, varg in i.items():
                    args = {}
                    for kargument, vargument in arguments.items():
                        if kargument != "args":
                            args.update({kargument: vargument})
                    args.update(i)
                    ret = self.pool.apply_async(func, [], args, error_callback=self.errorlog, callback=self.addResult)
                else:
                    ret = self.pool.apply_async(func, arguments, error_callback=self.errorlog, callback=self.addResult)
            else:
                ret = self.pool.apply_async(func, [], kwargs, error_callback=self.errorlog, callback=self.addResult)
        return ret

    def __del__(self):
        pass

    def getResults(self):
        """
            **Finalise the processing and get the results**

        :return: A dict of results
        """

        self.pool.close()
        self.pool.join()
        self.pool.terminate()
        return self._rets

    def addResult(self, result):
        """
            **This function is used by the threads to return results**

        :param result: The result to add
        :return: None
        """

        self.lock.acquire()
        self._rets.append(result)
        self.lock.release()

    def errorlog(self, e):
        """

        :param e: Error in multiprocessing
        :return: None
        """

        print(e)


##########################################################

lock = Lock()


class MultiProcess():
    def __init__(self, klass):
        self.pool = {}
        self.lock = None

    def initcall(self, key, klass):
        """
            **Initialises the Multi-processing**

        :param key: A unique-key to each calculation
        :param klass: Calls where the function to call is
        :return: None
        """
        self.pool[key] = Thread(klass)
        if self.lock is None:
            self.lock = self.pool[key].lock
        return key

    def call(self, key, klass, function, arguments, **kwargs):
        """
            **Calls a function for multi-processing**

        :param key: The unique-key of the calculation
        :param klass: Calls where the function to call is
        :param function: Name of the function to call
        :param arguments: Arguments of the function to call
        :param kwargs: Keyword arguments passed to the function to call
        :return: None if multi-processing support is activated, or the result of the call otherwise
        """

        res = None
        if multithrds:  # is not None:
            self.pool[key].call(function, arguments, **kwargs)
        else:
            func = getattr(klass, function)
            if not callable(func):
                func = klass.__getattribute__(function)
                res = klass.func(arguments, **kwargs)
            else:
                res = func(arguments, **kwargs)
        return res

    def getResults(self, key):
        """
            **Obtains the results for a calculation.**

        :param key: A unique-key identifying the calculation
        :return: The results of the previous call

        """

        ret = {}
        pool = self.pool.get(key)
        if pool is not None:
            res = pool.getResults()
            for item in res:
                ret.update({key: item})
        return ret


def key():
    """

    :return: a random number to be used as key for multiprocessing
    """
    return random.randint(0, 10000000)


def callThread(klass, function, arguments, **kwargs):
    """
        **Auxiliary function for calling a function with a single Thread**

        The use of Multiprocessing class is preferred, as its faster

    :param klass: Class where the function to call is
    :param function: Name of the function to call
    :param arguments: Arguments of the function to call
    :param kwargs: Keyword arguments passed to the function to call
    :return: None if multi-processing support is activated, or the result of the call otherwise
    """

    if multithrds:  # is not None:
        pool = Thread(klass)
        pool.call(function, arguments, **kwargs)
        tdict = pool.getResults()
        del pool
        return tdict
    else:
        func = klass.__getattribute__(function)
        res = func(**kwargs)
    return res


def call(klass, function, arguments, **kwargs):
    """
        **Auxiliary function for calling a function with a Multiprocessing**


    :param klass: Class where the function to call is
    :param function: Name of the function to call
    :param arguments: Arguments of the function to call
    :param kwargs: Keyword arguments passed to the function to call
    :return: the keyid and a pool for getting the results
    """
    pool = MultiProcess(klass)
    k = pool.initcall(key(), klass)
    pool.call(k, klass, function, arguments, **kwargs)
    return k, pool
