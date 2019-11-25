# -*- coding: utf-8 -*-
"""

@author: dragc

"""
import random
from multiprocessing import pool, Lock

from PyGFETdb import multithrds, numprocs

lock = Lock()


class Thread(pool.ThreadPool):
    def __init__(self, package, processes=numprocs):
        """

        :param package: Module where the function to process is

        """

        super(Thread, self).__init__(processes)
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
        del self._rets
        self.pool.close()
        self.pool.join()
        self.pool.terminate()

    def getResults(self):
        """
            **Finalise the processing and get the results**

        :return: A dict of results
        """
        return self._rets

    def addResult(self, result):
        """
            **This function is used by the threads to return results**

        :param result: The result to add
        :return: None
        """

        # self.lock.acquire()
        self._rets.append(result)
        # self.lock.release()

    def errorlog(self, e):
        """

        :param e: Error in multiprocessing
        :return: None
        """

        print(e)

    def end(self):
        self.pool.close()
        self.pool.join()
        self.pool.terminate()


########################################################################
#
# MULTIPROCESSING
#
########################################################################
class MultiProcess(object):
    def __init__(self, klass, processes=numprocs):
        self.pool = Thread(klass, processes)
        self.tasks = {}

    def __del__(self):
        del self.pool
        del self.tasks

    def initcall(self, key, klass):
        """
            **Initialises the Multi-processing**

        :param key: A unique-key to each calculation
        :param klass: Calls where the function to call is
        :return: None
        """
        # if self.lock is None:
        #    self.lock = self.pool.lock
        self.tasks[key] = {}
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
            ret = self.pool.call(function, arguments, **kwargs)
            self.tasks[key].update({key: ret})
            return ret
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
        temp = []
        res = self.pool.getResults()
        tasks = self.tasks.get(key)
        if tasks is not None:
            for t in tasks.values():
                temp.append(t.get())
        ret.update({key: temp})
        return ret

    def end(self):
        self.pool.end()


########################################################################
#
# MULTIPROCESSING UTILITY FUNCTIONS
#
########################################################################
def key():
    """

    :return: a random number to be used as key for multiprocessing
    """
    return random.randint(0, int(1e22))


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
