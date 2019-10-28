#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.
"""

import numpy as np


def Divide(Dividend, Divisor):
    """
    :param Dividend: the array-like to be divided.
    :param Divisor:  the array-like to divide by.

    :returns: Divides Dividend by Divisor avoiding division by zero.
    """

    ret = np.divide(Dividend, Divisor)
    if type(ret) is np.ndarray:
        return np.where(np.isfinite(ret), ret, Dividend)
    if np.isfinite(ret):
        return ret
    else:
        return Dividend
