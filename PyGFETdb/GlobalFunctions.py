# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.
"""

import numpy as np
import quantities as pq

Quantities = True
# Quantities = False
"""
 Activate python-quantities support
"""

DefaultUnits = {'Vds': pq.V,
                'Ud0': pq.V,
                'PSD': pq.A ** 2 / pq.Hz,
                'Fgm': pq.S,
                'gm': pq.S,
                'gmV': pq.S,
                'Vgs': pq.V,
                'Fpsd': pq.Hz,
                'Ig': pq.A,
                'Irms': pq.V,
                'Vrms': pq.V / pq.S,
                'Ids': pq.A,
                'Rds': pq.V / pq.A
                }
"""
 Dictionary with the association between PlotParameters and Units
"""


def isQuantityKey(key):
    return DefaultUnits.get(key) is not None


def createDefaultQuantity(key, value):
    """
    :param key: One of the keys defined in DefaultUnits
    :param value: The magnitude measured
    :return: The value with the default units
    """
    if not Quantities:
        return value
    unit = DefaultUnits.get(key)
    if unit is not None:
        value = pq.Quantity(value, unit)
    else:
        print('Error: No default unit for:', key)
        value = pq.Quantity(value, pq.V)
    return value


def returnQuantity(param, unitKey=None, **kwargs):
    """

    :param param: The argument to be evaluated if the Quatities support is activated
    :param unitKey: The key in Default Units.
    :param kwargs: The argument keywords in order to check unit rescaling.
    :return: The Quantity (rescaled if indicated in kwargs) or param if the Quantities
             support is not activated
    """
    if not Quantities:
        return param
    unit = DefaultUnits.get(unitKey)

    if type(param) is pq.Quantity or not unit:
        return param
    else:
        return pq.Quantity(param, unit)


def createQuantityList():
    """

    :return: An empty list if Quantities support is activated or
                 a empty numpy.array otherwise
    """
    return [] if Quantities else np.array([])


def appendQuantity(vals, val):
    """

    :param vals: A list of Quantities.
    :param val: A Quantity.
    :return: A new list with val appended to vals.
    """
    vals[len(vals):] = [val]
    return vals


def rescaleFromKey(qtylist, units):
    """

    :param qtylist: The input list of Quantities
    :param units: The units to rescale
    :param kwargs: Keyword arguments
    :return:
    """
    ret = createQuantityList()
    if qtylist and units:
        for enum in enumerate(qtylist):
            for qty in enumerate(enum[1]):
                ret = appendQuantity(ret, qty[1].rescale(units))
    return ret

def Divide(Dividend, Divisor):
    """
    :param Dividend: the array-like to be divided.
    :param Divisor:  the array-like to divide by.

    :returns: Divides Dividend by Divisor avoiding division by zero.
    """
    units = None
    ret = np.divide(Dividend, Divisor)
    if type(ret) is pq.Quantity:
        units = ret.units
    ret = np.where(np.isfinite(ret), ret, Dividend)

    if units:
        return pq.Quantity(ret, units)
    elif ret.ndim == 0:
        return ret.tolist()
    else:
        return ret
