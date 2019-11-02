# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Python Quantities Support Class
"""

import sys
from itertools import chain

import numpy as np
import quantities as pq


class QuantityClass(object):
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

    def __init__(self, Quantities=True):
        """

        :param Quantities: (De)activates python-quantities support

        """
        self.Quantities = Quantities

    def isActive(self):
        """

        :returns: if the Quantity support is active
        """
        return self.Quantities

    def setActive(self, Active=True):
        """

        :param Active: sets new state for the Quantity support
        :return: None
        """
        self.Quantities = Active

    def isDefaultQuantityKey(self, key: str) -> bool:
        """

        :param key: string indicating the plot parameter to check
        :return: if exist a default unit for the plot parameter
        """
        return self.Quantities and self.DefaultUnits.get(key) is not None

    def createDefaultQuantity(self, key, value):
        """
        :param key: One of the keys defined in DefaultUnits
        :param value: The magnitude measured
        :return: The value with the default units
        """
        if not self.Quantities:
            return value
        unit = self.DefaultUnits.get(key)
        if unit is not None:
            value = pq.Quantity(value, unit)

        return value

    def returnQuantity(self, param, unitKey=None, **kwargs):
        """

            :param param: The argument to be evaluated if the Quatities support is activated
            :param unitKey: The key in Default Units.
            :param kwargs: The argument keywords in order to check unit rescaling.
            :return: The Quantity or param if the Quantities support is not activated
        """
        if not self.Quantities:
            return param
        unit = self.DefaultUnits.get(unitKey)

        if not unit or type(param) is pq.Quantity:
            return param
        else:
            return pq.Quantity(param, unit)

    def createQuantityList(self):
        """

        :return: An empty list if Quantities support is activated or
                 a empty numpy.array otherwise
        """
        return [] if self.Quantities else np.array([])

    def appendQuantity(self, vals, val):
        """

        :param vals: A list of Quantities.
        :param val: A Quantity.
        :return: A new list with val appended to vals.
        """
        vals[len(vals):] = [val]
        return vals

    def rescaleFromKey(self, qtylist, units):
        """

        :param qtylist: The input Quantity-like
        :param units: The units to rescale
        :return: The input Quantity-like rescaled to the intented units
        """
        if not self.Quantities or not units or not qtylist: return qtylist
        if type(qtylist) is pq.Quantity:
            try:
                return qtylist.rescale(units)
            except:
                raise BaseException(sys.exc_info()[1])

        ret = self.createQuantityList()
        if units:
            for enum in enumerate(qtylist):
                for qty in enumerate(enum[1]):
                    try:
                        ret = self.appendQuantity(ret, qty[1].rescale(units))
                    except:
                        raise BaseException(sys.exc_info()[1])
        return ret

    def getQuantityUnits(self, qtylist):
        """

        :param qtylist: The Quantity-like to obtain the units from
        :return: A string with the proper units in latex format for easy plotting
        """
        ret = None
        if qtylist:
            if type(qtylist) is pq.Quantity:
                return qtylist.dimensionality.latex
            elif type(qtylist) is list:
                if len(qtylist):
                    if len(qtylist[0]):
                        return qtylist[0][0].dimensionality.latex
        return ret

    def toQuantity(self, listOfQuantities):
        if type(listOfQuantities) is pq.Quantity: return listOfQuantities

        ret = listOfQuantities
        if type(ret) is list:
            templist = list(chain.from_iterable(ret))
            if len(templist):
                ret = pq.Quantity(templist, templist[0].units)
            else:
                ret = pq.Quantity(templist)
        else:
            ret = pq.Quantity(ret)
        return ret

    def Divide(self, Dividend, Divisor):
        """
        :param Dividend: the Quantity-like to be divided.
        :param Divisor:  the Quantity-like to divide by.

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
