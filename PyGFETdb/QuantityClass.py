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
                    'Rds': pq.ohm
                    }
    """
    Dictionary with the association between PlotParameters and Units
    """

    def __init__(self, Quantities=True):
        """

        :param Quantities: (De)activates python-quantities support

        """
        self._Active = Quantities

    def isActive(self):
        """

        :returns: if the Quantity support is active
        """
        return self._Active

    def setActive(self, Active=True):
        """

        :param Active: sets new state for the Quantity support
        :return: None
        """
        self._Active = Active

    def isDefaultQuantityKey(self, key: str) -> bool:
        """

        :param key: string indicating the plot parameter to check
        :return: if exist a default unit for the plot parameter
        """
        return self.isActive() and self.DefaultUnits.get(key) is not None

    def createDefaultQuantity(self, key, value):
        """
        :param key: One of the keys defined in DefaultUnits
        :param value: The magnitude measured
        :return: The value with the default units
        """
        if not self.isActive():
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
        if not self.isActive():
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
        return [] if self.isActive() else np.array([])

    def appendQuantity(self, vals, val):
        """

        :param vals: A list of Quantities.
        :param val: A Quantity.
        :return: A new list with val appended to vals.
        """
        if vals is not None:
            vals.append(val)
        else:
            vals = []
        return vals

    def rescaleFromKey(self, qtylist, units):
        """

        :param qtylist: The input Quantity-like
        :param units: The units to rescale
        :return: The input Quantity-like rescaled to the intented units
        """
        if not self.isActive() or units is None or qtylist is None: return qtylist
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
        :return: The proper UnitQuantity
        """
        ret = None
        if self.isActive() and qtylist is not None:
            if type(qtylist) is pq.Quantity or type(qtylist) is pq.UnitQuantity:
                return qtylist.units
            elif type(qtylist) is list:
                if len(qtylist) > 0:
                    ret = self.getQuantityUnits(qtylist[0])

        return ret

    def toQuantity(self, value):
        """

        :param value: The value to be converted, a Quantity-like
        :return: A Quantity with consistent units depending on parameter value
        """
        ret = pq.Quantity(np.nan)
        if self.isActive() and value is not None:
            if type(value) is list:
                ret = self.getConsistentQuantityFromList(value)
            elif type(value) is pq.Quantity:
                ret = value
            else:  # is a number
                ret = pq.Quantity(value)
        return ret

    def getConsistentQuantityFromList(self, listQty: list):
        """

        :param listQty:
        :return: A list of Quantities if all the Quantities have consistent Units,
                 with pq.dimensionless in the case were incompatible,
                 or a 'nan' Quantity if it were not possible to convert them
        """
        ret = listQty

        if self.isActive() and listQty is not None:
            ret = pq.Quantity(np.nan)
            templist = list(chain.from_iterable(listQty))
            if len(templist):
                if type(templist) is list and type(templist[0]) is pq.Quantity:
                    retunits = templist[0].units
                    units = None
                    for i, quan in enumerate(templist):
                        if not units:
                            units = templist[i].units
                        elif quan.units != units:
                            retunits = pq.dimensionless
                    ret = pq.Quantity(templist, retunits)
                else:
                    ret = pq.Quantity(templist)
        return ret

    def flatten(self, quantity):
        """

        :param quantity: A Quantity-like variable
        :return: A flatted list of values with the magnitudes of quantity
        """
        tvals = []
        if type(quantity) is list and len(quantity) > 0:
            if type(quantity[0]) is list:
                Dat = list(chain.from_iterable(quantity))
            else:
                Dat = quantity
            for item in Dat:
                if type(item) is pq.Quantity and item.size > 1 or type(item) is np.ndarray and item.size > 0:
                    if len(item) > 0:
                        for item2 in item:
                            tvals.append(item2)
                    else:
                        tvals.append(item)
                else:
                    if type(item) is list:
                        # if len(item)>0 and type(item[0]) is list:
                        tvals = self.flatten(item)
                        # else:
                        #   tvals.append(item)
                    else:
                        tvals.append(item)
        else:
            tvals = quantity
        return tvals

    def flattenQuantity(self, quantity):
        """
        Flattens a Quantity with the units
        :param quantity: Quantity to flatten
        :return: the Quantity flattened
        """
        units = self.getQuantityUnits(quantity)
        ret = self.flatten(quantity)
        if units is None or units == pq.dimensionless:
            return ret
        else:
            return pq.Quantity(ret, units)

    def Divide(self, Dividend, Divisor):
        """
        :param Dividend: the Quantity-like to be divided.
        :param Divisor:  the Quantity-like to divide by.

        :returns: Divides Dividend by Divisor avoiding division by zero.
        """
        units = None
        try:
            ret = np.divide(Dividend, Divisor)
        except:
            print("Error in Divide:", sys.exc_info())
            raise ArithmeticError

        if type(ret) is pq.Quantity:
            units = ret.units
        ret = np.where(np.isfinite(ret), ret, Dividend)

        if units:
            return pq.Quantity(ret, units)
        elif ret.ndim == 0:
            return ret.tolist()

        return ret
