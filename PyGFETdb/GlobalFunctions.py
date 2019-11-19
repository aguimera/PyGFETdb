# -*- coding: utf-8 -*-
"""
@author: dragc 25/10/19 22:48.

Global Functions that do not fit in the previous files.

"""

import numpy as np


def updateDictOfLists(dict, key, value):
    """
        **Modifies a dictionary of lists, appending the value at the list obtained
        of applying the key to the dictionary**

    :param dict: A dictionary to update
    :param key: The key to update
    :param value: The value to append
    :return: None
    """
    k = dict.get(key)
    if k is None:
        dict[key] = [value]
    else:
        k.append(value)

def remove(Valx, index):
    """
        **Removes the value at the index specified from an array**

    :param Valx: array to remove the value
    :param index: index to remove the value
    :return: the array without the value at the index specified
    """
    if type(Valx) is list:
        Valx.remove(Valx[index])
    else:
        Valx = Valx.tolist()
        Valx.remove(Valx[index])
        Valx = np.array(Valx)
    return Valx
