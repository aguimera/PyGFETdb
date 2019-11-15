# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:27:23 2016

@author: aguimera
"""
from PyGFETdb import QuantitySupport as quantities

qty = quantities.QuantitySupport()  # Activates Quantities support globally

"""
# qty = quantities.QuantitySupport(False)   # Deactivates Quantities support globally

However, the preferable method to deactivate Quantity support is calling *locally* 
the method qty.setActive(False)

"""

multithrds = True
