# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:27:23 2016

@author: aguimera
"""
import PyGFETdb.QuantityClass as quantities

qty = quantities.QuantityClass()  # Activates Quantities support globally
"""
    The preferred method to deactivate Quantity support is calling *locally* the following method: qty.setActive(False)
"""

# qty = quantities.QuantityClass(False)  # Deactivates Quantities support globally
