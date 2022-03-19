#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 20:20:54 2022

@author: aguimera
"""

import datetime
from PyGFETdb.ws_json_class import ws_json
import pickle
import sys


class PyFETdb():
    @staticmethod
    def _execute(query, values):   
        return ws_json.call("pyfet", queryt % values)
    
    
    
if __name__ == '__main__':

   
        
    queryt = """ SELECT Wafers.Name, Trts.Name as Trts_Name, DCcharacts.Comments, DCcharacts.FuncStep, DCcharacts.Ph,
                DCcharacts.Data FROM DCcharacts
                INNER JOIN Trts ON Trts.idTrts = DCcharacts.Trt_id
                INNER JOIN TrtTypes ON TrtTypes.idTrtTypes = Trts.Type_id
                INNER JOIN Devices ON Devices.idDevices = Trts.Device_id
                INNER JOIN Wafers ON Wafers.idWafers = Devices.Wafer_id
                WHERE Trts.Name = '%s' """    

    res = PyFETdb._execute(queryt, ('B12744W3-Xip6NS-Ch12',))

    
    