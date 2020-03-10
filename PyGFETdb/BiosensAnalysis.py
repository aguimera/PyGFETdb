#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:01:03 2020

@author: aguimera
"""

import copy

def GetStepData(Data, Parkwargs):
    Data2Norm = copy.deepcopy(Data)

    for devn, devD in Data2Norm.items():        
        for trtn, trtD in devD.items():            
            for fn, fnD in trtD.items():
                dat = fnD.Get(**Parkwargs)
                trtD[fn] = dat
    
    return Data2Norm