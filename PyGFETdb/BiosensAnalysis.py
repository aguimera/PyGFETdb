#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:01:03 2020

@author: aguimera
"""

import copy

def GetStepData(Data, Parkwargs, RefSteps=None):
    Data2Norm = copy.deepcopy(Data)
    
    for devn, devD in Data2Norm.items():        
        for trtn, trtD in devD.items():            
            for fn, fnD in trtD.items():
                if type(fnD)==list:
                    trtD[fn] = fnD
                    continue
                dat = fnD.Get(**Parkwargs)
                trtD[fn] = dat.flatten()[0]
    
    if RefSteps is not None:
        for devn, devD in Data2Norm.items():
            ToRemove = []
            for trtn, trtD in devD.items():
                Ref = None
                for fn, fnD in trtD.items():
                    if fn in RefSteps:
                        Ref = fnD                    
                if Ref is None:
                    print('Ref not found ->>', trtn, list(trtD.keys()))
                    ToRemove.append(trtn)
                    continue
                for fn, fnD in trtD.items():
                    if type(fnD)==list:
                        trtD[fn] = fnD
                        continue
                    trtD[fn] = fnD - Ref
            
            if len(ToRemove):
                print('\n ************ REMOVING !!!!!!!!')
                print( ToRemove, '\n')
                for rm in ToRemove:
                    devD.pop(rm)

    return Data2Norm