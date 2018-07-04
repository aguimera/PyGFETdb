#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 00:22:59 2018

@author: aguimera
"""
import PhyREC.SignalProcess as RPro
import types
import inspect

Functions = [RPro.__dict__.get(a) for a in dir(RPro)
             if isinstance(RPro.__dict__.get(a), types.FunctionType)]


for f in Functions:
    print f.__code__.co_name
    print f.__code__.co_argcount
    for arg in inspect.getargspec(f)[0]:
        print arg

    