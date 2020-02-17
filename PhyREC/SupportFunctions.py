#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 20:23:10 2020

@author: aguimera
"""

import matplotlib.pyplot as plt
from scipy import signal
import numpy as np


class DebugFilterPlot():
    nFFT = 2**15
    
    def __init__(self):
        self.Fig = None
        
    def InitFigure(self):
        self.Fig, self.AxM = plt.subplots()
        self.AxP = plt.twinx(self.AxM)        

    def PlotResponse(self, a, b, Fs):
        if self.Fig is None:
            self.InitFigure()
        
        w, h = signal.freqz(b, a,
                            worN=self.nFFT,
                            fs=2*np.pi*Fs)
        ff = (w/(2*np.pi))[1:]
        self.AxM.semilogx(ff, np.abs(h[1:]))
        self.AxP.semilogx(ff, np.unwrap(np.angle(h)[1:]),'--')
        plt.show()


class ColorBarPlot():
    
    def __init__(self):
        self.Fig = None
        
    def InitFigure(self):
        self.Fig, self.AxM = plt.subplots()
        self.AxP = plt.twinx(self.AxM)        

    def PlotResponse(self, a, b, Fs):
        if self.Fig is None:
            self.InitFigure()
        
        w, h = signal.freqz(b, a,
                            worN=self.nFFT,
                            fs=2*np.pi*Fs)
        ff = (w/(2*np.pi))[1:]
        self.AxM.semilogx(ff, np.abs(h[1:]))
        self.AxP.semilogx(ff, np.unwrap(np.angle(h)[1:]),'--')
        plt.show()