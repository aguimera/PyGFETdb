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

    def PlotResponse(self, sos, Fs):
        if self.Fig is None:
            self.InitFigure()

        # w, h = signal.freqz(b, a,
        w, h = signal.sosfreqz(sos,
                               worN=self.nFFT,
                               fs=2*np.pi*Fs)
        ff = (w/(2*np.pi))[1:]
        self.AxM.semilogx(ff, np.abs(h[1:]))
        self.AxP.semilogx(ff, np.unwrap(np.angle(h)[1:]), '--')
        plt.show()


class ColorBarPlot():
    ImgDicts = {}

    def __init__(self):
        self.Fig = None

    def GenColorBars(self):
        self.Fig, Ax = plt.subplots(1, len(self.ImgDicts))
        try:
            ax = Ax[0]
        except:
            Ax = [Ax, ]

        for iax, (name, img) in enumerate(self.ImgDicts.items()):
            plt.colorbar(img, cax=Ax[iax])
            Ax[iax].set_title(name)

