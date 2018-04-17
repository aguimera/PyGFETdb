#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:11:53 2017

@author: aguimera
"""

import numpy as np
import matplotlib.pyplot as plt
import quantities as pq
from scipy import signal
import matplotlib.colors as colors
from collections import OrderedDict


class SpecSlot():
    def __init__(self, Signal, Units=None, Position=None, DispName=None,):
        self.Fmin = 0.5
        self.Fmax = 15
        self.TimeRes = 0.05
        self.MinPSD = 1e-13
        self.MaxPSD = 1e-10
        self.Cmap = 'jet'
        self.LogScale = False

        self.Signal = Signal

        self.Units = Units
        self.Position = Position
        if DispName is None:
            self.DispName = self.Signal.Name
        else:
            self.DispName = DispName

    def PlotSignal(self, Time, Units=None):
        if Units is not None:
            self.Units = Units
        sig = self.Signal.GetSignal(Time, Units=self.Units)

        nFFT = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/self.Fmin))+1))
        Ts = sig.sampling_period.magnitude
        noverlap = int((Ts*nFFT - self.TimeRes)/Ts)

        f, t, Sxx = signal.spectrogram(sig, sig.sampling_rate,
                                       window='hanning',
                                       nperseg=nFFT,
                                       noverlap=noverlap,
                                       scaling='density',
                                       axis=0)
        finds = np.where(f < self.Fmax)[0][1:]
        r, g, c = Sxx.shape
        S = Sxx.reshape((r, c))[finds][:]
        cax = self.Ax.pcolor(t*pq.s+sig.t_start, f[finds], S,
                             cmap=self.Cmap,
                             norm=colors.Normalize(self.MinPSD,
                                                   self.MaxPSD))
#        fig = plt.figure()
#        fig.colorbar(cax)
#                self.Fig.colorbar(cax, orientation='horizontal')

#                self.Ax.set_yscale('log')


class WaveSlot():
    UnitsInLabel = True

    def __init__(self, Signal, Units=None, Position=None, DispName=None,
                 Color='k', Line='-', Alpha=1, Ylim=None):
        self.Signal = Signal

        self.Units = Units
        self.Position = Position
        if DispName is None:
            self.DispName = self.Signal.Name
        else:
            self.DispName = DispName

        self.Ax = None
        self.Fig = None

        self.Color = Color
        self.Line = Line
        self.Alpha = Alpha
        self.Ylim = Ylim

        self.Name = self.Signal.Name

    def PlotSignal(self, Time, Units=None):
        if self.Ax is None:
            self.Fig, self.Ax = plt.subplots()
        if Units is not None:
            self.Units = Units

        sig = self.Signal.GetSignal(Time, self.Units)

        if self.UnitsInLabel is True:
            su = str(sig.units).split(' ')[-1]
            label = "{} [{}]".format(self.DispName, su)
        else:
            label = self.DispName

        self.Ax.plot(sig.times,
                     sig,
                     self.Line,
                     color=self.Color,
                     label=label,
                     alpha=self.Alpha)

        if self.Ylim is not None:
            self.Ax.set_ylim(self.Ylim)

    def PlotEvent(self, Time, color=None, alpha=0.5):
        if color is None:
            color = 'r--'
        self.Ax.plot((Time, Time), (1, -1), color, alpha=alpha)

    def GetSignal(self, Time, Units=None):
        return self.Signal.GetSignal(Time, Units)


class PlotSlots():
    LegNlabCol = 4  # Number of labels per col in legend
    LegFontSize = 'xx-small'

    def __init__(self, Slots, ShowNameOn='Axis', figsize=None,
                 ShowAxis=True, AutoScale=True):
        self.ShowNameOn = ShowNameOn  # 'Axis', 'Legend', None
        self.Slots = Slots
        self.ShowAxis = ShowAxis
        self.AutoScale = AutoScale

        Pos = []
        for isl, sl in enumerate(self.Slots):
            if sl.Position is None:
                sl.Position = isl
            Pos.append(sl.Position)

        self.Fig, A = plt.subplots(max(Pos) + 1, 1,
                                   sharex=True,
                                   figsize=figsize)
        if type(A).__module__ == np.__name__:
            self.Axs = []
            for a in A:
                self.Axs.append(a)
        else:
            self.Axs = [A, ]

        for sl in self.Slots:
            sl.Ax = self.Axs[sl.Position]

        self.SlotsInAxs = {}
        for ax in self.Axs:
            sll = []
            for sl in self.Slots:
                if sl.Ax == ax:
                    sll.append(sl)
            self.SlotsInAxs.update({ax: sll})

    def ClearAxes(self):
        for sl in self.Slots:
            while sl.Ax.lines:
                sl.Ax.lines[0].remove()

    def FormatFigure(self):
        if self.ShowAxis:
            for Ax in self.Axs:
    #            Ax.get_yaxis().set_visible(False)
                Ax.get_xaxis().set_visible(False)
                Ax.spines['top'].set_visible(False)
                Ax.spines['right'].set_visible(False)
    #            Ax.spines['left'].set_visible(False)
                Ax.spines['bottom'].set_visible(False)
                Ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

            TimeAx = self.Axs[-1]
            TimeAx.set_xlabel('Time [s]')
            TimeAx.get_xaxis().set_visible(True)
        else:
            for Ax in self.Axs:
                Ax.get_yaxis().set_visible(False)
                Ax.get_xaxis().set_visible(False)
                Ax.spines['top'].set_visible(False)
                Ax.spines['right'].set_visible(False)
                Ax.spines['left'].set_visible(False)
                Ax.spines['bottom'].set_visible(False)

        self.Fig.tight_layout()
        self.Fig.subplots_adjust(hspace=0)

    def AddLegend(self, Ax):
        if self.ShowNameOn is None:
            return
        elif self.ShowNameOn == 'Axis':
            lbls = []
            for sl in self.SlotsInAxs[Ax]:
                su = str(sl.Units).split(' ')[-1]
                lbls.append("{} [{}]".format(sl.DispName, su))
            sl.Ax.set_ylabel('\n'.join(set(lbls)))
        elif self.ShowNameOn == 'Legend':
            handles, labels = Ax.get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))

            nLines = len(by_label)
            nlc = self.LegNlabCol
            if nLines > nlc:
                ncol = (nLines / nlc) + ((nLines % nlc) > 0)
            else:
                ncol = 1
            Ax.legend(by_label.values(), by_label.keys(),
                      loc='best',
                      ncol=ncol,
                      fontsize=self.LegFontSize)

    def PlotChannels(self, Time, Units=None):
        self.ClearAxes()
        for sl in self.Slots:
            sl.PlotSignal(Time, Units=Units)
        if Time is not None:
            sl.Ax.set_xlim(Time[0], Time[1])

        for Ax in self.Axs:
            self.AddLegend(Ax)
            if self.AutoScale:
                Ax.autoscale(enable=True, axis='y')

        self.FormatFigure()

    def PlotEvents(self, Time, (EventRec, EventName), color=None, alpha=0.5):

        Te = EventRec.GetEventTimes(EventName, Time)
        for te in Te:
            for sl in self.Slots:
                sl.PlotEvent(te, color=color, alpha=alpha)

