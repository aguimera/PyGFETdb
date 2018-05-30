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
from scipy.interpolate import interp2d


class SpecSlot():
    def __init__(self, Signal, Units=None, Position=None, DispName=None,):
        self.Fres = 5.0
        self.TimeRes = 0.2
        self.Fmin = 1
        self.Fmax = 150
        self.MinPSD = None
        self.MaxPSD = None
        self.MaxPSDrange = None
        self.LogNormalize = True
        self.Cmap = 'jet'
        self.LogScale = False

        self.Ax = None
        self.CAx = None
        self.Fig = None

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

        nFFT = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/self.Fres))+1))
        Ts = sig.sampling_period.magnitude
        noverlap = int((Ts*nFFT - self.TimeRes)/Ts)

        f, t, Sxx = signal.spectrogram(sig, sig.sampling_rate,
                                       window='hanning',
                                       nperseg=nFFT,
                                       noverlap=noverlap,
                                       scaling='spectrum',
                                       axis=0)

        finds = np.where((self.Fmin < f) & (f < self.Fmax))[0][1:]
        r, g, c = Sxx.shape
        data = Sxx.reshape((r, c))[finds][:]

        if self.MaxPSD is None:
            MaxPSD = data.max()
        else:
            MaxPSD = self.MaxPSD

        if self.MinPSD is None:
            if self.MaxPSDrange is None:
                MinPSD = 10**(np.log10(MaxPSD)-5)
            else:
                MinPSD = 10**(np.log10(MaxPSD)-self.MaxPSDrange)
        else:
            MinPSD = self.MinPSD

        if self.LogNormalize:
            Norm = colors.LogNorm(MinPSD, MaxPSD)
        else:
            Norm = colors.Normalize(MinPSD, MaxPSD)

        if self.Ax is None:
            self.Fig, self.Ax = plt.subplots()

        x = t + sig.t_start.magnitude
        y = f[finds].magnitude
        img = self.Ax.imshow(data,
                             cmap='jet',
                             norm=Norm,
                             interpolation='quadric',
                             origin='lower',
                             aspect='auto',
                             extent=(np.min(x), np.max(x), np.min(y), np.max(y)))

        cbar = plt.colorbar(img, ax=self.CAx, fraction=0.8)
        cbar.ax.tick_params(length=1, labelsize='x-small')

        su = str(sig.units).split(' ')[-1]
        label = "[{}^2]".format(su)
        cbar.set_label(label, fontsize='x-small')

        if self.LogScale:
            self.Ax.set_yscale('log')


class WaveSlot():
    UnitsInLabel = True

    def __init__(self, Signal, Units=None, Position=None, DispName=None,
                 Color='k', Line='-', Alpha=1, Ylim=None, LineWidth=0.5):
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
        self.LineWidth = LineWidth

        self.Name = self.Signal.Name

    def GetSignal(self, Time, Units=None):
        if Units is not None:
            self.Units = Units

        sig = self.Signal.GetSignal(Time, self.Units)
        self.units = sig.units
        return sig

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
                     linewidth=self.LineWidth,
                     color=self.Color,
                     label=label,
                     alpha=self.Alpha)

        if self.Ylim is not None:
            self.Ax.set_ylim(self.Ylim)


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

        self.Fig, A = plt.subplots(max(Pos) + 1, 2,
                                   sharex=True,
                                   figsize=figsize,
                                   gridspec_kw={'width_ratios': (10, 1)})
        self.Axs = [a[0] for a in A]
        self.CAxs = [a[1] for a in A]

        for ca in self.CAxs:
            ca.axis('off')

        for sl in self.Slots:
            if isinstance(sl, SpecSlot):
                sl.CAx = self.CAxs[sl.Position]
            sl.Ax = self.Axs[sl.Position]
            sl.Fig = self.Fig

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
                Ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))

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

        self.Fig.subplots_adjust(top=1.0,
                                 bottom=0.05,
                                 left=0.1,
                                 right=0.94,
                                 hspace=0.0,
                                 wspace=0.0)
#        self.Fig.tight_layout()

    def AddLegend(self, Ax):
        if isinstance(self.SlotsInAxs[Ax][0], SpecSlot):
            self.SlotsInAxs[Ax][0].Ax.set_ylabel('Frequency [Hz]',
                                                 fontsize=self.LegFontSize)
        elif self.ShowNameOn is None:
            return
        elif self.ShowNameOn == 'Axis':
            lbls = []
            for sl in self.SlotsInAxs[Ax]:
                su = str(sl.Units).split(' ')[-1]
                lbls.append("{} [{}]".format(sl.DispName, su))
            sl.Ax.set_ylabel('\n'.join(set(lbls)),
                             fontsize=self.LegFontSize)
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

    def PlotEvents(self, Times, color='r', alpha=0.5):
        for ax in self.Axs:
            ylim = ax.get_ylim()
            ax.vlines(Times, ylim[0], ylim[1], color=color, alpha=alpha)
        
