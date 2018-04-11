#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:11:53 2017

@author: aguimera
"""

import neo
import numpy as np
import matplotlib.pyplot as plt
import quantities as pq
from scipy import signal
from scipy import interpolate
import pickle
from fractions import Fraction
from matplotlib import cm
import matplotlib.colors as colors


def threshold_detection(signal, threshold=0.0 * pq.mV, sign='above',
                        RelaxTime=None):
    """
    Returns the times when the analog signal crosses a threshold.
    Usually used for extracting spike times from a membrane potential.
    Adapted from version in NeuroTools.

    Parameters
    ----------
    signal : neo AnalogSignal object
        'signal' is an analog signal.
    threshold : A quantity, e.g. in mV
        'threshold' contains a value that must be reached
        for an event to be detected. Default: 0.0 * mV.
    sign : 'above' or 'below'
        'sign' determines whether to count thresholding crossings
        that cross above or below the threshold.
    format : None or 'raw'
        Whether to return as SpikeTrain (None)
        or as a plain array of times ('raw').

    Returns
    -------
    result_st : neo SpikeTrain object
        'result_st' contains the spike times of each of the events (spikes)
        extracted from the signal.
    """

    assert threshold is not None, "A threshold must be provided"

    if sign == 'above':
        cutout = np.where(signal > threshold)[0]
    elif sign == 'below':
        cutout = np.where(signal < threshold)[0]

    if len(cutout) <= 0:
        events = np.zeros(0)
    else:
        take = np.where(np.diff(cutout) > 1)[0] + 1
        take = np.append(0, take)

        time = signal.times
        events = time[cutout][take]

    if RelaxTime:
        outevents = []
        told = 0*pq.s
        for te in events:
            if (te-told) > RelaxTime:
                outevents.append(te)
                told = te
    else:
        outevents = events

    return outevents


class SpecSlot():
    def __init__(self):
        self.SpecFmin = 0.5
        self.SpecFmax = 15
        self.SpecTimeRes = 0.05
        self.SpecMinPSD = -3
        self.SpecCmap = 'jet'
        if self.PlotType == 'Spectrogram':
            sig = self.GetSignal(Time, Resamp=False)

            nFFT = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/self.SpecFmin))+1))
            Ts = sig.sampling_period.magnitude
            noverlap = int((Ts*nFFT - self.SpecTimeRes)/Ts)

            f, t, Sxx = signal.spectrogram(sig, sig.sampling_rate,
                                           window='hanning',
                                           nperseg=nFFT,
                                           noverlap=noverlap,
                                           scaling='density',
                                           axis=0)
            finds = np.where(f < self.SpecFmax)[0][1:]
            r, g, c = Sxx.shape
            S = Sxx.reshape((r, c))[finds][:]
            cax = self.Ax.pcolor(t*pq.s+sig.t_start, f[finds], S,
                                 cmap=self.SpecCmap,
                                 norm=colors.LogNorm(self.SpecMinPSD,
                                                     self.SpecMaxPSD))
            fig = plt.figure()
            fig.colorbar(cax)
#                self.Fig.colorbar(cax, orientation='horizontal')

#                self.Ax.set_yscale('log')


class WaveSlot():
    def __init__(self, Signal, Units=None, Position=None, DispName=None,
                 Color='k', Line='-', Alpha=1, Ymax=0, Ymin=0, AutoScale=True):
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
        self.Ymax = Ymax
        self.Ymin = Ymin
        self.AutoScale = AutoScale

    def PlotSignal(self, Time, Units=None):
        if self.Ax is None:
            self.Fig, self.Ax = plt.subplots()
        if Units is not None:
            self.Units = Units

        sig = self.Signal.GetSignal(Time, self.Units)
        self.Ax.plot(sig.times,
                     sig,
                     self.Line,
                     color=self.Color,
                     label=self.DispName,
                     alpha=self.Alpha)

        if self.AutoScale:
            ylim = (np.min(sig), np.max(sig))
            self.Ax.set_ylim(ylim)
        else:
            ylim = (self.Ymin, self.Ymax)
            self.Ax.set_ylim(ylim)

    def PlotEvent(self, Time, color=None, alpha=0.5):
        if color is None:
            color = 'r--'
        self.Ax.plot((Time, Time), (1, -1), color, alpha=alpha)


class PlotSlots():
    LegNlabCol = 4  # Number of labels per col in legend
    LegFontSize = 'x-small'

    def __init__(self, Slots, ShowLegend=False, figsize=None,
                 ShowAxis=True):
        self.ShowLegend = ShowLegend
        self.Slots = Slots
        self.ShowAxis = ShowAxis

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

        for sl in self.Slots:  # init slots axis
            sl.Ax = self.Axs[sl.Position]
            sl.Fig = self.Fig

            if not ShowLegend:  # Add labels to axis
                lb = sl.Ax.get_ylabel()
                sl.Ax.set_ylabel(lb + ' ' + sl.DispName + '\n',
                                 rotation='horizontal',
                                 ha='center')
        self.FormatFigure()

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
        nLines = len(Ax.lines)
        nlc = self.LegNlabCol
        if nLines > nlc:
            ncol = (nLines / nlc) + ((nLines % nlc) > 0)
        else:
            ncol = 1
        Ax.legend(loc='best',
                  fancybox=True,
                  shadow=True,
                  ncol=ncol,
                  fontsize=self.LegFontSize)

    def PlotChannels(self, Time, Units=None):
        self.ClearAxes()
        for sl in self.Slots:
            sl.PlotSignal(Time, Units=Units)

        if Time is not None:
            sl.Ax.set_xlim(Time[0], Time[1])

        if self.ShowLegend:
            for (Ax, AutoScale) in self.Axs:
                self.AddLegend(Ax)
                if AutoScale:
                    Ax.autoscale(enable=True, axis='y', tight=True)

    def PlotEvents(self, Time, (EventRec, EventName), color=None, alpha=0.5):

        Te = EventRec.GetEventTimes(EventName, Time)
        for te in Te:
            for sl in self.Slots:
                sl.PlotEvent(te, color=color, alpha=alpha)

