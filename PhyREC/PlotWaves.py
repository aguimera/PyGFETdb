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
from matplotlib.widgets import Slider, Button


def DrawBarScale(Ax, Location='Bottom Left',
                 xsize=None, ysize=None, xoff=0.1, yoff=0.1,
                 xlabelpad=-0.04, ylabelpad=-0.04,
                 xunit='sec', yunit='mV', LineWidth=5, Color='k',
                 FontSize=None):

    # calculate length of the bars
    xmin, xmax, ymin, ymax = Ax.axis()
    AxTrans = Ax.transAxes
    if xsize is None:
        xsize = (xmax - xmin)/5
        xsize = int(np.round(xsize, 0))
    if ysize is None:
        ysize = (ymax - ymin)/5
        ysize = int(np.round(ysize, 0))
    xlen = 1/((xmax - xmin)/xsize)  # length in axes coord
    ylen = 1/((ymax - ymin)/ysize)

    # calculate locations
    if Location == 'Bottom Rigth':
        xoff = 1 - xoff
        ylabelpad = - ylabelpad
        xlen = - xlen
    elif Location == 'Top Left':
        yoff = 1 - yoff
        ylen = - ylen
        xlabelpad = -xlabelpad
    elif Location == 'Top Rigth':
        xoff = 1 - xoff
        ylabelpad = - ylabelpad
        xlen = - xlen
        yoff = 1 - yoff
        ylen = - ylen
        xlabelpad = -xlabelpad
    xdraw = xoff + xlen
    ydraw = yoff + ylen

    # Draw lines
    Ax.hlines(yoff, xoff, xdraw,
              Color,
              linewidth=LineWidth,
              transform=AxTrans,
              clip_on=False)

    Ax.text(xoff + xlen/2,
            yoff + xlabelpad,
            str(xsize) + ' ' + xunit,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=FontSize,
            transform=AxTrans)

    Ax.vlines(xoff, yoff, ydraw,
              Color,
              linewidth=LineWidth,
              transform=AxTrans,
              clip_on=False)

    Ax.text(xoff + ylabelpad,
            yoff + ylen/2,
            str(ysize) + ' ' + yunit,
            horizontalalignment='center',
            verticalalignment='center',
            rotation='vertical',
            fontsize=FontSize,
            transform=AxTrans)


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

        self.units = Units
        self.Position = Position
        if DispName is None:
            if Signal is not None:
                self.DispName = self.Signal.name
        else:
            self.DispName = DispName

    def GetSignal(self, Time, Units=None):
        if Units is None:
            _Units = self.units
        else:
            _Units = Units
        sig = self.Signal.GetSignal(Time, _Units)
        self.units = sig.units
        return sig

    def PlotSignal(self, Time, Units=None):
        sig = self.GetSignal(Time, Units)

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
        cbar.ax.tick_params(length=1, labelsize='xx-small')

        su = str(sig.units).split(' ')[-1]
        label = "[{}^2]".format(su)
        cbar.set_label(label, fontsize='xx-small')

        if self.LogScale:
            self.Ax.set_yscale('log')


class WaveSlot():

    def __init__(self, Signal, Units=None, Position=None, DispName=None,
                 Color='k', Line='-', Alpha=1, Ylim=None, LineWidth=0.5,
                 Ax=None, Fig=None, UnitsInLabel=True, clip_on=True):
        self.Signal = Signal

        self.units = Units
        self.Position = Position
        if DispName is None:
            self.DispName = self.Signal.name
        else:
            self.DispName = DispName

        self.Ax = Ax
        self.Fig = Fig

        self.Color = Color
        self.Line = Line
        self.Alpha = Alpha
        self.Ylim = Ylim
        self.LineWidth = LineWidth
        self.clip_on = clip_on

        self.Name = self.Signal.name
        self.name = self.Signal.name

        self.UnitsInLabel = UnitsInLabel

    def GetSignal(self, Time, Units=None):
        if Units is None:
            _Units = self.units
        else:
            _Units = Units
        sig = self.Signal.GetSignal(Time, _Units)
        self.units = sig.units
        return sig

    def PlotSignal(self, Time, Units=None):
        if self.Ax is None:
            self.Fig, self.Ax = plt.subplots()

        sig = self.GetSignal(Time, Units)

        if self.UnitsInLabel is True:
            su = str(sig.units).split(' ')[-1]
            label = "{} [{}]".format(self.DispName, su)
        else:
            label = self.DispName

        self._PlotSignal(sig, label)

    def _PlotSignal(self, sig, label):
        self.Ax.plot(sig.times,
                     sig,
                     self.Line,
                     linewidth=self.LineWidth,
                     color=self.Color,
                     label=label,
                     alpha=self.Alpha,
                     clip_on=self.clip_on)

        if self.Ylim is not None:
            self.Ax.set_ylim(self.Ylim)

    def CalcAvarage(self, TimeAvg, TimesEvent, Units=None, PltStd=False, 
                    StdAlpha=0.2, PlotTrials=False, TrialsColor='k', TrialsAlpha=0.01):
        avsig = self.GetSignal(None, Units)
        avg = np.array([])

        Ts = avsig.sampling_period
        nSamps = int((TimeAvg[1]-TimeAvg[0])/Ts)
        t = np.arange(nSamps)*Ts + TimeAvg[0]

        for et in TimesEvent:
            start = et+TimeAvg[0]
            stop = et+TimeAvg[1]

            st = np.array(avsig.GetSignal((start, stop))[:nSamps])
            try:
                avg = np.hstack([avg, st]) if avg.size else st
                if PlotTrials:
                    self.Ax.plot(t, st,
                                 color=TrialsColor,
                                 alpha=TrialsAlpha,
                                 clip_on=self.clip_on)
            except:
                print 'Error', nSamps, et, avg.shape, st.shape

        MeanT = np.mean(avg, axis=1)

        MeanTsig = avsig.duplicate_with_new_array(signal=MeanT*avsig.units)
        MeanTsig.t_start = TimeAvg[0]
        MeanTsig.name = MeanTsig.name + ' Avg'
        self._PlotSignal(MeanTsig, label=self.DispName + ' Avg')

        if PltStd:
            StdT = np.std(avg, axis=1)
            self.Ax.fill_between(t, MeanT+StdT, MeanT-StdT,
                                 alpha=StdAlpha,
                                 facecolor=self.Color,
                                 edgecolor=None,
                                 clip_on=self.clip_on)

        ylim = self.Ax.get_ylim()
        self.Ax.vlines((0,), ylim[0], ylim[1], 'r', 'dashdot', alpha=0.5)

class ControlFigure():

    def __init__(self, pltSL, figsize=(20*0.394, 5*0.394)):

        self.pltSL = pltSL
        
        TMax = np.max([sl.Signal.t_stop for sl in pltSL.Slots])
        TMin = np.min([sl.Signal.t_start for sl in pltSL.Slots])
        
        self.Fig, ax = plt.subplots(2, 1, figsize=figsize)
        self.sTstart = Slider(ax[0],
                              label='TStart [s]',
                              valmax=TMax,
                              valmin=TMin,
                              valinit=TMin)

        self.sTshow = Slider(ax[1],
                             label='TShow [s]',
                             valmax=TMax-TMin,
                             valmin=0,
                             valinit=(TMax-TMin)/10)

        self.sTshow.on_changed(self.Update)
        self.sTstart.on_changed(self.Update)

    def Update(self, val):
        twind = (self.sTstart.val * pq.s,
                 self.sTstart.val * pq.s + self.sTshow.val * pq.s)
        self.pltSL.PlotChannels(twind)
        self.pltSL.Fig.canvas.draw()


class PlotSlots():
    LegNlabCol = 4  # Number of labels per col in legend
    LegFontSize = 'xx-small'
    ScaleBarKwargs = {'Location': 'Bottom Left',
                      'xsize': None,
                      'ysize': None,
                      'xoff': 0.1,
                      'yoff': 0.1,
                      'xlabelpad': -0.04,
                      'ylabelpad': -0.04,
                      'xunit': 'sec',
                      'yunit': None,
                      'LineWidth': 5,
                      'Color': 'k',
                      'FontSize': None}

    def __init__(self, Slots, ShowNameOn='Axis', figsize=None,
                 ShowAxis='All', AutoScale=True, Fig=None,
                 ScaleBarAx=None, LiveControl=False):

        self.ShowNameOn = ShowNameOn  # 'Axis', 'Legend', None
        self.Slots = Slots
        self.ShowAxis = ShowAxis      # 'All', int, None
        self.AutoScale = AutoScale
        self.ScaleBarAx = ScaleBarAx

        for sl in self.Slots:
            sig = sl.Signal
            sl.Signal = sig.GetSignal(None)

        if LiveControl:
            self.CtrFig = ControlFigure(self)

        if Fig is not None:
            self.Fig = Fig
            self.Axs = []
            self.CAxs = []
            for sl in self.Slots:
                self.Axs.append(sl.Ax)
                sl.Ax.set_facecolor('#FFFFFF00')
            self.SortSlotsAx()
            return

        Pos = []
        for isl, sl in enumerate(self.Slots):
            if sl.Position is None:
                sl.Position = isl
            Pos.append(sl.Position)

        self.Fig, A = plt.subplots(max(Pos) + 1, 2,
                                   sharex=True,
                                   figsize=figsize,
                                   gridspec_kw={'width_ratios': (10, 1)})
        if len(A.shape) == 1:
            A = A[:, None].transpose()
        self.Axs = [a[0] for a in A]
        self.CAxs = [a[1] for a in A]

        for ca in self.CAxs:
            ca.axis('off')

        for sl in self.Slots:
            if isinstance(sl, SpecSlot):
                sl.CAx = self.CAxs[sl.Position]
            sl.Ax = self.Axs[sl.Position]
            sl.Fig = self.Fig
            sl.Ax.set_facecolor('#FFFFFF00')
# Chech if ax is empty here
        self.SortSlotsAx()

    def SortSlotsAx(self):
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
        if self.ShowAxis == 'All':
            for Ax in self.Axs:
                Ax.get_xaxis().set_visible(False)
                Ax.spines['top'].set_visible(False)
                Ax.spines['right'].set_visible(False)
                Ax.spines['bottom'].set_visible(False)
                if Ax.yaxis.get_scale() == 'linear':
                    Ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))
        elif self.ShowAxis is None:        
            for Ax in self.Axs:
                Ax.get_yaxis().set_visible(False)
                Ax.get_xaxis().set_visible(False)
                Ax.spines['top'].set_visible(False)
                Ax.spines['right'].set_visible(False)
                Ax.spines['left'].set_visible(False)
                Ax.spines['bottom'].set_visible(False)
        else:
            for nAx, Ax in enumerate(self.Axs):
                Ax.spines['top'].set_visible(False)
                Ax.spines['bottom'].set_visible(False)
                Ax.spines['right'].set_visible(False)
                if nAx == self.ShowAxis:
                    continue
                Ax.get_yaxis().set_visible(False)
                Ax.get_xaxis().set_visible(False)
                Ax.spines['left'].set_visible(False)

        if self.ShowAxis is not None:
            TimeAx = self.Axs[-1]
            TimeAx.set_xlabel('Time [s]', fontsize=self.LegFontSize)
            TimeAx.get_xaxis().set_visible(True)

        if self.ScaleBarAx is not None:
            if self.ScaleBarKwargs['yunit'] is None:
                sl = self.SlotsInAxs[self.Axs[self.ScaleBarAx]][0]
                su = str(sl.units).split(' ')[-1]
                self.ScaleBarKwargs['yunit'] = su
            DrawBarScale(self.Axs[self.ScaleBarAx], **self.ScaleBarKwargs)

        if len(self.CAxs) == 0:
            self.Fig.tight_layout()
        else:
            self.Fig.subplots_adjust(top=0.975,
                                     bottom=0.095,
                                     left=0.1,
                                     right=0.95,
                                     hspace=0.0,
                                     wspace=0.0)

    def AddLegend(self, Ax):
        if len(self.SlotsInAxs[Ax]) == 0:
            print 'empty Ax'
            return
        if isinstance(self.SlotsInAxs[Ax][0], SpecSlot):
            self.SlotsInAxs[Ax][0].Ax.set_ylabel('Freq. [Hz]',
                                                 fontsize=self.LegFontSize)
        elif self.ShowNameOn is None:
            return
        elif self.ShowNameOn == 'Axis':
            lbls = []
            for sl in self.SlotsInAxs[Ax]:
                su = str(sl.units).split(' ')[-1]
                lbls.append("{} [{}]".format(sl.DispName, su))
            sl.Ax.set_ylabel('\n'.join(set(lbls)),
                             fontsize=self.LegFontSize,
#                             rotation='horizontal',
                             )
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
                      ncol=ncol,
                      fontsize=self.LegFontSize)

    def PlotChannels(self, Time, Units=None, FormatFigure=True):
        self.ClearAxes()
        for sl in self.Slots:
            sl.PlotSignal(Time, Units=Units)
        if Time is not None:
            if Time[0] is not None:
                sl.Ax.set_xlim(left=Time[0].magnitude)
            if Time[1] is not None:
                sl.Ax.set_xlim(right=Time[1].magnitude)

        for Ax in self.Axs:
            self.AddLegend(Ax)
            if self.AutoScale:
                Ax.autoscale(enable=True, axis='y')

        if FormatFigure:
            self.FormatFigure()

    def PlotEvents(self, Times, color='r', alpha=0.5,
                   Labels=None, lAx=0, fontsize='xx-small', LabPosition='top'):

        self.Texts = []
        if Labels is not None:
            for ilbl, lbl in enumerate(Labels):
                for ax in self.Axs:
                    ylim = ax.get_ylim()
                    ax.vlines(Times[ilbl], ylim[0], ylim[1],
                              color=color,
                              alpha=alpha)
                lax = self.Axs[lAx]
                if LabPosition == 'top':
                    ylim = lax.get_ylim()[1]
                else:
                    ylim = lax.get_ylim()[0]
                txt = lax.text(Times[ilbl], ylim, lbl, fontsize=fontsize)
                self.Texts.append(txt)
            return

        EventLines = []
        for ax in self.Axs:
            ylim = ax.get_ylim()
            lines = ax.vlines(Times, ylim[0], ylim[1], color=color, alpha=alpha)
#            EventLines.append(lines[0])

        return EventLines

    def PlotEventAvarage(self, TimeAvg, TimesEvent, Units=None, PltStd=False, 
                         StdAlpha=0.2,
                         PlotTrials=False, TrialsColor='k', TrialsAlpha=0.01,
                         ClearAxes=True, AvgColor=None):

        if ClearAxes:
            self.ClearAxes()

        for sl in self.Slots:
            if AvgColor is not None:
                sl.Color = AvgColor
            sl.CalcAvarage(TimeAvg, TimesEvent,
                           Units=Units,
                           PlotTrials=PlotTrials,
                           TrialsColor=TrialsColor,
                           TrialsAlpha=TrialsAlpha,
                           PltStd=PltStd,
                           StdAlpha=StdAlpha)

        sl.Ax.set_xlim(left=TimeAvg[0].magnitude)
        sl.Ax.set_xlim(right=TimeAvg[1].magnitude)
        
        self.FormatFigure()
