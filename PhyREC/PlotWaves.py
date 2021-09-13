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
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.artist import ArtistInspector
import PhyREC.SignalProcess as Spro
from . import SpectColBars
from sklearn.impute import KNNImputer
from copy import deepcopy


class ImageSequence(pq.Quantity):
   
    def __new__(cls, signal, sampling_rate,  t_start=0*pq.s, units=None, dtype=None, copy=True,
                name=None, **kwargs):
       
        signal = cls._rescale(cls, signal=signal, units=units)
        obj = pq.Quantity(signal, units=units, dtype=dtype, copy=copy).view(cls)                                         
        
        obj.sampling_rate = sampling_rate
        obj.t_start = t_start
        obj.name = name
        obj.annotations = {}
        return obj
    
    def __array_finalize__(self, obj):
        if obj is None:
            return
        super(ImageSequence, self).__array_finalize__(obj)
        self.sampling_rate = getattr(obj, 'sampling_rate', None)
        self.t_start = getattr(obj, 't_start', None)
        self.name = getattr(obj, 'name', None)
        self.annotations = getattr(obj, 'annotations', {})

    
    def __repr__(self):
        return("{cls} {frame} frames with {width} pixels of width and {height} pixels of height \n"
                "units {units} \n"
                "sampling_rate {sampling_rate} \n"
                "Duration {t_start} --  {t_stop}\n"
                "datatype {dtype} \n"
                "name {name} ".format(cls=self.__class__.__name__,
                                            frame=self.shape[0],
                                            height=self.shape[1],
                                            width=self.shape[2],
                                            units=self.units,
                                            dtype=self.dtype,
                                            sampling_rate=self.sampling_rate,
                                            t_start=self.t_start,
                                            t_stop=self.t_stop,
                                            name=self.name))

    def __deepcopy__(self, memo):
        """
            Create a deep copy of the data object.
            All attributes and annotations are also deep copied.
            References to parent objects are not kept, they are set to None.


            :param memo: (dict) Objects that have been deep copied already
            :return: (DataObject) Deep copy of the input DataObject
        """
        cls = self.__class__
        necessary_attrs = {'signal':self,
                           'units':self.units,
                           'sampling_rate': self.sampling_rate,
                           'name': self.name,
                           't_start': self.t_start}
        # Create object using constructor with necessary attributes
        new_obj = cls(**necessary_attrs)
        # Add all attributes

        return new_obj

    def _rescale(self, signal, units=None):
        '''
        Check that units are present, and rescale the signal if necessary.
        This is called whenever a new signal is
        created from the constructor. See :meth:`__new__' in
        :class:`AnalogSignal` and :class:`IrregularlySampledSignal`
        '''
        if units is None:
            if not hasattr(signal, "units"):
                raise ValueError("Units must be specified")
        elif isinstance(signal, pq.Quantity):
            # This test always returns True, i.e. rescaling is always executed if one of the units
            # is a pq.CompoundUnit. This is fine because rescaling is correct anyway.
            if pq.quantity.validate_dimensionality(units) != signal.dimensionality:
                signal = signal.rescale(units)
        return signal

    def rescale(self, units):
        obj = super(ImageSequence, self).rescale(units)

        obj.sampling_rate = self.sampling_rate
        obj.t_start = self.t_start
        obj.name = self.name
        obj.annotations = self.annotations        
        return obj.view(ImageSequence)
    
    def time_index(self, t):
        """Return the array index corresponding to the time `t`"""
        i = (t - self.t_start) * self.sampling_rate
        i = int(np.rint(i.simplified.magnitude))
        return i

    def time_slice(self, t_start, t_stop):
        '''
        Creates a new AnalogSignal corresponding to the time slice of the
        original AnalogSignal between times t_start, t_stop. Note, that for
        numerical stability reasons if t_start does not fall exactly on
        the time bins defined by the sampling_period it will be rounded to
        the nearest sampling bin. The time bin for t_stop will be chosen to
        make the duration of the resultant signal as close as possible to
        t_stop - t_start. This means that for a given duration, the size
        of the slice will always be the same.
        '''

        # checking start time and transforming to start index
        if t_start is None:
            i = 0
            t_start = 0 * pq.s
        else:
            i = self.time_index(t_start)

        # checking stop time and transforming to stop index
        if t_stop is None:
            j = len(self)
        else:
            delta = (t_stop - t_start) * self.sampling_rate
            j = i + int(np.rint(delta.simplified.magnitude))

        if (i < 0) or (j > len(self)):
            raise ValueError('t_start, t_stop have to be within the analog \
                              signal duration')

        # Time slicing should create a deep copy of the object
        obj = deepcopy(self[i:j, :, :])

        obj.t_start = self.t_start + i * self.sampling_period      

        return obj

    @property
    def duration(self):
        '''
        Signal duration

        (:attr:`size` * :attr:`sampling_period`)
        '''
        return self.shape[0] / self.sampling_rate

    @property
    def t_stop(self):
        '''
        Time when signal ends.

        (:attr:`t_start` + :attr:`duration`)
        '''
        return self.t_start + self.duration

    @property
    def times(self):
        '''
        The time points of each sample of the signal

        (:attr:`t_start` + arange(:attr:`shape`)/:attr:`sampling_rate`)
        '''
        return self.t_start + np.arange(self.shape[0]) / self.sampling_rate
   
    @property
    def sampling_period(self):
        '''
        Interval between two samples.

        (1/:attr:`sampling_rate`)
        '''
        return 1. / self.sampling_rate
    
    def annotate(self, **annotations):
        """
        Add annotations (non-standardized metadata) to a Neo object.

        Example:

        >>> obj.annotate(key1=value0, key2=value1)
        >>> obj.key2
        value2
        """        
        self.annotations.update(annotations)
    
    def __array_ufunc__(self, ufunc, method, *inputs, out=None, **kwargs):
        args = []
        in_no = []
        for i, input_ in enumerate(inputs):
            if isinstance(input_, ImageSequence):
                in_no.append(i)
                args.append(input_.view(pq.Quantity))
            else:
                args.append(input_)

        results = super(ImageSequence, self).__array_ufunc__(ufunc, method, *args, **kwargs)
        if results is NotImplemented:
            return NotImplemented

        results = results.view(ImageSequence)
        results.sampling_rate = inputs[0].sampling_rate
        results.t_start = inputs[0].t_start
        results.name = inputs[0].name
        results.annotations = inputs[0].annotations
        
        return results
        
        # if method == 'at':
        #     if isinstance(inputs[0], A):
        #         inputs[0].info = info
        #     return

        # if ufunc.nout == 1:
        #     results = (results,)

        # results = tuple((np.asarray(result).view(A)
        #                  if output is None else output)
        #                 for result, output in zip(results, outputs))
        # if results and isinstance(results[0], A):
        #     results[0].info = info

        # return results[0] if len(results) == 1 else results
    def duplicate_with_new_data(self, signal, units=None):
        '''
        Create a new signal with the same metadata but different data.
        Required attributes of the signal are used.
        Note: Array annotations can not be copied here because length of data can change
        '''
        if units is None:
            units = self.units
        # else:
        #     units = pq.quantity.validate_dimensionality(units)

        # signal is the new signal
        necessary_attrs = {'signal': signal,
                           'units':self.units,
                           'sampling_rate': self.sampling_rate,
                           'name': self.name,
                           't_start': self.t_start}
        new = self.__class__(**necessary_attrs)
        new.annotations.update(self.annotations)
        # Note: Array annotations are not copied here, because it is not ensured
        # that the same number of signals is used and they would possibly make no sense
        # when combined with another signal
        return new

#from NeoInterface import NeoTrain

def DrawBarScale(Ax, Location='Bottom Left',
                 xsize=None, ysize=None, xoff=0.1, yoff=0.1,
                 xlabelpad=-0.04, ylabelpad=-0.04,
                 xunit='sec', yunit='mV', LineWidth=5, Color='k',
                 FontSize=None, ylabel=None, xlabel=None):

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

    if xlabel is None:
        xlabel = str(xsize) + ' ' + xunit

    Ax.text(xoff + xlen/2,
            yoff + xlabelpad,
            xlabel,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=FontSize,
            transform=AxTrans)

    Ax.vlines(xoff, yoff, ydraw,
              Color,
              linewidth=LineWidth,
              transform=AxTrans,
              clip_on=False)

    if ylabel is None:
        ylabel = str(ysize) + ' ' + yunit

    Ax.text(xoff + ylabelpad,
            yoff + ylen/2,
            ylabel,
            horizontalalignment='center',
            verticalalignment='center',
            rotation='vertical',
            fontsize=FontSize,
            transform=AxTrans)


def UpdateTreeDictProp(obj, prop):
    ains = ArtistInspector(obj)
    validp = ains.get_setters()
    for p in prop.keys():
        if p in validp:
            obj.set(**{p: prop[p]})
        else:
            obj2 = getattr(obj, 'get_' + p)()
            UpdateTreeDictProp(obj2, prop[p])


class SpecSlot():
    """Comment"""
    DefspecKwargs = {
                  'Fmax': 100*pq.Hz,
                  'Fmin': 0.5*pq.Hz,
                  'Fres': 0.5*pq.Hz,
                  'TimeRes': 0.1*pq.s,
                  'Zscored': True,
                 }
    
    DefAvgSpectKwargs = {'SpecArgs': {
                                      'Fmax': 100*pq.Hz,
                                      'Fmin': 0.5*pq.Hz,
                                      'Fres': 0.5*pq.Hz,
                                      'TimeRes': 0.1*pq.s,
                                      'Zscored': False,
                                      },
                         'AvgSpectNorm': 'Zscore',
                         'AvgSpectNormTime': None,
                         }

    DefAxKwargs = {
                   'ylabel': 'Freq [Hz]',
                   'xaxis': {'visible': False,
                             },
                   'yaxis': {'visible': True,
                             },
                   }

    DefImKwargs = {
                    'norm': colors.Normalize(-3, 3),
                    'cmap': 'seismic',
                    'interpolation': 'bilinear',
                    }

    def UpdateLineKwargs(self, LineKwargs):
        pass
#        self.LineKwargs.update(LineKwargs)
#        UpdateTreeDictProp(self.Line, self.LineKwargs)

    def UpdateAxKwargs(self, AxKwargs):
        # pass
        self.AxKwargs.update(AxKwargs)
        UpdateTreeDictProp(self.Ax, self.AxKwargs)

    def __init__(self, Signal, Units=None, Position=None, imKwargs=None,
                 specKwargs=None, AxKwargs=None, Ax=None, MaxPoints=10000,
                 AvgSpectKwargs=None):

        self.MaxPoints = MaxPoints
        self.specKwargs = self.DefspecKwargs.copy()
        self.AxKwargs = self.DefAxKwargs.copy()
        self.imKwargs = self.DefImKwargs.copy()
        self.AvgSpectKwargs = self.DefAvgSpectKwargs.copy()

        if specKwargs is not None:
            self.specKwargs.update(specKwargs)

        self.Signal = Signal
        self.name = self.Signal.name

        self.units = Units
        self.Position = Position

        self.Ax = Ax
        if AxKwargs is not None:
            self.AxKwargs.update(AxKwargs)
        if imKwargs is not None:
            self.imKwargs.update(imKwargs)
        if AvgSpectKwargs is not None:
            self.AvgSpectKwargs.update(AvgSpectKwargs)
        if self.Ax is not None:
            UpdateTreeDictProp(self.Ax, self.AxKwargs)

    def CheckTime(self, Time):
        if Time is None:
            return (self.Signal.t_start, self.Signal.t_stop)

        if len(Time) == 1:
            Time = (Time[0], Time[0] + self.Signal.sampling_period)

        if Time[0] is None or Time[0] < self.Signal.t_start:
            Tstart = self.Signal.t_start
        else:
            Tstart = Time[0]

        if Time[1] is None or Time[1] > self.Signal.t_stop:
            Tstop = self.Signal.t_stop
        else:
            Tstop = Time[1]

        return (Tstart, Tstop)

    def GetSignal(self, Time, Units=None):
        if Units is None:
            _Units = self.units
        else:
            _Units = Units

        Time = self.CheckTime(Time)
        sig = self.Signal.time_slice(Time[0], Time[1])

        if _Units is not None:
            sig = sig.rescale(_Units)
        self.units = sig.units
        return sig

    def PlotSignal(self, Time, Units=None):
        sig = self.GetSignal(Time, Units)

        if 'spec' in sig.annotations:
            spec = sig
        else:
            spec = Spro.Spectrogram(sig, **self.specKwargs)

        f = spec.annotations['Freq']
        data = Spro.Resample(spec, MaxPoints=self.MaxPoints)
        img = self.Ax.imshow(np.array(data).astype(np.float).transpose(),
                             origin='lower',
                             aspect='auto',
                             extent=(spec.t_start, spec.t_stop,
                                     np.min(f), np.max(f)),
                             **self.imKwargs,
                             )
        self.img = img

    def CalcAvarage(self, TimeAvg, TimesEvent, Units=None, **Kwargs):

        sig = self.GetSignal((None, None), Units)

        spect = Spro.AvgSpectrogram(sig, 
                                    TimesEvent=TimesEvent,
                                    TimeAvg=TimeAvg,
                                    **self.AvgSpectKwargs)

        f = spect.annotations['Freq']
        img = self.Ax.imshow(np.array(spect).transpose(),
                             origin='lower',
                             aspect='auto',
                             extent=(spect.t_start, spect.t_stop,
                                     np.min(f), np.max(f)),
                             **self.imKwargs,
                             )
        self.img = img
        return spect


class SpikeSlot():
    DefLineKwargs = {'color': 'r',
                     'linestyle': '-.',
                     'alpha': 0.5,
                     'linewidth': 0.5,
                     }
    DefAxKwargs = {}

    def UpdateLineKwargs(self, LineKwargs):
        self.LineKwargs.update(LineKwargs)
        UpdateTreeDictProp(self.Line, self.LineKwargs)

    def UpdateAxKwargs(self, AxKwargs):
        self.AxKwargs.update(AxKwargs)
        UpdateTreeDictProp(self.Ax, self.AxKwargs)

    def __init__(self, Signal, Units='s',
                 Position=None, Ax=None, AxKwargs=None,
                 **LineKwargs):

        self.LineKwargs = self.DefLineKwargs.copy()
        self.AxKwargs = self.DefAxKwargs.copy()
        self.Signal = Signal
#        self.Signal.__class__ = NeoTrain
        self.name = self.Signal.name
        self.Position = Position
        self.Ax = Ax
        self.units = Units

        if AxKwargs is not None:
            self.AxKwargs.update(AxKwargs)

        if self.Ax is not None:
            UpdateTreeDictProp(self.Ax, self.AxKwargs)
        self.LineKwargs.update(LineKwargs)

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

        xmin, xmax, ymin, ymax = self.Ax.axis()

        self.Lines = self.Ax.vlines(sig,
                                    ymin, ymax,
                                    **self.LineKwargs
                                    )


class WaveSlot():

    DefTrialLineKwargs = {'color': 'k',
                          'linestyle': '-',
                          'alpha': 0.05,
                          'linewidth': 0.5,
                          'clip_on': True,
                          }

    DefLineKwargs = {'color': 'k',
                     'linestyle': '-',
                     'alpha': 1,
                     'linewidth': 0.5,
                     'clip_on': True,
                     }
    DefAxKwargs = {}

    def UpdateLineKwargs(self, LineKwargs):
        self.LineKwargs.update(LineKwargs)
        UpdateTreeDictProp(self.Line, self.LineKwargs)

    def UpdateAxKwargs(self, AxKwargs):
        self.AxKwargs.update(AxKwargs)
        UpdateTreeDictProp(self.Ax, self.AxKwargs)

    def __init__(self, Signal, Units=None, UnitsInLabel=False,
                 Position=None, Ax=None, AxKwargs=None, TrialProcessChain=None,
                 **LineKwargs):

        self.TrialLineKwargs = self.DefTrialLineKwargs.copy()
        self.LineKwargs = self.DefLineKwargs.copy()
        self.AxKwargs = self.DefAxKwargs.copy()
        
        self.TrialProcessChain = TrialProcessChain

        self.Signal = Signal
        self.name = self.Signal.name

        self.units = Units

        self.Position = Position
        self.UnitsInLabel = UnitsInLabel

        self.Ax = Ax
        if AxKwargs is not None:
            self.AxKwargs.update(AxKwargs)

        if self.Ax is not None:
            UpdateTreeDictProp(self.Ax, self.AxKwargs)

        self.LineKwargs.update(LineKwargs)
        if 'label' not in self.LineKwargs:
            self.LineKwargs.update({'label': self.name})
        else:
            self.name = self.LineKwargs['label']

        self.TrialLineKwargs['color'] = self.LineKwargs['color']
        
        self.current_time = None

    def CheckTime(self, Time):
        if Time is None:
            return (self.Signal.t_start, self.Signal.t_stop)

        if len(Time) == 1:
            Time = (Time[0], Time[0] + self.Signal.sampling_period)

        if Time[0] is None or Time[0] < self.Signal.t_start:
            Tstart = self.Signal.t_start
        else:
            Tstart = Time[0]

        if Time[1] is None or Time[1] > self.Signal.t_stop:
            Tstop = self.Signal.t_stop
        else:
            Tstop = Time[1]

        return (Tstart, Tstop)

    def GetSignal(self, Time, Units=None):
        if Units is None:
            _Units = self.units
        else:
            _Units = Units
        Time = self.CheckTime(Time)
        sig = self.Signal.time_slice(Time[0], Time[1])
        if _Units is not None:
            sig = sig.rescale(_Units)
        self.units = sig.units
        return sig

    def PlotSignal(self, Time, Units=None):
        if self.Ax is None:
            self.Fig, self.Ax = plt.subplots()

        sig = self.GetSignal(Time, Units)

        if self.UnitsInLabel is True:
            su = str(sig.units).split(' ')[-1]
            label = "{} [{}]".format(self.name, su)
        else:
            label = self.name
        self.LineKwargs.update({'label': label})

        self._PlotSignal(sig)

    def _PlotSignal(self, sig):
        self.Lines = self.Ax.plot(sig.times.rescale('s'),
                                  sig,
                                  **self.LineKwargs
                                  )
        self.current_time = (sig.t_start.rescale('s'),
                             sig.t_stop.rescale('s'))
        self.Ax.set_xlim(left=self.current_time[0],
                         right=self.current_time[1])
        self.Line = self.Lines[0]

    def CalcAvarage(self, TimeAvg, TimesEvent, Units=None,
                    PlotMean=True, PlotStd=False, PlotTrials=False,
                    TrialLineKwargs=None, StdAlpha=0.2, TrialProcessChain=None,
                    **kwargs):

        if TrialLineKwargs is not None:
            self.TrialLineKwargs.update(TrialLineKwargs)
            
        if TrialProcessChain is not None:
            self.TrialProcessChain = TrialProcessChain
        
        sig = self.GetSignal((None, None), Units)
        avg = Spro.TrigAveraging(sig,
                                 TimesEvent=TimesEvent,
                                 TimeAvg=TimeAvg,
                                 TrialProcessChain=self.TrialProcessChain)

        if PlotMean:
            self._PlotSignal(avg)

        if PlotTrials:
            acc = avg.annotations['acc']
            self.Ax.plot(acc.times,
                         acc,
                         **self.TrialLineKwargs)

        if PlotStd:
            std = avg.annotations['std']
            self.Ax.fill_between(std.times,
                                 np.array(avg+std).flatten(),
                                 np.array(avg-std).flatten(),
                                 alpha=StdAlpha,
                                 facecolor=self.LineKwargs['color'],
                                 edgecolor=None,
                                 clip_on=False)
        return avg


class ControlFigure():

    def __init__(self, pltSL, figsize=(20*0.394, 5*0.394)):

        self.pltSL = pltSL

        self.MapSlots = []
        for sl in pltSL.Slots:
            if hasattr(sl, 'Map'):
                self.MapSlots.append(sl)
        
        TMax = []
        TMin = []
        for sl in pltSL.Slots:
            if not hasattr(sl, 'Map'):
                TMax.append(sl.Signal.t_stop.rescale('s'))
                TMin.append(sl.Signal.t_start.rescale('s'))
                
        TMax = np.max(TMax)
        TMin = np.min(TMin)

        self.Fig, ax = plt.subplots(10, 1, figsize=figsize)
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

        self.TextStart = TextBox(ax[2],
                                 'Start time [s]',
                                 initial='0')
        self.TextStart.on_submit(self.submit_start)

        self.TextStop = TextBox(ax[3],
                                 'Stop time [s]',
                                 initial='10')
        self.TextStop.on_submit(self.submit_stop)
        
        self.Refresh = True
        self.OldStart = 0
        self.OldStop = 0
        
        # self.bStart = Button(ax[4],
        #                      label='Start')        
        # self.bStart.on_clicked(self.StartAnimation)
        # self.Timer = None
        # self.TextInterval = TextBox(ax[5],
        #                             'Interval [ms]',
        #                             initial='2000')

        self.bStartMapAni = Button(ax[6],
                                   label='Annimate Maps START')    
        self.TextMapAniSpeed = TextBox(ax[7],
                                       'Map Speed [1x]',
                                       initial='1')   
        self.TextMapAniPoints = TextBox(ax[8],
                                        'Map Points [n]',
                                        initial='500') 
        self.bStartMapAni.on_clicked(self.StartMapAnimation)
        self.TimerMap = None

        self.bStartSetZero = Button(ax[9],
                                   label='Set Zero at Start Time')
        self.bStartSetZero.on_clicked(self.BtSetZero)
    
    def BtSetZero(self, val):
        Twind = (self.sTstart.val * pq.s,
                 (self.sTstart.val+5) * pq.s)        
           
        for sl in self.pltSL.Slots:
            if 'LiveZero' in sl.Signal.annotations:
                if sl.Signal.annotations['LiveZero']:
                    # print(sl.Signal.name)
                    sl.Signal = Spro.SetZero(sl.Signal, TWind=Twind)
        self.Update(None)

    def StartMapAnimation(self, val):
        if self.TimerMap is not None:
            self.bStartMapAni.label.set_label('Annimate Maps START')
            self.TimerMap.stop()
            self.TimerMap = None
            return

        twind = (self.sTstart.val * pq.s,
                 self.sTstart.val * pq.s + self.sTshow.val * pq.s)

        points = int(self.TextMapAniPoints.text)
        speed = float(self.TextMapAniSpeed.text)
        interval = ((twind[1] - twind[0])/points)/speed
        self.MapCount = 0
        self.MapTimes = np.linspace(twind[0], twind[1], points)
        print('Start', interval)
        
        self.TimerMap = self.Fig.canvas.new_timer(interval=interval.rescale('ms'))
        self.TimerMap.add_callback(self.UpdateMapAnimation)
        self.TimerMap.start()
        self.bStartMapAni.label.set_label('Annimate Maps STOP')

    def UpdateMapAnimation(self):
        for sl in self.MapSlots:
            sl.PlotSignal((self.MapTimes[self.MapCount],
                           self.MapTimes[self.MapCount]+1*pq.s))

        if self.MapCount >= (len(self.MapTimes)-1):
            self.MapCount = 0
        else:
            self.MapCount += 1

    def StartAnimation(self, val):         
        if self.Timer is not None:
            self.bStart.label.set_label('Start')
            self.Timer.stop()
            self.Timer = None
            return
            
        try:
            interval = float(self.TextInterval.text)
        except:
            return

        self.Timer = self.Fig.canvas.new_timer(interval=interval)
        self.Timer.add_callback(self.UpdateAnimation)
        self.Timer.start()
        self.bStart.label.set_label('Stop')

    def UpdateAnimation(self):        
        self.sTstart.set_val(self.sTstart.val + self.sTshow.val/2)        

    def Update(self, val):
        twind = (self.sTstart.val * pq.s,
                 self.sTstart.val * pq.s + self.sTshow.val * pq.s)

        if self.Refresh:
            self.UpdateGraph(twind)

    def UpdateGraph(self, twind):
        self.pltSL.PlotChannels(twind)
        self.pltSL.Fig.canvas.draw()

    def submit_start(self, text):
        try:
            val = float(text)
        except:
            print('bad number')
            return

        if self.OldStart == val:
            return

        if val > self.OldStop:
            self.TextStart.set_val(str(self.sTstart.val))
            return
    
        self.OldStart = val
        self.sTstart.set_val(float(text))

    def submit_stop(self, text):
        try:
            val = float(text)
        except:
            print('bad number')
            return

        if self.OldStop == val:
            return
        
        if val < self.sTstart.val:
            self.TextStop.set_val(str(self.sTstart.val + self.sTshow.val))
            return
        
        self.OldStop = val
    
        show = val - self.sTstart.val 
        self.sTshow.set_val(show)

    def SetTimes(self, twind):
        self.Refresh = False
        self.TextStop.set_val(str(np.array(twind[1])))
        self.TextStart.set_val(str(np.array(twind[0])))
        self.TextStop.set_val(str(np.array(twind[1])))
        self.Refresh = True


class PlotSlots():
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

    RcGeneralParams = {
#                       'axes.spines.left': False,
#                       'axes.spines.bottom': False,
#                       'axes.spines.top': False,
#                       'axes.spines.right': False,
                       }

    FigKwargs = {}

    gridspec_Kwargs = {'width_ratios': (15, 1)}

    TimeAxisProp = {'xaxis': {'visible': True, },
                    'xlabel': 'Time [s]',
                    }

    LegendKwargs = {'fontsize': 'xx-small',
                    'ncol': 5,
                    'loc': 'upper right',
                    'frameon': False}

    def UpdateFigKwargs(self, FigKwargs):
        self.FigKwargs.update(FigKwargs)
        UpdateTreeDictProp(self.Fig, self.FigKwargs)

    def _GenerateFigure(self):

        Pos = []
        for isl, sl in enumerate(self.Slots):
            if sl.Position is None:
                sl.Position = isl
            Pos.append(sl.Position)

        self.Fig, A = plt.subplots(max(Pos) + 1, 2,
                                   sharex=True,
                                   gridspec_kw=self.gridspec_Kwargs,
                                   )

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
            UpdateTreeDictProp(sl.Ax, sl.AxKwargs)

    def __init__(self, Slots, Fig=None, FigKwargs=None, RcGeneralParams=None,
                 AxKwargs=None, TimeAxis=-1,
                 ScaleBarAx=None, LiveControl=False):

        if RcGeneralParams is not None:
            self.RcGeneralParams.update(RcGeneralParams)
        plt.rcParams.update(self.RcGeneralParams)

        if FigKwargs is not None:
            self.FigKwargs.update(FigKwargs)

        self.Slots = Slots

        self.ScaleBarAx = ScaleBarAx

        if LiveControl:
            self.CtrFig = ControlFigure(self)
        else:
            self.CtrFig = None

        if Fig is None:
            self._GenerateFigure()
        else:
            self.Fig = Fig
            self.Axs = []
            for sl in self.Slots:
                self.Axs.append(sl.Ax)

        for sl in self.Slots:
            if AxKwargs is not None:
                sl.UpdateAxKwargs(AxKwargs)

        self.TimeAxis = TimeAxis
        if self.TimeAxis is not None:
            sl = self.Slots[TimeAxis]
            sl.UpdateAxKwargs(self.TimeAxisProp)

        UpdateTreeDictProp(self.Fig, self.FigKwargs)
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

        if self.ScaleBarAx is not None:
            if self.ScaleBarKwargs['yunit'] is None:
                sl = self.SlotsInAxs[self.Axs[self.ScaleBarAx]][0]
                su = str(sl.units).split(' ')[-1]
                self.ScaleBarKwargs['yunit'] = su
            DrawBarScale(self.Axs[self.ScaleBarAx], **self.ScaleBarKwargs)

    def AddLegend(self, **LegendKwargs):
        self.LegendKwargs.update(LegendKwargs)
        for Ax in self.Axs:
            Ax.legend(**self.LegendKwargs)

    def PlotChannels(self, Time, Units=None, FormatFigure=True):
        self.ClearAxes()
        print('plot channels')
        SpectColBars.ImgDicts = {}
        for isl, sl in enumerate(self.Slots):
            sl.PlotSignal(Time, Units=Units)
            if hasattr(sl, 'img'):
                SpectColBars.ImgDicts.update({isl: sl.img})
#        if Time is not None:
#            if Time[0] is not None:
#                sl.Ax.set_xlim(left=Time[0].magnitude)
#            if Time[1] is not None:
#                sl.Ax.set_xlim(right=Time[1].magnitude)
            if sl.current_time is not None:
                if not hasattr(sl, 'Map'):
                    self.current_time = sl.current_time

        if self.CtrFig is not None:            
            self.CtrFig.SetTimes(self.current_time)

    def PlotEvents(self, Times, Labels=None, lAx=0, fontsize='xx-small',
                   LabPosition='top', duration=None, **kwargs):

        xlim = self.Axs[0].get_xlim()

        Times = Times.rescale('s')

        self.Texts = []

        if Labels is not None:
            for ilbl, lbl in enumerate(Labels):
                for ax in self.Axs:
                    ylim = ax.get_ylim()
                    if duration is not None:
                        ax.vspan(Times[ilbl],Times[ilbl]+duration,ylim[0], ylim[1], **kwargs )                        
                    else:
                        ax.vlines(Times[ilbl], ylim[0], ylim[1], **kwargs)
                lax = self.Axs[lAx]
                if LabPosition == 'top':
                    ylim = lax.get_ylim()[1]
                else:
                    ylim = lax.get_ylim()[0]
                txt = lax.text(Times[ilbl], ylim, lbl, fontsize=fontsize)
                self.Texts.append(txt)
            return

        # EventLines = []
        for ax in self.Axs:
            ylim = ax.get_ylim()
            if duration is not None:
                lines= ax.axvspan(Times,Times+duration,ylim[0], ylim[1], **kwargs )                        
            else:
                lines = ax.vlines(Times, ylim[0], ylim[1], **kwargs)
#            EventLines.append(lines[0])
        # return EventLines
        self.Axs[0].set_xlim(xlim)

    def PlotEventAvarage(self, TimeAvg, TimesEvent, Units=None, ClearAxes=True,
                         **Avgkwargs):

        if ClearAxes:
            self.ClearAxes()

        MeanSigs = []
        for isl, sl in enumerate(self.Slots):
            print('Calculating Avg ', sl.name, isl)
            MeanSig = sl.CalcAvarage(TimeAvg, TimesEvent, Units=Units,
                                     **Avgkwargs)
            MeanSigs.append(MeanSig)
            if hasattr(sl, 'img'):
                SpectColBars.ImgDicts.update({isl: sl.img})

        sl.Ax.set_xlim(left=TimeAvg[0].magnitude)
        sl.Ax.set_xlim(right=TimeAvg[1].magnitude)

        self.FormatFigure()
        return MeanSigs
    
    
class ImgSlot():
    DefAxKwargs = {}
    DefImKwargs = {
                    'vmin': -10,
                    'vmax': 10,
                    'cmap': 'seismic',
                    'interpolation': 'bicubic',
                    }
    
    def __init__(self, Signal, Ax=None,
                 AxKwargs=None, Units=None, imKwargs=None,):
        
        self.Signal = Signal
        self.Ax = Ax
        
        self.Map=True

        self.current_time = None
        self.units = Units
        
        self.imKwargs = self.DefImKwargs.copy()
        if imKwargs is not None:
            self.imKwargs.update(imKwargs)

        self.Img = self.Ax.imshow(np.array(self.Signal[0,:,:]),
                                  **self.imKwargs,
                                  )

    def CheckTime(self, Time):
        if Time is None:
            return (self.Signal.t_start,
                    self.Signal.t_start + self.Signal.sampling_period)

        if len(Time) == 1:
            return (Time[0], Time[0] + self.Signal.sampling_period)

        if Time[0] is None or Time[0] < self.Signal.t_start:
            Tstart = self.Signal.t_start
        else:
            Tstart = Time[0]

        if Time[0] > self.Signal.t_stop:
            Tstart = self.Signal.t_stop - 2*self.Signal.sampling_period

        return (Tstart, Tstart+ self.Signal.sampling_period)

    def GetSignal(self, Time, Units=None):
        if Units is None:
            _Units = self.units
        else:
            _Units = Units
        Time = self.CheckTime(Time)
        sig = self.Signal.time_slice(Time[0], Time[1])
        if _Units is not None:
            sig = sig.rescale(_Units)
        self.units = sig.units
        return sig

        
    def PlotSignal(self, Time, Units=None):
        sig = self.GetSignal(Time, Units)
    
        self.Img.set_array(np.array(sig[0, :, :]))        
        self.Ax.set_title(str(sig.t_start))
        self.current_time = (sig.t_start.rescale('s'),
                             sig.t_stop.rescale('s'))
        self.Ax.figure.canvas.draw()

    def UpdateAxKwargs(self, AxKwargs):
        pass
  
    
    