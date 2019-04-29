#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 11:27:48 2019

@author: aguimera
"""

from __future__ import print_function


import os
from PyQt5 import Qt
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtGui import QPen

import numpy as np
import time

from PhyREC.NeoInterface import NeoSegment, NeoSignal
import PhyREC.PlotWaves as Rplt
import PhyREC.SignalProcess as Rpro
import PhyREC.SignalAnalysis as Ran
import matplotlib.pyplot as plt
import quantities as pq
import h5py
import numpy as np

import pyqtgraph as pg

import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
from itertools import  cycle
import copy
from scipy.signal import welch


params = [
            {'name': 'Sampling Simulation',
             'type': 'group',
             'children': [{'name': 'Sampling Rate',
                           'type': 'float',
                           'value': 1e4,
                           'step': 100,
                           'siPrefix': True,
                           'suffix': 'Hz'},
                          {'name': 'Interval samples',
                           'type': 'float',
                           'value': 1,
                           'step': 0.1,
                           'siPrefix': True,
                           'suffix': 's'},
                          {'name': 'Channels Number',
                           'type': 'int',
                           'value': 1,
                           'limits': (1, 128),
                           'step': 1},
                          {'name': 'Refresh time',
                           'type': 'str',
                           'value': '',
                           'readonly': True},
                          ]
             },
         ]


GenFsPar = {'name': 'Fs',
            'tip': 'FsKw.Fs',
            'type': 'float',
            'value': 1e4,
            'step': 100,
            'siPrefix': True,
            'suffix': 'Hz'}

GenIntTimePar = {'name': 'IntervalTime',
                 'tip': 'FsKw.Fs',
                 'type': 'float',
                 'value': 1,
                 'step': 0.1,
                 'siPrefix': True,
                 'suffix': 's'}

GenIntSamplesPar = {'name':'nSamples',
                    'tip': 'Interval samples',
                    'type': 'int',
                    'value': 1e4,
                    }

GenNChannelsPar = {'name': 'nChannels',
                   'tip': 'Channels Number',
                   'type': 'int',
                   'value': 16,
                   'limits': (1, 128),
                   'step': 1}


class DataGeneratorParameters(pTypes.GroupParameter):
    def __init__(self, **kwargs):
#        kwargs['type'] = 'bool'
#        kwargs['value'] = True
        pTypes.GroupParameter.__init__(self, **kwargs)

        self.addChild(GenFsPar)
        self.addChild(GenIntTimePar)
        self.addChild(GenIntSamplesPar)
        self.addChild(GenNChannelsPar)
        self.Fs = self.param(GenFsPar['name'])
        self.IntTime = self.param(GenIntTimePar['name'])
        self.IntSamples = self.param(GenIntSamplesPar['name'])
        self.NChannels = self.param(GenNChannelsPar['name'])

        self.Fs.sigValueChanged.connect(self.on_Fs_Changed)
        self.IntTime.sigValueChanged.connect(self.on_Time_Changed)
        self.IntSamples.sigValueChanged.connect(self.on_Samples_Changed)

    def on_Fs_Changed(self):
        Fs = self.Fs.value()
        Ts = 1/Fs
        nSamps = self.IntSamples.value()

        nTime = nSamps * Ts
        self.IntTime.setValue(nTime, blockSignal=self.on_Time_Changed)

    def on_Time_Changed(self):
        Fs = self.Fs.value()
        Ts = 1/Fs
        nTime = self.IntTime.value()

        nSamps = int(nTime/Ts)
        self.IntSamples.setValue(nSamps, blockSignal=self.on_Samples_Changed)

    def on_Samples_Changed(self):
        Fs = self.Fs.value()
        Ts = 1/Fs
        nSamps = self.IntSamples.value()

        nTime = nSamps * Ts
        self.IntTime.setValue(nTime, blockSignal=self.on_Time_Changed)


class DataSamplingThread(Qt.QThread):
    ''' Data generator '''
    NewSample = Qt.pyqtSignal()

    def __init__(self, Fs, nChannels, nSamples, IntervalTime):
        super(DataSamplingThread, self).__init__()
        
        self.Timer = Qt.QTimer()
        self.Timer.moveToThread(self)
        self.Timer.timeout.connect(self.GenData)

        self.Fs = float(Fs)
        self.nChannels = int(nChannels)
        self.nSamples = int(nSamples)
        self.OutData = np.ndarray((self.nSamples, self.nChannels))
        self.IntervalTime = IntervalTime*1000
        self.Timer.setInterval(self.IntervalTime)
        
        Pcycle = np.round(self.Fs/1)
        Fsig = Fs/Pcycle

        Ts =  1/self.Fs
        tstop = Ts*(Pcycle)
        t = np.arange(0, tstop, Ts)

        samples = np.sin(2*np.pi*Fsig*t)
        self.InSamples = cycle(samples)
        self.chFacts = np.linspace(0, nChannels/10, nChannels)
        
        for isamp in range(self.nSamples):
            samps = self.chFacts * next(self.InSamples)
            self.OutData[isamp, :] = samps

    def run(self, *args, **kwargs):        
        self.Timer.start()
        loop = Qt.QEventLoop()
        loop.exec_()
        
#        while True:
#            Qt.QThread.msleep(self.IntervalTime)
#            self.OutData = np.random.sample(self.OutData.shape)
#            self.NewSample.emit()
#    
    def GenData(self):
#        for isamp in range(self.nSamples):
#            samps = self.chFacts * next(self.InSamples)
#            self.OutData[isamp, :] = samps
#        self.OutData = self.OutData + np.random.sample(self.OutData.shape)            
#        self.OutData = np.random.sample(self.OutData.shape)
        self.NewSample.emit()


labelStyle = {'color': '#FFF',
              'font-size': '7pt',
              'bold': True}


class Buffer2D(np.ndarray):
    def __new__(subtype, shape, dtype=float, buffer=None, offset=0,
                strides=None, order=None, info=None):
        # Create the ndarray instance of our type, given the usual
        # ndarray input arguments.  This will call the standard
        # ndarray constructor, but return an object of our type.
        # It also triggers a call to InfoArray.__array_finalize__
        obj = super(Buffer2D, subtype).__new__(subtype, shape, dtype,
                                               buffer, offset, strides,
                                               order)
        # set the new 'info' attribute to the value passed
        obj.bufferind = 0
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        self.bufferind = getattr(obj, 'bufferind', None)

    def AddData(self, NewData):
        newsize = NewData.shape[0]
        self[0:-newsize, :] = self[newsize:, :]
        self[-newsize:, :] = NewData
        self.bufferind += newsize

    def IsFilled(self):
        return self.bufferind >= self.shape[0]

    def Reset(self):
        self.bufferind = 0


class FFTThread(Qt.QThread):
    def __init__(self, nChannels, ChannelConf):
        super(FFTThread, self).__init__()

        self.NewData = None

        self.nFFT = 2**15
        self.nChannels = nChannels
        self.Fs = 10e3
        self.BufferSize = self.nFFT * 5
        self.Buffer = Buffer2D((self.BufferSize, self.nChannels))

        self.Plots = [None]*nChannels
        self.Curves = [None]*nChannels
        for wind, chs in ChannelConf.items():
            p = wind.addPlot()            
            for ch in chs:
                self.Plots[ch['Input index']] = p
                c = p.plot(pen=pg.mkPen(ch['color'],
                                        width=ch['width']))
                self.Curves[ch['Input index']] = c

    def run(self, *args, **kwargs):
        while True:
            if self.NewData is not None:
                self.Buffer.AddData(self.NewData)
                self.NewData = None
                if self.Buffer.IsFilled():
                    ff, psd = welch(self.Buffer,
                                    fs=self.Fs,
                                    nperseg=self.nFFT,
                                    axis=0)
                    self.Buffer.Reset()
                    for i in range(self.nChannels):
                        self.Curves[i].setData(ff, psd[:, i])
            else:
                Qt.QThread.msleep(10)

    def AddData(self, NewData):
        if self.NewData is not None:
            print('Error FFT !!!!')
            return
        self.NewData = NewData


class PlottingThread(Qt.QThread):

    def __init__(self, nChannels, ChannelConf):
        super(PlottingThread, self).__init__()

        self.NewData = None

        self.winds = []
        self.nChannels = nChannels
        self.Plots = [None]*nChannels
        self.Curves = [None]*nChannels

        self.Fs = 10e3
        self.Ts = 1/float(self.Fs)
        self.ViewTime = 10
        self.BufferSize = int(self.ViewTime/self.Ts)
        self.Buffer = Buffer2D((self.BufferSize, self.nChannels))

        self.Winds = []
        for win, chs in ChannelConf.items():
            wind = PlotWindow()
            self.Winds.append(wind)
            xlink = None
            for ch in chs:
                wind.plot.nextRow()
                p = wind.plot.addPlot()
                p.hideAxis('bottom')
                p.setLabel('left',
                           ch['Name'],
                           units='A',
                           **labelStyle)
                c = p.plot(pen=pg.mkPen(ch['color'],
                                        width=ch['width']),
                           autoDownsample=True)
#                c = p.plot()
                self.Plots[ch['Input index']] = p
                self.Curves[ch['Input index']] = c

                if xlink is not None:
                    p.setXLink(xlink)
                xlink = p
            p.showAxis('bottom')
            p.setLabel('bottom', 'Time', units='s', **labelStyle)

    def run(self, *args, **kwargs):
        while True:
            if self.NewData is not None:
                self.Buffer.AddData(self.NewData)
                self.NewData = None
                for i in range(self.nChannels):
                    self.Curves[i].setData(self.Buffer[:, i])
                    self.Plots[i].setXRange(self.BufferSize/10,
                                            self.BufferSize)
            else:
#                pg.QtGui.QApplication.processEvents()
                Qt.QThread.msleep(10)

    def AddData(self, NewData):
        if self.NewData is not None:
            print('Error plotting!!!!')
        self.NewData = NewData


class PlotWindow(Qt.QWidget):
    ''' Main Window '''

    def __init__(self):
        super(PlotWindow, self).__init__()
        
        layout = Qt.QVBoxLayout(self) 

        self.plot = pg.GraphicsLayoutWidget()
#        self.LayOut = self.plot
        
        layout.addWidget(self.plot)
        self.show()


class Plotting():

    def __init__(self, nChannels, ChannelConf):
        self.winds = []
        self.nChannels = nChannels
        self.Plots = [None]*nChannels
        self.Curves = [None]*nChannels

        self.Fs = 10e3
        self.Ts = 1/float(self.Fs)
        self.ViewTime = 10
        self.BufferSize = int(self.ViewTime/self.Ts)
        self.Buffer = Buffer2D((self.BufferSize, self.nChannels))

        self.Winds = []
        for win, chs in ChannelConf.items():
#            wind = pg.GraphicsWindow(title="Real Time Plot")
            wind = PlotWindow()
            self.Winds.append(wind)
            xlink = None
            for ch in chs:
                wind.plot.nextRow()
                p = wind.plot.addPlot()
                p.hideAxis('bottom')
                p.setLabel('left',
                           ch['Name'],
                           units='A',
                           **labelStyle)
                p.setDownsampling(mode='peak')
                p.setClipToView(True)
                c = p.plot(pen=pg.mkPen(ch['color'],
                                        width=ch['width']))
#                c = p.plot()
                self.Plots[ch['Input index']] = p
                self.Curves[ch['Input index']] = c

                if xlink is not None:
                    p.setXLink(xlink)
                xlink = p
            p.showAxis('bottom')
            p.setLabel('bottom', 'Time', units='s', **labelStyle)

    def AddData(self, NewData):
        self.Buffer.AddData(NewData)        
        for i in range(self.nChannels):
            self.Curves[i].setData(self.Buffer[:, i])
#            self.Curves[i].setData(NewData[:, i])
            self.Plots[i].setXRange(self.BufferSize/10,
                                    self.BufferSize)


class FileBuffer():
    def __init__(self, FileName, nChannels):
        self.FileName = FileName
        self.nChannels = nChannels
        self.h5File = h5py.File(FileName, 'w')
        self.Dset = self.h5File.create_dataset('data',
                                               shape=(0, nChannels),
                                               maxshape=(None, nChannels),
                                               compression="gzip")

        self.Sigs = []
        for i in range(nChannels):
            self.Sigs.append(NeoSignal(signal=self.Dset[:, i],
                                       units='V',
                                       sampling_rate=10e3*pq.Hz,
                                       t_start=0*pq.s,
                                       copy=False,
                                       name='ch{}'.format(i)))

    def AddSample(self, Sample):
        nSamples = Sample.shape[0]
        FileInd = self.Dset.shape[0]
        self.Dset.resize((FileInd + nSamples, self.nChannels))
        self.Dset[FileInd:, :] = Sample
        self.h5File.flush()


class DataSavingThread(Qt.QThread):
    def __init__(self, FileName, nChannels, MaxSize=None):
        super(DataSavingThread, self).__init__()
        self.NewData = None
        self.FileBuff = FileBuffer(FileName,
                                   nChannels)

    def run(self, *args, **kwargs):
        while True:
            if self.NewData is not None:
                self.FileBuff.AddSample(self.NewData)
                self.NewData = None
            else:
                Qt.QThread.msleep(10)

    def AddData(self, NewData):
        if self.NewData is not None:
            print('Error Saving !!!!')
        self.NewData = NewData


ChPars = {'name': 'Ch01',
          'type': 'group',
          'children': [{'name': 'Name',
                        'type': 'str',
                        'value': 'Ch10'},
                       {'name': 'color',
                        'type': 'color',
                        'value': "FFF"},
                       {'name': 'width',
                        'type': 'float',
                        'value': 0.5},
                       {'name': 'Plot Window',
                        'type': 'int',
                        'value': 1,},
                       {'name': 'Input index',
                        'type': 'int',
                        'readonly': True,
                        'value': 1,}
                       ]}



SaveFilePars = [{'name': 'Save File',
                 'type': 'action'},
                {'name': 'File Path',
                 'type': 'str',
                 'value': ''},
                {'name': 'MaxSize',
                 'type': 'int',
                 'siPrefix': True,
                 'suffix': 'B',
                 'limits': (1e6, 1e9),
                 'step': 1e6,
                 'value': 50e6}
                ]

       

class MainWindow(Qt.QWidget):
    ''' Main Window '''

    def __init__(self):
        super(MainWindow, self).__init__()

        layout = Qt.QVBoxLayout(self) 

        self.btnGen = Qt.QPushButton("Start Gen!")
        layout.addWidget(self.btnGen)

        self.DataGenConf = DataGeneratorParameters(name='Data Generator')
        self.pars = Parameter.create(name='params',
                                     type='group',
                                     children=(self.DataGenConf,))

        self.pars.sigTreeStateChanged.connect(self.on_pars_changed)
        self.treepar = ParameterTree()
        self.treepar.setParameters(self.pars, showTop=False)
        self.treepar.setWindowTitle('pyqtgraph example: Parameter Tree')

        layout.addWidget(self.treepar)

#        self.setGeometry(550, 10, 300, 700)
        self.setWindowTitle('MainWindow')

        self.btnGen.clicked.connect(self.on_btnGen)
        self.threadGen = None

        self.FileParams = Parameter.create(name='File Params',
                                           type='group',
                                           children=SaveFilePars)
        self.pars.addChild(self.FileParams)
        self.FileParams.param('Save File').sigActivated.connect(self.FileDialog)

        self.GenChannelsViewParams(nChannels=self.DataGenConf.NChannels.value(),
                                   nWindows=1)

    def FileDialog(self):
        RecordFile, _ = QFileDialog.getSaveFileName(self,
                                                    "Recording File",
                                                    "",
                                                    )
        if RecordFile:
            if not RecordFile.endswith('.h5'):
                RecordFile = RecordFile + '.h5'
            self.FileParams.param('File Path').setValue(RecordFile)
            
    def on_pars_changed(self, param, changes):
        print("tree changes:")
        for param, change, data in changes:
            path = self.pars.childPath(param)
            if path is not None:
                childName = '.'.join(path)
            else:
                childName = param.name()
        print('  parameter: %s'% childName)
        print('  change:    %s'% change)
        print('  data:      %s'% str(data))
        print('  ----------')
       
        if childName == 'Data Generator.nChannels':
            self.pars.removeChild(self.Parch)
            self.GenChannelsViewParams(nChannels=data,
                                       nWindows=2)
        
        if childName == 'Channels View.Plot Windows':
            self.pars.removeChild(self.Parch)
            self.GenChannelsViewParams(nChannels=self.DataGenConf.NChannels.value(),
                                       nWindows=data)
            
            #        print(param, changes)
#        if changes[0].name() == 'nChannels':
#            print('New nCh', changes[1].name())
#        
#        print('Fs         --> ', self.DataGenConf.Fs.value())
#        print('IntTime    --> ', self.DataGenConf.IntTime.value())
#        print('IntSamples --> ', self.DataGenConf.IntSamples.value())
#        print('NChannels  --> ', self.DataGenConf.NChannels.value())

    def GenChannelsViewParams(self, nChannels, nWindows):
        ChParams = [{'name': 'Plot Windows',
                     'type': 'int',
                     'value': nWindows},]

        chPWind = int(nChannels/nWindows)                
        for i in range(nChannels):
            Ch = copy.deepcopy(ChPars)
            pen = pg.mkPen((i,1.3*nChannels))
            chn = 'Ch{0:02}'.format(i)
            Ch['name'] = chn
            Ch['children'][0]['value'] = chn
            Ch['children'][1]['value'] = pen.color()
            Ch['children'][3]['value'] = int(i/chPWind)
            Ch['children'][4]['value'] = i
            ChParams.append(Ch)
        
        self.Parch = Parameter.create(name='Channels View',
                                     type='group',
                                     children=ChParams)        
        self.pars.addChild(self.Parch)

    
    def GetChannelsPars(self):
        channelspars = {}
        for i in range(self.Parch.param('Plot Windows').value()):
            channelspars[i] = []
        
        for p in self.pars.child('Channels View').children():
            if not p.hasChildren():
                continue
            chp = {}
            for pp in p.children():
                chp[pp.name()] = pp.value()
            channelspars[chp['Plot Window']].append(chp.copy())
        return channelspars

    def on_btnGen(self):
        self.GetChannelsPars()
        
        GenKwargs = {}
        for p in self.pars.child('Data Generator').children():
            GenKwargs[p.name()] = p.value()
        print(GenKwargs)    

        if self.threadGen is None:

            self.threadGen = DataSamplingThread(**GenKwargs)

            self.threadGen.NewSample.connect(self.on_NewSample)
            self.threadGen.start()

            self.btnGen.setText("Stop Gen")
            self.OldTime = time.time()
#
            plotpars = self.GetChannelsPars()
            self.Plots = Plotting(GenKwargs['nChannels'],
                                       plotpars)
#            self.Plots = PlottingThread(GenKwargs['nChannels'],
#                                       plotpars)
#            self.Plots.start()

#            self.Plot = PlotWindow()
#            self.Plot.show()
#            fftplotopts = {}
#            for iw, ch in plotpars.items():
#                wind = pg.GraphicsWindow(title="Real Time Plot")
#                fftplotopts[wind] = ch            
#            self.threadPlotFFT = FFTThread(GenKwargs['nChannels'],
#                                           fftplotopts)
#            self.threadPlotFFT.start()


            FileName = self.FileParams.param('File Path').value()
            if os.path.isfile(FileName):
                print('Remove File')
                os.remove(FileName)
            self.threadSave = DataSavingThread(FileName=FileName,
                                               nChannels=self.threadGen.nChannels)
            self.threadSave.start()

        else:
            self.threadGen.terminate()
            self.threadGen = None

#            self.threadPlot.terminate()
#            self.threadPlot = None

            self.threadSave.terminate()
            self.threadSave = None

            self.btnGen.setText("Start Gen")

    def on_NewSample(self):
        ''' Visualization of streaming data-WorkThread. '''
        Ts = time.time() - self.OldTime
        self.OldTime = time.time()
        print('Sample time', Ts)
        self.threadSave.AddData(self.threadGen.OutData)
        self.Plots.AddData(self.threadGen.OutData)
#        self.threadPlotFFT.AddData(self.threadGen.OutData)
        
        

#        if self.threadSave.NewData is not None:
#            print('Error New data not clear')



if __name__ == '__main__':
    app = Qt.QApplication([])
    mw  = MainWindow()
    mw.show()
    app.exec_() 