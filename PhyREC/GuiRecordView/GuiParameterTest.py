#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 11:27:48 2019

@author: aguimera
"""


import os
from PyQt5 import Qt
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


channelparams = [
                 {'name': 'ChName',
                  'type': 'str',
                  'value': 'Chxx' },
                 {'name': 'Color',
                  'type': 'color',
                  'value': "000" },
                 {'name': 'Color',
                  'type': 'color',
                  'value': "000" },
                ]


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


GenFsPar = {'name': 'Sampling Rate',
            'type': 'float',
            'value': 1e4,
            'step': 100,
            'siPrefix': True,
            'suffix': 'Hz'}

GenIntTimePar = {'name': 'Interval time',
                 'type': 'float',
                 'value': 1,
                 'step': 0.1,
                 'siPrefix': True,
                 'suffix': 's'}

GenIntSamplesPar = {'name': 'Interval samples',
                    'type': 'int',
                    'value': 1e4,
                    }

GenNChannelsPar = {'name': 'Channels Number',
                   'type': 'int',
                   'value': 16,
                   'limits': (1, 128),
                   'step': 1}


class DataGeneratorParameters(pTypes.GroupParameter):
    def __init__(self, **kwargs):
        kwargs['type'] = 'bool'
        kwargs['value'] = True
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
        super().__init__()
        self.Fs = float(Fs)
        self.nChannels = int(nChannels)
        self.nSamples = int(nSamples)
        self.OutData = np.ndarray((self.nSamples, self.nChannels))
        self.IntervalTime = IntervalTime

    def run(self, *args, **kwargs):
        while True:
            Qt.QThread.msleep(self.IntervalTime)
            self.OutData = np.random.sample(self.OutData.shape)
            self.NewSample.emit()



class MainWindow(Qt.QWidget):
    ''' Main Window '''

    def __init__(self):
        super().__init__()

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

        self.setGeometry(550, 65, 300, 200)
        self.setWindowTitle('MainWindow')

        self.btnGen.clicked.connect(self.on_btnGen)
        self.threadGen = None

    def on_pars_changed(self, param, changes):
#        print("tree changes:")
#        for param, change, data in changes:
#            path = self.pars.childPath(param)
#            if path is not None:
#                childName = '.'.join(path)
#            else:
#                childName = param.name()
#        print('  parameter: %s'% childName)
#        print('  change:    %s'% change)
#        print('  data:      %s'% str(data))
#        print('  ----------')
#        
        print('Fs         --> ', self.DataGenConf.Fs.value())
        print('IntTime    --> ', self.DataGenConf.IntTime.value())
        print('IntSamples --> ', self.DataGenConf.IntSamples.value())
        print('NChannels  --> ', self.DataGenConf.NChannels.value())

    def on_btnGen(self):
        if self.threadGen is None:
            Fs = self.DataGenConf.Fs.value()
            IntTime = self.DataGenConf.IntTime.value()
            nSamples = self.DataGenConf.IntSamples.value()
            nChannels = self.DataGenConf.NChannels.value()           
            
            self.threadGen = DataSamplingThread(Fs=Fs,
                                                nChannels=nChannels,
                                                nSamples=nSamples,
                                                IntervalTime=IntTime)

            self.threadGen.NewSample.connect(self.on_NewSample)
            self.threadGen.start()

            self.btnGen.setText("Stop Gen")
            self.OldTime = time.time()
        else:
            self.threadGen.terminate()
            self.threadGen = None
            self.btnGen.setText("Start Gen")
#    
#        self.pars2 = Parameter.create(name='params',
#                                     type='group',
#                                     children=params)
#        self.treepar.addParameters(self.pars2, showTop=False)
        
    def on_NewSample(self):
        ''' Visualization of streaming data-WorkThread. '''
        Ts = time.time() - self.OldTime
        print('Sample time', Ts, 1/Ts, (1/(Ts/self.threadGen.nSamples)), self.threadGen.OutData.shape)
#        if self.threadSave.NewData is not None:
#            print('Error New data not clear')



if __name__ == '__main__':
    app = Qt.QApplication([])
    mw  = MainWindow()
    mw.show()
    app.exec() 