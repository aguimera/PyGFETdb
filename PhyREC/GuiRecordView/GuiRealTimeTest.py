#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 11:44:05 2019

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

class DataSamplingThread(Qt.QThread):
    ''' Data generator '''
    NewSample = Qt.pyqtSignal()

    def __init__(self, Fs, nChannels, nSamples):
        super().__init__()
        self.Fs = float(Fs)
        self.Ts = 1/self.Fs
        self.nChannels = nChannels
        self.nSamples = nSamples
        self.OutData = np.ndarray((nSamples, nChannels))
        self.IntervalTime = (self.Ts * self.nSamples)*1000

    def run(self, *args, **kwargs):
        while True:
            Qt.QThread.msleep(self.IntervalTime)
            self.OutData = np.random.sample(self.OutData.shape)
            self.NewSample.emit()


class FileBuffer():
    def __init__(self, FileName, nChannels):
        os.remove(FileName)
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
    def __init__(self):
        super().__init__()
        self.NewData = None
        self.FileBuff = FileBuffer('test.h5',
                                   64)

    def run(self, *args, **kwargs):
        while True:
            if self.NewData is not None:
                self.FileBuff.AddSample(self.NewData)
                self.NewData = None
                print('Saved')
            else:
                Qt.QThread.msleep(100)

    def AddData(self, NewData):
        self.NewData = NewData


class PlottingThread(Qt.QThread):
    def __init__(self):
        super().__init__()
        self.NewData = None
        self.win = pg.GraphicsWindow()
        self.Plots = []
        self.Curves = []
        for i in range(64):
            p = self.win.addPlot()
            self.Plots.append(p)
            self.Curves.append(p.plot())
            
#        self.Fig, self.Ax = plt.subplots()

    def run(self, *args, **kwargs):        
        while True:                              
            if self.NewData is not None:                
                for i in range(64):
                    self.Curves[i].setData(self.NewData[:,i])
#                self.Ax.clear()
#                self.Ax.plot(self.NewData)
#                self.Fig.canvas.draw()
                self.NewData = None
                print('Plot')
            else:
                Qt.QThread.msleep(100)

    def AddData(self, NewData):
        self.NewData = NewData


class MainWindow(Qt.QWidget):    
    ''' Main Window '''

    def __init__(self):
        super().__init__()

        layout = Qt.QVBoxLayout(self) 
        self.labelMain = Qt.QLabel("The result of the Main task: ")
        layout.addWidget(self.labelMain)
        self.labelThread = Qt.QLabel("The result of the Thread task: ")
        layout.addWidget(self.labelThread)

        self.btnGen = Qt.QPushButton("Start Gen!")
        layout.addWidget(self.btnGen)
#        self.btnMain = Qt.QPushButton("Start Main!")
#        layout.addWidget(self.btnMain)
        self.setGeometry(550, 65, 300, 200)
        self.setWindowTitle('MainWindow')

        self.btnGen.clicked.connect(self.on_btnGen)
#        self.btnMain.clicked.connect(self.on_btnMain)

#        self.msg = MsgBox()  
        self.threadGen     = None
#        self.threadMain = None  

    def on_btnGen(self):
        ''' Starting or Stopping an Additional Stream-WorkThread from the main window '''

        if self.threadGen is None:
            self.threadSave = DataSavingThread()
            self.threadSave.start()
            self.threadPlot = PlottingThread()
            self.threadPlot.start()

            self.threadGen = DataSamplingThread(Fs=10e3,
                                                nChannels=64,
                                                nSamples=10000)

            self.threadGen.NewSample.connect(self.on_NewSample)
            self.threadGen.start()

            self.btnGen.setText("Stop Gen")
            self.OldTime = time.time()
        else:
            self.threadGen.terminate()
            self.threadGen = None
            self.threadSave.terminate()
            self.threadSave = None
            self.btnGen.setText("Start Gen")

    def on_NewSample(self):
        ''' Visualization of streaming data-WorkThread. '''
        Ts = time.time() - self.OldTime
        print('Sample time', Ts, 1/Ts, (1/(Ts/self.threadGen.nSamples)))
        if self.threadSave.NewData is not None:
            print('Error New data not clear')
        self.threadSave.AddData(self.threadGen.OutData)
        self.threadPlot.AddData(self.threadGen.OutData)
        self.OldTime = time.time()
        
        

#        self.msg.label.setText(str(value))
#        self.labelThread.setText("The result of the Thread task: " + str(value)) # We show also in the main window
#
#        # We restore the rendering of the stream window if it was closed. The flow is working.
#        if not self.msg.isVisible():        
#            self.msg.show()

#
#    def on_btnMain(self):
#        ''' Starting or Stopping the Main Thread-WorkThreadMain '''
#
#        cM = random.randrange(1, 100)
#        if self.threadMain is None:
#            self.threadMain = WorkThreadMain(cM)
#            self.threadMain.threadSignalMain.connect(self.on_threadSignalMain)
#            self.threadMain.start()
#            self.btnMain.setText("Stop Main")
#        else:
#            self.threadMain.terminate()         
#            self.threadMain = None
#            self.btnMain.setText("Start Main")
#
#    def on_threadSignalMain(self, value):
#        ''' Visualization of streaming data WorkThreadMain '''
#
#        self.labelMain.setText("The result of the Main task: " + str(value)) 



#class WorkThreadMain(Qt.QThread):
#    ''' Streaming Main task '''
#
#    threadSignalMain = Qt.pyqtSignal(int)
#    def __init__(self, startParm):
#        super().__init__()
#        self.startParm = startParm
#
#    def run(self, *args, **kwargs):
#        c = self.startParm                            
#        while True:
#            Qt.QThread.msleep(1000)
#            c += 1
#            self.threadSignalMain.emit(c)
#
#
#class MsgBox(Qt.QDialog):
#    """ Window initialization class for visualizing an additional stream
#         and a button to close the stream window if the thread is stopped! """
#
#    def __init__(self):
#        super().__init__()
#
#        layout     = Qt.QVBoxLayout(self)
#        self.label = Qt.QLabel("")
#        layout.addWidget(self.label)
#
#        close_btn  = Qt.QPushButton("Close thread")
#        layout.addWidget(close_btn)
#
#        close_btn.clicked.connect(self.close) 
#
#        self.setGeometry(900, 65, 400, 80)
#        self.setWindowTitle('MsgBox for WorkThread')
#



if __name__ == '__main__':
    app = Qt.QApplication([])
    mw  = MainWindow()
    mw.show()
    app.exec() 