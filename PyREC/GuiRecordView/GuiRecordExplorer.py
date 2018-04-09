#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 18:09:07 2017

@author: aguimera
"""

import os
import neo
import numpy as np
import matplotlib.pyplot as plt
import quantities as pq
import sys
from qtpy.QtWidgets import QHeaderView
from qtpy import QtWidgets, uic
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QFileDialog, QColorDialog, QInputDialog, QAction
from PyREC.RecordCore import NeoRecord
from PyREC.RecordPlot import PltSlot, PlotRecord, threshold_detection


class RecordExplorer(QtWidgets.QMainWindow):
    OutFigFormats = ('svg', 'png')
    DevDCVals = None
    # 'SlotProperty' (Header, position, string Conv, Editable)
    ChFields = {'Position': ('Position', 0, 1, 1),
                'DispName': ('DispName', 1, 0, 1),
                'SigName': ('SigName', 2, 0, 0),
                'PlotType': ('PlotType', 3, 0, 0),
                'TStart': ('TStart', 4, 2, 1),
                'OutType': ('Out', 5, 0, 1),
                'IVGain': ('IVGain', 6, 3, 1),
                'Gm': ('Gm', 7, 3, 1),
                'GmSignal': ('GmSig', 8, 3, 1),
                'FiltType': ('FiltType', 9, 4, 1),
                'FiltOrder': ('FiltOrder', 10, 5, 1),
                'FiltF1': ('FiltF1', 11, 5, 1),
                'FiltF2': ('FiltF2', 12, 5, 1),
                'RecN': ('RecNumber', 13, 1, 0),
                'SpecFmin': ('SpecFmin', 14, 3, 1),
                'SpecFmax': ('SpecFmax', 15, 3, 1),
                'SpecTimeRes': ('SpecTimeRes', 16, 3, 1),
                'SpecMinPSD': ('SpecMinPSD', 17, 3, 1),
                'SpecCmap': ('SpecCmap', 18, 0, 1),
                'Ymax': ('Y Max', 19, 3, 1),
                'Ymin': ('Y min', 20, 3, 1),
                'FileName': ('FileName', 21, 0, 0)}

    def InitMenu(self):
        mainMenu = self.menuBar
        fileMenu = mainMenu.addMenu('File')

        SaveRec1Action = QAction('Save Rec1 as', self)
        SaveRec1Action.setStatusTip('Save a copy of Rec1 and Evens')
        SaveRec1Action.triggered.connect(self.SaveRec1)
        fileMenu.addAction(SaveRec1Action)

        SaveJoinRecAction = QAction('Join Rec1 and Rec2 as', self)
        SaveJoinRecAction.setStatusTip('Save a copy of Rec1 and Rec2 with Tstart')
        SaveJoinRecAction.triggered.connect(self.SaveJoinRec)
        fileMenu.addAction(SaveJoinRecAction)

        SaveFigAction = QAction('Save Figures', self)
        SaveFigAction.setShortcut('Ctrl+s')
        SaveFigAction.setStatusTip('Save all open figures')
        SaveFigAction.triggered.connect(self.SaveFigures)
        fileMenu.addAction(SaveFigAction)

        CloseFigsAction = QAction('Close Figures', self)
        CloseFigsAction.setStatusTip('Close all open figures')
        CloseFigsAction.triggered.connect(self.CloseFigures)
        fileMenu.addAction(CloseFigsAction)

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        uipath = os.path.join(os.path.dirname(__file__), 'GuiRecordExplorer.ui')
        uic.loadUi(uipath, self)
        
        self.InitMenu()

        self.setWindowTitle('Record Explorer')

        self.ButFileRec.clicked.connect(self.ButFileRecClick)
        self.ButADDRec.clicked.connect(self.ButADDRecClick)
        self.ButDelete.clicked.connect(self.ButDeleteClick)

        self.ButFileRec2.clicked.connect(self.ButFileRec2Click)
        self.ButADDRec2.clicked.connect(self.ButADDRec2Click)

        self.SpnTStart.valueChanged.connect(self.SpnTStartChanged)
        self.SpnTShow.valueChanged.connect(self.SpnTShowChanged)
        self.SLTStart.valueChanged.connect(self.SLTStartChanged)
        self.SLTStop.valueChanged.connect(self.SLTStopChanged)
        self.ButNextEvent.clicked.connect(self.ButNextEventClick)

        self.ButPSD.clicked.connect(self.ButButPSDClick)
        self.ButEventAvg.clicked.connect(self.ButEventAvgClick)

        self.ButSaveView.clicked.connect(self.ButSaveViewClick)

        self.ButThres.clicked.connect(self.ButThresClick)
        self.ButPSDSNR.clicked.connect(self.ButPSDSNRClick)

        self.ButLoadDC.clicked.connect(self.ButLoadDCClick)
        self.ButViewDC.clicked.connect(self.ButViewDCClick)

        self.ButColor.clicked.connect(self.ButColorClick)

        # Init Tables
        self.TblChs.setColumnCount(len(self.ChFields))
        for c in self.ChFields.values():
            self.TblChs.setHorizontalHeaderItem(c[1],
                                                QtWidgets.QTableWidgetItem(c[0]))
        header = self.TblChs.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeToContents)

        self.ButPlot.clicked.connect(self.ButPlotClick)

    def CloseFigures(self):
        plt.close('all')

    def SaveFigures(self):
        Dir = QFileDialog.getExistingDirectory(self)
        Prefix, okPressed = QInputDialog.getText(self,
                                                 'Prefix',
                                                 'Prefix for files',
                                                 text='Figure')

        if Dir and okPressed:
            for i in plt.get_fignums():
                plt.figure(i)
                for ext in self.OutFigFormats:
                    fileOut = Dir + '/' + Prefix + '{}.' + ext
                    print fileOut
                    plt.savefig(fileOut.format(i))

    def ButColorClick(self):
        color = QColorDialog.getColor()
        if color.isValid():
            print color.name(), color.rgba64()

    def ButFileRecClick(self):

        RecordFile, _ = QFileDialog.getOpenFileName(self,
                                                    "Recording File", "",
                                                    "NixIO (*.h5);;Spike2 (*.smr)")

        print RecordFile
        if RecordFile:
            Gain = self.SpnRecGain.value()
            self.rec = NeoRecord(RecordFile, UnitGain=Gain)

            self.LblRecFile.setText(RecordFile.split('/')[-1])
            self.UpdateTimeControl()

            self.LstChRec.clear()
            for ch in sorted(self.rec.SigNames.keys()):
                self.LstChRec.addItem(ch)

            self.LstEvents.clear()
            for ch in sorted(self.rec.EventNames.keys()):
                self.LstEvents.addItem(ch)

    def ButFileRec2Click(self):
        RecordFile, _ = QFileDialog.getOpenFileName(self,
                                                    "Recording File", "",
                                                    "NixIO (*.h5);;Spike2 (*.smr)")

        if RecordFile:
            Gain = self.SpnRec2Gain.value()
            self.rec2 = NeoRecord(RecordFile, UnitGain=Gain)

            self.LblRecFile2.setText(RecordFile.split('/')[-1])

            self.LstChRec2.clear()
            for ch in sorted(self.rec2.SigNames.keys()):
                self.LstChRec2.addItem(ch)

            self.LstEvents2.clear()
            for ch in sorted(self.rec2.EventNames.keys()):
                self.LstEvents2.addItem(ch)

    def ButDeleteClick(self):
        try:
            sel = self.TblChs.selectedRanges()
            self.TblChs.removeRow(sel[0].topRow())
        except:
            print 'ERROR deleting row'

    def ButADDRec2Click(self):
        self.TblChs.setSortingEnabled(False)

        chs = []
        for s in self.LstChRec2.selectedItems():
            chs.append(s.text())

        rc = len(chs)
        offset = self.TblChs.rowCount()
        self.TblChs.setRowCount(rc+offset)
        index = -1

        for ir, ch in enumerate(chs):
            index += 1
            Slinf = PltSlot()
            Slinf.TStart = self.SpnTstart2.value()*pq.s
            Slinf.RecN = 2
            Slinf.DispName = ch
            Slinf.SigName = ch
            Slinf.Position = index+offset
            Slinf.FileName = self.rec2.RecFile.filename
            Slinf.PlotType = self.CmbWaveType.currentText()
            Slinf.OutType = 'V'

            self.UpdateRownTableChannels(Slinf,
                                         Row=index+offset)

        self.TblChs.setSortingEnabled(True)

    def ButADDRecClick(self):
        self.TblChs.setSortingEnabled(False)

        chs = []
        for s in self.LstChRec.selectedItems():
            chs.append(s.text())

        rc = len(chs)
        offset = self.TblChs.rowCount()
        self.TblChs.setRowCount(rc+offset)
        index = -1

        for ir, ch in enumerate(chs):
            index += 1
            Slinf = PltSlot()
            Slinf.TStart = self.rec.GetTstart(ch)
            Slinf.RecN = 1
            Slinf.DispName = ch
            Slinf.SigName = ch
            Slinf.Position = index+offset
            Slinf.FileName = self.rec.RecFile.filename
            Slinf.PlotType = self.CmbWaveType.currentText()
            if ch.startswith('T'):
                Slinf.OutType = 'I'
            else:
                Slinf.OutType = 'V'

            self.UpdateRownTableChannels(Slinf,
                                         Row=index+offset)

        self.TblChs.setSortingEnabled(True)

    def UpdateRownTableChannels(self, Slinf, Row):

        for f, c in self.ChFields.iteritems():
            item = QtWidgets.QTableWidgetItem()
            if c[3] == 0:
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)

            val = getattr(Slinf, f)
            if c[2] == 1 or c[2] == 2 or c[2] == 3:
                item.setData(Qt.DisplayRole, str(val))
            else:
                item.setData(Qt.DisplayRole, val)

            self.TblChs.setItem(Row, c[1], item)

    def ButPlotClick(self):

        Rows = self.TblChs.rowCount()

        slots = []

        for r in range(Rows):
            sl = PltSlot()

            for f, c in self.ChFields.iteritems():
                v = self.TblChs.item(r, c[1]).text()
                if f == 'GmSignal':
                    if v != '':
                        setattr(sl, f, (self.rec2, v))
                    continue

                if c[2] == 1:
                    setattr(sl, f, int(v))
                elif c[2] == 2:
                    vv = v.split(' ')
                    setattr(sl, f, float(vv[0])*pq.s)
                elif c[2] == 3:
                    setattr(sl, f, float(v))
                elif c[2] == 4:
                    vl = v.split(',')
                    print vl
                    setattr(sl, f, vl)
                elif c[2] == 5:
                    vl = v.split(',')
                    try:
                        vf = [float(i) for i in vl]
                        setattr(sl, f, vf)
                    except:
                        print 'Error in input data ', v, vl
                else:
                    setattr(sl, f, v)

            if sl.RecN == 1:
                setattr(sl, 'rec', self.rec)
            else:
                setattr(sl, 'rec', self.rec2)

#            if not sl.FiltType[0] == '':
#                sl.DispName = sl.DispName+'F'

            slots.append(sl)

        self.pltrec = PlotRecord()
        self.pltrec.CreateFig(slots)
        plt.show()

        self.LstSlots.clear()
        for i, sl in enumerate(self.pltrec.Slots):
            self.LstSlots.addItem(sl.DispName)

    def UpdateTimeControl(self):
        starts = []
        stops = []

        for sig in self.rec.Seg.analogsignals:
            starts.append(sig.t_start)
            stops.append(sig.t_stop)

        Ctrs = (self.SLTStart,
                self.SLTStop,
                self.SpnTStart,
                self.SpnAvgTimeStart,
                self.SpnAvgTimeStop,
                self.SpnPSDTimeStart,
                self.SpnPSDTimeStop)

        print stops, starts
        self.SpnTShow.setMaximum(max(stops)-min(starts))

        self.UpdateTstart(Ctrs, min(starts))
        self.UpdateTstop(Ctrs, max(stops))
        self.LblTstartMax.setText(str(max(stops)))
        self.LblStopMax.setText(str(max(stops)))

    def UpdateTstart(self, Controls, Value):
        for ctr in Controls:
            ctr.setMinimum(Value)

    def UpdateTstop(self, Controls, Value):
        for ctr in Controls:
            ctr.setMaximum(Value)

    def SpnTStartChanged(self):
        self.SLTStart.setValue(self.SpnTStart.value())
        self.SLTStop.setValue(self.SpnTStart.value()+self.SpnTShow.value())

    def SpnTShowChanged(self):
        self.SLTStop.setValue(self.SpnTStart.value()+self.SpnTShow.value())

    def SLTStartChanged(self):
        self.SpnTStart.setValue(self.SLTStart.value())
        self.SLTStop.setValue(self.SpnTStart.value()+self.SpnTShow.value())

    def SLTStopChanged(self):
        if (self.SLTStop.value() < self.SLTStart.value()):
            self.SLTStop.setValue(self.SLTStart.value()+1)
        self.SpnTShow.setValue(self.SLTStop.value()-self.SLTStart.value())
        self.UpdatePlot()

    def UpdatePlot(self):
        Tstart = self.SLTStart.value()
        Tstop = self.SLTStop.value()

        self.SpnAvgTimeStart.setValue(Tstart)
        self.SpnAvgTimeStop.setValue(Tstop)
        self.SpnPSDTimeStart.setValue(Tstart)
        self.SpnPSDTimeStop.setValue(Tstop)

        resampPoints = self.SpnResampPoints.value()
        resamp = self.ChckResamp.isChecked()

        try:
            self.pltrec.ClearAxes()
            self.pltrec.PlotChannels((Tstart*pq.s, Tstop*pq.s),
                                     ResampPoints=resampPoints,
                                     Resamp=resamp)
            self.pltrec.Fig.canvas.draw()
        except:
            print "Unexpected error:", sys.exc_info()[0]

        (EventRec, EventName) = self.GetSelectedEvent()
        if EventName:
            try:
                self.pltrec.PlotEvents((Tstart*pq.s, Tstop*pq.s),
                                       (EventRec, EventName))
            except:
                print 'Event error'
            self.pltrec.Fig.canvas.draw()

    def ButNextEventClick(self):
        EventName = self.GetSelectedEvent()

        if EventName:
            etimes = EventName[0].GetEventTimes(EventName[1]) 
            Tstart = self.SLTStart.value()
            
            ne = np.where(etimes > Tstart)
            print etimes[ne[0][1]]          
            
            tsw = self.SpnTShow.value()
            self.SLTStart.setValue(etimes[ne[0][1]].magnitude-tsw/2)
        
    def ButButPSDClick(self):
        Tstart = self.SLTStart.value()
        Tstop = self.SLTStop.value()
        nFFT = 2**self.SpnPSDViewnFFT.value()

        self.pltrec.PlotPSD((Tstart*pq.s, Tstop*pq.s), nFFT=nFFT)
        plt.show()
        self.pltrec.FigFFT.canvas.draw()

    def ButEventAvgClick(self):

        Tstart = self.SpnAvgTimeStart.value()
        Tstop = self.SpnAvgTimeStop.value()

        Estart = self.SpnAvgEventStart.value()
        Estop = self.SpnAvgEventStop.value()

        EventName = self.GetSelectedEvent()

        if EventName:
            self.pltrec.PlotEventAvg(EventName,
                                     Time=(Tstart*pq.s, Tstop*pq.s),
                                     TimeWindow=(Estart*pq.s, Estop*pq.s),
                                     OverLap=self.ChkEventOverlap.isChecked(),
                                     Std=self.ChkEventStd.isChecked(),
                                     Spect=self.ChkEventSpect.isChecked())

    def GetSelectedEvent(self):
        EventName = None
        EventRec = None

        if self.ChckEvent.isChecked():
            events = []
            for s in self.LstEvents.selectedItems():
                events.append(s.text())
            if len(events) > 0:
                EventName = events[0]
                EventRec = self.rec

            if not EventName:
                for s in self.LstEvents2.selectedItems():
                    events.append(s.text())
                if len(events) > 0:
                    EventName = events[0]
                    EventRec = self.rec2

        return EventRec, EventName

    def SaveJoinRec(self):
        RecordFile, _ = QFileDialog.getSaveFileName(self,
                                                    "Rec1 file", "",
                                                    "NixIO (*.h5)")
        print RecordFile
        if RecordFile:
            for sig in self.rec2.Seg.analogsignals:
                print sig.name
                sig.t_start = self.SpnTstart2.value()*pq.s
                self.rec.AddSignal(sig*self.rec.UnitGain)
            self.rec.SaveRecord(RecordFile)

    def SaveRec1(self):
        RecordFile, _ = QFileDialog.getSaveFileName(self,
                                                    "Rec1 file", "",
                                                    "NixIO (*.h5)")
        if RecordFile:
            self.rec.SaveRecord(RecordFile)

    def ButSaveViewClick(self):  # TODO finish save

        if self.ChkSaveRec1.isChecked():
            RecordFile, _ = QFileDialog.getSaveFileName(self,
                                                        "Rec1 file", "",
                                                        "NixIO (*.h5)")
            if RecordFile:
                out_f = neo.io.NixIO(filename=RecordFile)
                out_f.write_block(self.rec.Block)
                out_f.close()

        if self.ChkSaveRec2.isChecked():
            RecordFile, _ = QFileDialog.getSaveFileName(self,
                                                        "Rec2 file", "",
                                                        "NixIO (*.h5)")
            if RecordFile:
                out_f = neo.io.NixIO(filename=RecordFile)
                out_f.write_block(self.rec2.Block)
                out_f.close()


        if self.ChkSaveView.isChecked():
            RecordFile, _ = QFileDialog.getSaveFileName(self,
                                                        "Current View file", "",
                                                        "NixIO (*.h5)")
            if RecordFile:
                Tstart = self.SLTStart.value()*pq.s
                Tstop = self.SLTStop.value()*pq.s

                out_seg = neo.Segment(name='NewSeg')

                for sl in self.pltrec.Slots:
                    SName = sl.DispName
                    print SName
                    sig = sl.rec.Signal(SName)
                    out_seg.analogsignals.append(sig.time_slice(Tstart, Tstop))

                out_bl = neo.Block(name='NewBlock')
                out_bl.segments.append(out_seg)
                out_f = neo.io.NixIO(filename='Test_out.h5')
                out_f.write_block(out_bl)
                out_f.close()

    def ButThresClick(self):  # TODO move to PltSlot
        
        sl = self.pltrec.Slots[self.LstSlots.currentRow()]
        sig = sl.GetSignal((sl.Signal().t_start, sl.Signal().t_stop),
                           Resamp=False)

#        sign = str(self.CmbThSign.currentText())
#        print sign, type(sign)
        eve = threshold_detection(sig,
                                  threshold=float(self.TxvThVal.text())*sig.units,
                                  sign=self.CmbThSign.currentText(),
                                  RelaxTime=float(self.TxvThRTime.text())*pq.s)

        for e in eve:
            print e

        sl.rec.AddEvent(eve, self.TxvThName.text())

        if sl.RecN == 1:
            self.LstEvents.clear()
            for ch in sorted(self.rec.EventNames.keys()):
                self.LstEvents.addItem(ch)
        elif sl.RecN == 2:
            self.LstEvents2.clear()
            for ch in sorted(self.rec2.EventNames.keys()):
                self.LstEvents2.addItem(ch)

#        self.LstEvents.clear()
#        for ch in sorted(self.rec.EventNames.keys()):
#            self.LstEvents.addItem(ch)

    def ButLoadDCClick(self):
        DCFile, _ = QFileDialog.getOpenFileName(self,
                                                "DC charact file", "",
                                                "(*.h5);; (*.*)")

        if DCFile:
            self.DevDCVals, _ = FETData.LoadDataFromFile(DCFile)
            self.LblDCFile.setText(DCFile.split('/')[-1])
            self.ButViewDCClick()

    def ButViewDCClick(self):
        if self.DevDCVals:
            pltDC = FETplt.PyFETPlot()
            pltDC.AddAxes(('Ids', 'Gm', 'Rds'))
            pltDC.PlotDataCh(self.DevDCVals)
            pltDC.AddLegend()
            plt.show()
            pltDC.Fig.canvas.draw()
        else:
            self.ButLoadDCClick()

    def ButPSDSNRClick(self):
        if not self.DevDCVals:
            self.ButLoadDCClick()
            if not self.DevDCVals:
                return

        Tstart = self.SpnPSDTimeStart.value()*pq.s
        Tstop = self.SpnPSDTimeStop.value()*pq.s
        Estart = self.SpnPSDEventStart.value()*pq.s
        Estop = self.SpnPSDEventStop.value()*pq.s
        nFFT = 2**self.SpnPSDnFFT.value()

        EventName = self.GetSelectedEvent()
        if EventName:
            self.pltrec.PlotPSDSNR(EventName,
                                     TDelay=Estart,
                                     TEval=Estop,
                                     DevDCVals=self.DevDCVals,
                                     Time=(Tstart, Tstop),
                                     nFFT=nFFT)


def main():
    import argparse
    import pkg_resources

    # Add version option
    __version__ = pkg_resources.require("PyGFET")[0].version
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(
                            version=__version__))
    parser.parse_args()

    app = QtWidgets.QApplication(sys.argv)
    w = RecordExplorer()
    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
