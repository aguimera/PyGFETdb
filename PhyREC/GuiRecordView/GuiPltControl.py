#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 17:06:29 2019

@author: aguimera
"""

import sys
import os

from qtpy import QtWidgets, uic
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QTreeWidgetItem

from PhyREC.NeoInterface import NeoSegment, NeoSignal
import PhyREC.PlotWaves as Rplt
import PhyREC.SignalProcess as Rpro
import PhyREC.SignalAnalysis as Ran
import matplotlib.pyplot as plt
import quantities as pq

from matplotlib.artist import ArtistInspector

AxesProp = ('ylim',
            'ylabel',
            'xlabel',
            'visible',
            'autoscaley_on')

AxisProp = ('visible',
            'scale',
            )


class GuiPltControl(QtWidgets.QMainWindow):

    def __init__(self, PlotSlot):
        QtWidgets.QMainWindow.__init__(self)
        uipath = os.path.join(os.path.dirname(__file__),
                              'PltControl.ui')
        uic.loadUi(uipath, self)

        self.setWindowTitle('Plot Control')
        
        self.TreeCtr.setColumnCount(3)
        self.TreeCtr.setColumnWidth(0, 150)
        self.TreeCtr.setColumnHidden(2, True)

        self.PlotSlot = PlotSlot

        self.FillTreeView()
        self.InitTimeControl()

        self.TreeCtr.itemChanged.connect(self.ItemChanged)

        self.SlideTStart.valueChanged.connect(self.SlideTStartChange)
        self.SlideTStop.valueChanged.connect(self.SlideTStopChange)

    def SlideTStartChange(self):
        ts = self.SlideTStart.value()
        Tshow = self.PlotSlot.current_time[1] - self.PlotSlot.current_time[0]
        self.SlideTStop.setValue(ts+Tshow)
#        self.PlotSlot.PlotChannels((ts*pq.s, (ts+Tshow)*pq.s))
#        self.PlotSlot.Fig.canvas.draw()

    def SlideTStopChange(self):
        tst = self.SlideTStop.value()
        ts = self.SlideTStart.value()
        if ts > tst:
            self.SlideTStop.setValue(ts+1)
            return
        self.PlotSlot.PlotChannels((ts*pq.s, tst*pq.s))
        self.PlotSlot.Fig.canvas.draw()
        

    def InitTimeControl(self):
        tstarts = [sl.Signal.t_start for sl in self.PlotSlot.Slots]
        tstops = [sl.Signal.t_stop for sl in self.PlotSlot.Slots]
        self.SlideTStart.setMinimum(min(tstarts))
        self.SlideTStart.setMaximum(max(tstops))
        self.SlideTStop.setMinimum(min(tstarts))
        self.SlideTStop.setMaximum(max(tstops))

        self.SlideTStart.setValue(self.PlotSlot.current_time[0])
        self.SlideTStop.setValue(self.PlotSlot.current_time[1])

    def FillTreeView(self):
        ItemP = QTreeWidgetItem(self.TreeCtr)
        ItemP.setText(0, 'Figure')
        font = ItemP.font(0)
        font.setBold(True)
        font.setItalic(True)
        ItemP.setFont(0, font)
        for k, v in self.PlotSlot.FigKwargs.items():
            ItemCh = QTreeWidgetItem(ItemP)
            ItemCh.setFlags(ItemCh.flags() | Qt.ItemIsEditable)
            self.AddItemData(ItemCh, k, v)

        ItemP = QTreeWidgetItem(self.TreeCtr)
        ItemP.setText(0, 'Legend')
        font = ItemP.font(0)
        font.setBold(True)
        font.setItalic(True)
        ItemP.setFont(0, font)
        ItemP.setFlags(ItemP.flags() | Qt.ItemIsEditable)
        ItemP.setData(1, Qt.EditRole, True)
        for k, v in self.PlotSlot.LegendKwargs.items():
            ItemCh = QTreeWidgetItem(ItemP)
            ItemCh.setFlags(ItemCh.flags() | Qt.ItemIsEditable)
            self.AddItemData(ItemCh, k, v)

        for Axi, Ax in enumerate(self.PlotSlot.Axs):
            Slots = self.PlotSlot.SlotsInAxs[Ax]

            ItemP = QTreeWidgetItem(self.TreeCtr)
            ItemP.setText(0, 'Axes')
            font = ItemP.font(0)
            font.setBold(True)
            font.setItalic(True)
            ItemP.setFont(0, font)
            Values = Slots[0].AxKwargs
            ItemP.setData(2, Qt.EditRole, Axi)
            self.FillDictValues(ItemParent=ItemP,
                                PltObj=Ax,
                                Kwargs=Values)

            for isl, sl in enumerate(Slots):
                ItemP = QTreeWidgetItem(ItemP)
                ItemP.setText(0, 'Wave')
                font = ItemP.font(0)
                font.setBold(True)
                font.setItalic(True)
                ItemP.setFont(0, font)                
                ItemP.setData(2, Qt.EditRole, isl)
                Values = sl.LineKwargs
                self.FillDictValues(ItemParent=ItemP,
                                    PltObj=sl.Line,
                                    Kwargs=Values)

    def FillDictValues(self, ItemParent, PltObj, Kwargs):
        if PltObj is not None:
            ains = ArtistInspector(PltObj)
            validp = ains.get_setters()
        else:
            validp = list(Kwargs.keys())

        for p in Kwargs.keys():
            v = getattr(PltObj, 'get_' + p)()
            if p in validp:
                ItemV = QTreeWidgetItem(ItemParent)
                ItemV.setFlags(ItemV.flags() | Qt.ItemIsEditable)
                self.AddItemData(ItemV, p, v)
            else:
                ItemV = QTreeWidgetItem(ItemParent)
                ItemV.setText(0, p)
                self.FillDictValues(ItemParent=ItemV,
                                    PltObj=v,
                                    Kwargs=Kwargs[p])

    def AddItemData(self, Item, DataName, DataValue):
        Item.setText(0, DataName)
        if type(DataValue) == tuple:
            Item.setText(2, type(DataValue).__name__)
            for iv, val in enumerate(DataValue):
                Child = QTreeWidgetItem(Item)
                Child.setFlags(Child.flags() | Qt.ItemIsEditable)
                self.AddItemData(Child, str(iv), val)
            return
        elif type(DataValue) == str:
            D = DataValue
        elif type(DataValue) == bool:
            D = DataValue
        else:
            D = float(DataValue)
        Item.setData(1, Qt.EditRole, D)

    def ItemChanged(self, item, column):
        Parent = item.parent()
        Parents = []
        while Parent is not None:
            Parents.append((Parent.text(0), Parent))
            Parent = Parent.parent()

        if len(Parents) == 0:
            if item.text(0) == 'Legend':
                if item.data(1, Qt.EditRole) is False:
                    for ax in self.PlotSlot.Axs:
                        ax.get_legend().remove()
                    self.PlotSlot.Fig.canvas.draw()
            return

        firstp = Parents[0][1]
        if firstp.text(2) == 'tuple':
            data = []
            firstp.sortChildren(0, Qt.AscendingOrder)
            for ic in range(firstp.childCount()):
                data.append(firstp.child(ic).data(1, Qt.EditRole))
            PropName = Parents[0][0]
            Parents.pop(0)
        else:
            data = item.data(1, Qt.EditRole)
            PropName = item.text(0)

        DataDict = {PropName: data}
        for p in Parents[:-1]:
            DataDict = {p[0]: DataDict}
#        print(DataDict)

        kwtype = Parents[-1][0]
        iAx = Parents[-1][1].data(2, Qt.EditRole)
        if kwtype == 'Axes':
            Ax = self.PlotSlot.Axs[iAx]
            if Parents[0][0] == 'Wave':
                iSl = Parents[0][1].data(2, Qt.EditRole)
                Sl = self.PlotSlot.SlotsInAxs[Ax][iSl]
                Sl.UpdateLineKwargs(DataDict['Wave'])
            else:
                Slots = self.PlotSlot.SlotsInAxs[Ax]
                for Sl in Slots:
                    Sl.UpdateAxKwargs(DataDict)
        elif kwtype == 'Figure':
            self.PlotSlot.UpdateFigKwargs(DataDict)
        elif kwtype == 'Legend':
            if Parents[-1][1].data(1, Qt.EditRole):
                self.PlotSlot.AddLegend(**DataDict)

        self.PlotSlot.Fig.canvas.draw()


#def IsNumber(s):
#    try:
#        int(s)
#        return True
#    except ValueError:
#        return False


def main(PlotSlot):
    import argparse
    import pkg_resources

    # Add version option
    __version__ = pkg_resources.require("PhyREC")[0].version
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    parser.parse_args()

    app = QtWidgets.QApplication(sys.argv)
    w = GuiPltControl(PlotSlot)
    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":

    FileIn = 'TestSig.h5'
    Twind = (1*pq.s, None)
    Range = (-50, 50)

    aiChannels = (
                  'Ch05',
                  'Ch06',
                  'Ch07',
                  'Ch08',
                  'Ch13',
                  'Ch14',
                  'Ch15',
                  'Ch16',
                  )
    
    doColumns = ('Col1',
                 'Col2',
                 'Col3',
                 'Col4',
                 'Col5',
                 'Col6',
                 'Col7',
                 'Col8',
                 )  


    SigCond = ({'function': Rpro.Filter, 'args': {'Type':'bandstop',
                                                  'Order':2,
                                                  'Freqs':(49, 51),}},
               {'function': Rpro.Filter, 'args': {'Type':'highpass',
                                                  'Order':2,
                                                  'Freqs':(1, ),}},   
                                                  )

    Chorder = {}
    for irow, row in enumerate(aiChannels):
        for icol, col in enumerate(doColumns):
            Chorder[row+col] = (irow, icol)
            
    Rec = NeoSegment(FileIn)

    AxesProp = {
                'ylim': Range,
                'facecolor': '#FFFFFF00',
                'autoscaley_on': False,
                'xaxis': {'visible': False,
                          },
                'yaxis': {'visible': False,
                          },
                'ylabel':'',
                'title':None,
                }

    FigProp = {'tight_layout': True,
               'size_inches': (10,5),
#               'facecolor': '#FFFFFF00',
               }

#    fig, axs = plt.subplots(len(aiChannels), len(doColumns),
#                            sharex=True)
    Slots = []
    isig = 0
    for sig in Rec.Signals()[0:10]:
        if not sig.name.endswith('AC'):
            continue        
        chname = sig.name.split('_')[0]
        sig.ProcessChain = SigCond
        Slots.append(Rplt.WaveSlot(sig,
                                   Units='nA',
                                   Position=isig,
#                                   Ax=axs[Chorder[chname]],
#                                   AxKwargs=AxesProp,
                                   color='r',
                                   alpha=0.5))
        isig += 1

    splt = Rplt.PlotSlots(Slots,
                          AxKwargs=AxesProp,
                          FigKwargs=FigProp)

    splt.PlotChannels(Time=Twind)

    splt.AddLegend()

#
#test = Rplt.WaveSlot(sig,
#                                   Units='nA',
#                                   Position=isig,
##                                   Ax=axs[Chorder[chname]],
##                                   AxKwargs=AxesProp,
#                                   color='r',
#                                   alpha=0.5)
##
##
##    ains = ArtistInspector(fig)
##    validprop = ains.get_setters()
##    
#

    main(splt)







