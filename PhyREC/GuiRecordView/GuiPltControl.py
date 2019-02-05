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
        self.TreeCtr.itemChanged.connect(self.ItemChanged)

        self.PlotSlot = PlotSlot

        ItemP = QTreeWidgetItem(self.TreeCtr)
        ItemP.setText(0, 'Figure')
        for k, v in self.PlotSlot.FigKwargs.items():
            ItemCh = QTreeWidgetItem(ItemP)
            ItemCh.setFlags(ItemCh.flags() | Qt.ItemIsEditable)
            self.AddItemData(ItemCh, k, v)

        for Ax, Slots in self.PlotSlot.SlotsInAxs.items():
            ItemP = QTreeWidgetItem(self.TreeCtr)
            ItemP.setText(0, 'Axes')
            Values = Slots[0].AxKwargs
            self.FillDictValues(ItemParent=ItemP,
                                PltObj=Ax,
                                Kwargs=Values)
            for sl in Slots:
                ItemP = QTreeWidgetItem(ItemP)
                ItemP.setText(0, 'Wave')
                Values = sl.LineKwargs
                self.FillDictValues(ItemParent=ItemP,
                                    PltObj=sl.Line,
                                    Kwargs=Values)

    def FillDictValues(self, ItemParent, PltObj, Kwargs):
        ains = ArtistInspector(PltObj)
        validp = ains.get_setters()
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
        pass
#        Parent = item.parent()
#        ParentText = Parent.text(0)
#        if ParentText == 'Axes':            
#            Axi = Parent.data(2, Qt.EditRole)
#            prop = item.text(0)
#            propv = item.data(1, Qt.EditRole)
#            self.PlotSlot.Axs[Axi].set(**{prop: propv})
#            print(prop, propv, Axi)
#        elif ParentText == 'ylim':            
#            vals = []
#            for i in range(Parent.childCount()):
#                ch = Parent.child(i)
#                vals.append(ch.data(1, Qt.EditRole))
#            
#            AxParent = Parent.parent()
#            Axi = AxParent.data(2, Qt.EditRole)
##            print ('ylim', vals, Axi, self.PlotSlot.Axs[Axi])
#            self.PlotSlot.Axs[Axi].set(**{'ylim': vals})
#        elif ParentText == 'XAxis':
#            AxParent = Parent.parent()
#            Axi = AxParent.data(2, Qt.EditRole)
#            prop = item.text(0)
#            propv = item.data(1, Qt.EditRole)
#            xax = self.PlotSlot.Axs[Axi].get_xaxis()
#            xax.set(**{prop: propv})
#        elif ParentText == 'YAxis':
#            AxParent = Parent.parent()
#            Axi = AxParent.data(2, Qt.EditRole)
#            prop = item.text(0)
#            propv = item.data(1, Qt.EditRole)
#            xax = self.PlotSlot.Axs[Axi].get_yaxis()
#            xax.set(**{prop: propv})
#        elif ParentText.startswith('WaveSlot'):
#            print(ParentText)
#            
        

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
#
    main(splt)







