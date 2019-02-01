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

        self.PlotSlot = PlotSlot

#            parent.setFlags(parent.flags() | Qt.ItemIsTristate | Qt.ItemIsUserCheckable)
        for iax, ax in enumerate(PlotSlot.Axs):
            AxParent = QTreeWidgetItem(self.TreeCtr)
            AxParent.setText(0, "Axes")
            AxParent.setText(1, PlotSlot.SlotsInAxs[ax][0].name)
            AxParent.setData(2, Qt.EditRole, iax)
            
            self.AddAxProp(ParentItem=AxParent,
                           PropsToAdd=AxesProp,
                           Properties=ax.properties())

            AxYParent = QTreeWidgetItem(AxParent)
            AxYParent.setText(0, "YAxis")
            self.AddAxProp(ParentItem=AxYParent,
                           PropsToAdd=AxisProp,
                           Properties=ax.get_yaxis().properties())
            
            AxXParent = QTreeWidgetItem(AxParent)
            AxXParent.setText(0, "XAxis")
            self.AddAxProp(ParentItem=AxXParent,
                           PropsToAdd=AxisProp,
                           Properties=ax.get_yaxis().properties())
            
            for isl, sl in enumerate(PlotSlot.SlotsInAxs[ax]):
                WvChild = QTreeWidgetItem(AxParent)
                WvChild.setText(0, str(type(sl)).split('.')[-1])
                WvChild.setText(1, sl.name)
                WvChild.setData(2, Qt.EditRole, isl)
                for sp, sv in sl.LineKwargs.items():
                    WvPChild = QTreeWidgetItem(WvChild)
                    WvPChild.setFlags(WvPChild.flags() | Qt.ItemIsEditable)
                    self.AddItemData(WvPChild, sp, sv)

        self.TreeCtr.itemChanged.connect(self.ItemChanged)

    def AddAxProp(self, ParentItem, PropsToAdd, Properties):
        for p in PropsToAdd:
            AxChild = QTreeWidgetItem(ParentItem)
            AxChild.setFlags(AxChild.flags() | Qt.ItemIsEditable)
            self.AddItemData(AxChild, p, Properties[p])

    def AddItemData(self, Item, DataName, DataValue):
        Item.setText(0, DataName)
        if type(DataValue)==tuple:
            for iv, val in enumerate(DataValue):
                Child = QTreeWidgetItem(Item)
                Child.setFlags(Child.flags() | Qt.ItemIsEditable)
                self.AddItemData(Child, str(iv), val)
            return
        elif type(DataValue)==str:
            D = DataValue
        elif type(DataValue)==bool:
            D = DataValue
        else:
            D = float(DataValue)
        Item.setData(1, Qt.EditRole, D)

    def ItemChanged(self, item, column):
        Parent = item.parent()
        ParentText = Parent.text(0)
        if ParentText == 'Axes':            
            Axi = Parent.data(2, Qt.EditRole)
            prop = item.text(0)
            propv = item.data(1, Qt.EditRole)
            self.PlotSlot.Axs[Axi].set(**{prop: propv})
            print(prop, propv, Axi)
        elif ParentText == 'ylim':            
            vals = []
            for i in range(Parent.childCount()):
                ch = Parent.child(i)
                vals.append(ch.data(1, Qt.EditRole))
            
            AxParent = Parent.parent()
            Axi = AxParent.data(2, Qt.EditRole)
#            print ('ylim', vals, Axi, self.PlotSlot.Axs[Axi])
            self.PlotSlot.Axs[Axi].set(**{'ylim': vals})
        elif ParentText == 'XAxis':
            AxParent = Parent.parent()
            Axi = AxParent.data(2, Qt.EditRole)
            prop = item.text(0)
            propv = item.data(1, Qt.EditRole)
            xax = self.PlotSlot.Axs[Axi].get_xaxis()
            xax.set(**{prop: propv})
        elif ParentText == 'YAxis':
            AxParent = Parent.parent()
            Axi = AxParent.data(2, Qt.EditRole)
            prop = item.text(0)
            propv = item.data(1, Qt.EditRole)
            xax = self.PlotSlot.Axs[Axi].get_yaxis()
            xax.set(**{prop: propv})
        elif ParentText.startswith('WaveSlot'):
            print(ParentText)
            
        

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
    
    
    a = splt.Axs[0]
    prop = a.properties()
    xax = a.get_xaxis()




    

#
#
#    ains = ArtistInspector(fig)
#    validprop = ains.get_setters()
#    


#    main(splt)







