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


AxProp = ('ylim',
              'ylabel',
              'visible',
              'autoscaley_on')

AxXprop = ('visible', )

AxYprop = ('visible', )

class GuiPltControl(QtWidgets.QMainWindow):
    def __init__(self, PlotSlot):
        QtWidgets.QMainWindow.__init__(self)
        uipath = os.path.join(os.path.dirname(__file__),
                              'PltControl.ui')
        uic.loadUi(uipath, self)

        self.setWindowTitle('Plot Control')
        self.TreeCtr.setColumnCount(2)

#            parent.setFlags(parent.flags() | Qt.ItemIsTristate | Qt.ItemIsUserCheckable)
        for iax, ax in enumerate(PlotSlot.Axs):
            AxParent = QTreeWidgetItem(self.TreeCtr)
            AxParent.setText(0, "Axes {}".format(iax))

            props = ax.properties()
            for p in AxProp:
                AxChild = QTreeWidgetItem(AxParent)
                AxChild.setFlags(AxChild.flags() | Qt.ItemIsEditable)
                AxChild.setText(0, p)
                pval = props[p]
                if type(pval) == bool:
                    AxChild.setCheckState(1, pval)
                else:
                    AxChild.setData(1, Qt.EditRole, pval)
            
            for sl in PlotSlot.SlotsInAxs[ax]:
                WvChild = QTreeWidgetItem(AxParent)
                WvChild.setText(0, sl.name)
                for sp, sv in sl.LineKwargs.items():
                    WvPChild = QTreeWidgetItem(WvChild)
                    WvPChild.setFlags(WvPChild.flags() | Qt.ItemIsEditable)
                    
                    WvPChild.setText(0, sp)
#                    WvPChild.setText(1, str(sv))
#                    if type(sv) == bool:
#                        WvPChild.setCheckState(1, sv)
#                    else:
                    WvPChild.setData(1, Qt.EditRole, sv)

        self.TreeCtr.itemChanged.connect(self.ItemChanged)

    def ItemChanged(self, item, column):
        print(item.text(0), item.data(1, Qt.EditRole))
        

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

    fig, axs = plt.subplots(len(aiChannels), len(doColumns),
                            sharex=True)

    Slots = []
    for sig in Rec.Signals():
        if not sig.name.endswith('AC'):
            continue
        chname = sig.name.split('_')[0]
        sig.ProcessChain = SigCond
        Slots.append(Rplt.WaveSlot(sig,
                                   Units='nA',
                                   Ax=axs[Chorder[chname]],
                                   Fig=fig,
                                   Ylim=Range,
                                   color='r'))


    splt = Rplt.PlotSlots(Slots,
                          Fig=fig,
                          ShowAxis=None,
   #                      ShowNameOn='Legend',
                          AutoScale=False)

    splt.PlotChannels(Time=Twind)

    main(splt)







