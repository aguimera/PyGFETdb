# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import sys

import matplotlib.pyplot as plt
import pickle

from qtpy.QtWidgets import QHeaderView, QMessageBox
from qtpy.QtWidgets import QFileDialog, QAction, QInputDialog
from qtpy import QtWidgets, uic
from qtpy.QtCore import Qt, QItemSelectionModel

import PyGFETdb.DBCore as PyFETdb
import PyGFETdb.AnalyzeData as PyFETData
import PyGFETdb.PlotDataClass as PyFETpl


class DBViewApp(QtWidgets.QMainWindow):
    OutFigFormats = ('svg', 'png')
    # ('IdTrt',0,0) (Header, position, string Conv, Editable)
    TrtFields = {'Trts.idTrts': ('IdTrt', 0, 0, 0),
                 'Trts.Name': ('Name', 1, 0, 0),
                 'TrtTypes.Name': ('Trt Type', 2, 0, 0),
                 'TrtTypes.Length': ('Lenght', 3, 0, 0),
                 'TrtTypes.Width': ('Width', 4, 0, 0),
                 'VTrts.DCMeas': ('DCMeas', 5, 0, 0),
                 'VTrts.ACMeas': ('ACMeas', 6, 0, 0),
                 'TrtTypes.Contact': ('Contact', 7, 0, 0),
                 'TrtTypes.Pass': ('Pass', 8, 0, 0),
                 'TrtTypes.Area': ('Area', 9, 0, 0),
                 'Trts.Comments': ('T-Comments', 10, 0, 1),
                 'Devices.Comments': ('D-Comments', 11, 0, 1),
                 'Wafers.Comments': ('W-Comments', 12, 0, 1),
                 'Devices.State': ('D-State', 13, 0, 1),
                 'Devices.ExpOK': ('D-ExpOK', 14, 0, 1),
                 'Devices.idDevices': ('idDevices', 15, 0, 0),
                 'Wafers.idWafers': ('idWafers', 16, 0, 0),
                 'VTrts.GMeas': ('GMeas', 17, 0, 0)}

    # ('IdTrt',0,0) (Column : (Table, Field, condfield, col cond val))
    TrtsUpdateFields = {10: ('Trts', 'Comments', 'idTrts=', 0),
                        11: ('Devices', 'Comments', 'idDevices=', 15),
                        12: ('Wafers', 'Comments', 'idWafers=', 16),
                        13: ('Devices', 'State', 'idDevices=', 15),
                        14: ('Devices', 'ExpOK', 'idDevices=', 15)}

    TrtSearchFields = list(TrtFields.keys())
    TrtSearchFields.append('Devices.Name')
    TrtSearchFields.append('TrtTypes.Name')
    TrtSearchFields.append('Wafers.Name')
    TrtSearchFields.append('Wafers.Substrate')
    TrtSearchFields.append('Wafers.Masks')

    # ('IdTrt',0,0) (Header, position, string Conv, Editable)
    DCFields = {'DCcharacts.idDCcharacts': ('IdDC', 0, 0, 0),
                'Trts.Name': ('Name', 1, 0, 0),
                'DCcharacts.MeasDate': ('Date', 2, 1, 0),
                'DCcharacts.IsOK': ('IsOK', 3, 0, 1),
                'DCcharacts.IsCmp': ('IsCmp', 4, 0, 1),
                'DCcharacts.Ph': ('Ph', 5, 0, 1),
                'DCcharacts.IonStrength': ('IonStrength', 6, 0, 1),
                'DCcharacts.FuncStep': ('FuncStep', 7, 0, 1),
                'DCcharacts.AnalyteCon': ('AnalyteCon', 8, 1, 1),
                'DCcharacts.Comments': ('Comments', 9, 0, 1),
                'DCcharacts.UpdateDate': ('UpdateDate', 10, 1, 0),
                'DCcharacts.FileName': ('FileName', 11, 0, 0)}

    # ('IdTrt',0,0) (Column : Field)
    DCUpdateFields = {3: 'IsOK',
                      4: 'IsCmp',
                      5: 'Ph',
                      6: 'IonStrength',
                      7: 'FuncStep',
                      8: 'AnalyteCon',
                      9: 'Comments'}

    # ('IdTrt',0,0) (Header, position, string Conv, Editable)
    ACFields = {'ACcharacts.idACcharacts': ('IdAC', 0, 0, 0),
                'Trts.Name': ('Name', 1, 0, 0),
                'ACcharacts.MeasDate': ('Date', 2, 1, 0),
                'ACcharacts.IsOK': ('IsOK', 3, 0, 1),
                'ACcharacts.IsCmp': ('IsCmp', 4, 0, 1),
                'ACcharacts.Ph': ('Ph', 5, 0, 1),
                'ACcharacts.IonStrength': ('IonStrength', 6, 1, 1),
                'ACcharacts.FuncStep': ('FuncStep', 7, 0, 1),
                'ACcharacts.AnalyteCon': ('AnalyteCon', 8, 0, 1),
                'ACcharacts.Comments': ('Comments', 9, 0, 1),
                'ACcharacts.UpdateDate': ('UpdateDate', 10, 1, 0),
                'ACcharacts.DC_id': ('IdDC', 11, 0, 0),
                'ACcharacts.FileName': ('FileName', 12, 0, 0)}

    # ('IdTrt',0,0) (Column : Field)
    ACUpdateFields = {3: 'IsOK',
                      4: 'IsCmp',
                      5: 'Ph',
                      6: 'IonStrength',
                      7: 'FuncStep',
                      8: 'AnalyteCon',
                      9: 'Comments'}

    ViewAxsAC = ('GmMag', 'IdsPoly', 'GMPoly')
    ViewAxsDC = ('Gm', 'Ids', 'Rds')

    def InitMenu(self):
        mainMenu = self.menubar
        fileMenu = mainMenu.addMenu('File')

        SaveFigAction = QAction('Save Figures', self)
        SaveFigAction.setShortcut('Ctrl+s')
        SaveFigAction.setStatusTip('Save all open figures')
        SaveFigAction.triggered.connect(self.SaveFigures)
        fileMenu.addAction(SaveFigAction)

        CloseFigsAction = QAction('Close Figures', self)
        CloseFigsAction.setStatusTip('Close all open figures')
        CloseFigsAction.triggered.connect(self.CloseFigures)
        fileMenu.addAction(CloseFigsAction)

        toolsMenu = mainMenu.addMenu('Tools')
        DevRepAction = QAction('Xls Device Report', self)
        DevRepAction.setStatusTip('Generate Device Report')
        DevRepAction.triggered.connect(self.DevReport)
        toolsMenu.addAction(DevRepAction)

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)

        uipath = os.path.join(os.path.dirname(__file__), 'GuiDBView.ui')
        uic.loadUi(uipath, self)

        self.setWindowTitle('PyFETdb Viewer')

        self.DB = PyFETdb.PyFETdb()

        self.InitMenu()
        self.ConnectLst()

        self.ButResetSearch.clicked.connect(self.ButResetSearchClick)
        self.ButSetData.clicked.connect(self.ButSetDataClick)

        self.ButViewDC.clicked.connect(self.ButViewDCClick)
        self.ButExportDC.clicked.connect(self.ButExportDCClick)

        self.ButAnalyzeAC.clicked.connect(self.ButAnalyzeACClick)
        self.ButExportAC.clicked.connect(self.ButExportACClick)
        self.ButViewAC.clicked.connect(self.ButViewACClick)

        self.ButDeleteAC.clicked.connect(self.ButDeleteACClick)
        self.TblAC.cellChanged.connect(self.TblACCellChanged)
        self.TblDC.cellChanged.connect(self.TblDCCellChanged)
        self.TblTrts.cellChanged.connect(self.TblTrtsCellChanged)

        self.ButDeleteDC.clicked.connect(self.ButDeleteDCClick)
        self.ButAnalyzeDC.clicked.connect(self.ButAnalyzeDCClick)

        # Init Tables
        self.TblTrts.setColumnCount(len(self.TrtFields))
        for c in self.TrtFields.values():
            self.TblTrts.setHorizontalHeaderItem(c[1],
                                                 QtWidgets.QTableWidgetItem(c[0]))
        header = self.TblTrts.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeToContents)

        # Init Tables
        self.TblDC.setColumnCount(len(self.DCFields))
        for c in self.DCFields.values():
            self.TblDC.setHorizontalHeaderItem(c[1],
                                               QtWidgets.QTableWidgetItem(c[0]))
        header = self.TblDC.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeToContents)

        # Init Tables
        self.TblAC.setColumnCount(len(self.ACFields))
        for c in self.ACFields.values():
            self.TblAC.setHorizontalHeaderItem(c[1],
                                               QtWidgets.QTableWidgetItem(c[0]))
        header = self.TblAC.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeToContents)

        # Init
        self.ButResetSearchClick()

    def ConnectLst(self):
        self.LstWafers.itemSelectionChanged.connect(self.LstWafersChange)
        self.LstSubstrates.itemSelectionChanged.connect(self.LstSubstratesChange)
        self.LstMasks.itemSelectionChanged.connect(self.LstMasksChange)
        self.LstDevices.itemSelectionChanged.connect(self.LstDevicesChange)
        self.LstTypes.itemSelectionChanged.connect(self.LstTypesChange)
        self.LstL.itemSelectionChanged.connect(self.LstLChange)
        self.LstW.itemSelectionChanged.connect(self.LstWChange)
        self.LstContact.itemSelectionChanged.connect(self.LstContactChange)
        self.LstPass.itemSelectionChanged.connect(self.LstPassChange)
        self.LstArea.itemSelectionChanged.connect(self.LstAreaChange)

    def DisconnectLst(self):
        self.LstWafers.itemSelectionChanged.disconnect(self.LstWafersChange)
        self.LstSubstrates.itemSelectionChanged.disconnect(self.LstSubstratesChange)
        self.LstMasks.itemSelectionChanged.disconnect(self.LstMasksChange)
        self.LstDevices.itemSelectionChanged.disconnect(self.LstDevicesChange)
        self.LstTypes.itemSelectionChanged.disconnect(self.LstTypesChange)
        self.LstL.itemSelectionChanged.disconnect(self.LstLChange)
        self.LstW.itemSelectionChanged.disconnect(self.LstWChange)
        self.LstContact.itemSelectionChanged.disconnect(self.LstContactChange)
        self.LstPass.itemSelectionChanged.disconnect(self.LstPassChange)
        self.LstArea.itemSelectionChanged.disconnect(self.LstAreaChange)

    def DevReport(self):        
        Vals = []
        for r in self.Trts:
            Vals.append(r['Wafers.Name'])

        Prefix, okPressed = QInputDialog.getText(self,
                                                 'Prefix',
                                                 'Prefix for files',
                                                 text='Figure')

        print (list(set(Vals)))
        
        FileName = QFileDialog.getSaveFileName(self,
                                               "Output File", "",
                                               "Excel (*.xlsx)")

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
                    print (fileOut)
                    plt.savefig(fileOut.format(i))

    def ButResetSearchClick(self):
        self.Cond = {}
        self.DisconnectLst()
        self.UpdateWafersList()
        self.UpdateSearchList(Table=False)
        self.ConnectLst()

    def UpdateWafersList(self, Substrates=True, Masks=True, Wafers=True):
        if len(self.Cond) == 0:
            self.Trts = self.DB.GetTrtsInfo({'Trts.idTrts>': (0, )},
                                            self.TrtSearchFields)
        else:
            self.Trts = self.DB.GetTrtsInfo(self.Cond,
                                            self.TrtSearchFields)

        if Substrates:
            self.FillList(self.LstSubstrates, 'Wafers.Substrate')
        if Masks:
            self.FillList(self.LstMasks, 'Wafers.Masks')
        if Wafers:
            self.FillList(self.LstWafers, 'Wafers.Name')

    def LstSubstratesChange(self):
        sel = self.LstSubstrates.selectedItems()
        if len(sel) == 0:
            return

        self.Cond['Wafers.Substrate='] = []
        for s in sel:
            a = s.text()
            if s.text() == 'NONE':
                a = 'NULL'
            self.Cond['Wafers.Substrate='].append(a)
        self.UpdateWafersList(Substrates=False)

    def LstMasksChange(self):
        sel = self.LstMasks.selectedItems()
        if len(sel) == 0:
            return

        self.Cond['Wafers.Masks='] = []
        for s in sel:
            a = s.text()
            if s.text() == 'NONE':
                a = 'NULL'
            self.Cond['Wafers.Masks='].append(a)
        self.UpdateWafersList(Masks=False)

    def LstWafersChange(self):
        sel = self.LstWafers.selectedItems()
        if len(sel) == 0:
            return

        self.Cond = {}
        self.Cond['Wafers.Name='] = []
        for s in sel:
            self.Cond['Wafers.Name='].append(s.text())
        self.UpdateSearchList()

    def UpdateSearchList(self, Devices=True, Types=True, Length=True,
                         Width=True, Pass=True, Contact=True, Area=True,
                         Table=True):

        if len(self.Cond) == 0:
            self.Trts = self.DB.GetTrtsInfo({'Trts.idTrts>': (0, )},
                                            self.TrtSearchFields)
        else:
            self.Trts = self.DB.GetTrtsInfo(self.Cond,
                                            self.TrtSearchFields)

        if Devices:
            self.FillList(self.LstDevices, 'Devices.Name')
        if Types:
            self.FillList(self.LstTypes, 'TrtTypes.Name')
        if Length:
            self.FillList(self.LstL, 'TrtTypes.Length')
        if Width:
            self.FillList(self.LstW, 'TrtTypes.Width')
        if Pass:
            self.FillList(self.LstPass, 'TrtTypes.Pass')
        if Contact:
            self.FillList(self.LstContact, 'TrtTypes.Contact')
        if Area:
            self.FillList(self.LstArea, 'TrtTypes.Area')

        if Table:
            self.ChkUpdateTrts.setChecked(False)
            self.FillTable(self.TblTrts, self.Trts, self.TrtFields)

    def FillList(self, Lst, Field):
        Vals = []
        for r in self.Trts:
            Vals.append(r[Field])
        Lst.clear()
        try:
            for d in sorted(set(Vals)):
                Lst.addItem(str(d))
        except:
            for d in set(Vals):
                Lst.addItem(str(d))

    def LstDevicesChange(self):
        sel = self.LstDevices.selectedItems()
        if len(sel) == 0:
            return

        self.Cond['Devices.Name='] = []
        for s in sel:
            self.Cond['Devices.Name='].append(s.text())
        self.UpdateSearchList(Devices=False)

    def LstTypesChange(self):
        sel = self.LstTypes.selectedItems()
        if len(sel) == 0:
            return

        self.Cond['TrtTypes.Name='] = []
        for s in sel:
            self.Cond['TrtTypes.Name='].append(s.text())
        self.UpdateSearchList(Devices=False, Types=False)

    def LstLChange(self):
        sel = self.LstL.selectedItems()
        if len(sel) == 0:
            return

        self.Cond['TrtTypes.Length='] = []
        for s in sel:
            self.Cond['TrtTypes.Length='].append(s.text())
        self.UpdateSearchList(Devices=False, Length=False)

    def LstWChange(self):
        sel = self.LstW.selectedItems()
        if len(sel) == 0:
            return

        self.Cond['TrtTypes.Width='] = []
        for s in sel:
            self.Cond['TrtTypes.Width='].append(s.text())
        self.UpdateSearchList(Devices=False, Width=False)

    def LstContactChange(self):
        sel = self.LstContact.selectedItems()
        if len(sel) == 0:
            return

        self.Cond['TrtTypes.Contact='] = []
        for s in sel:
            self.Cond['TrtTypes.Contact='].append(s.text())
        self.UpdateSearchList(Devices=False, Contact=False)

    def LstPassChange(self):
        sel = self.LstPass.selectedItems()
        if len(sel) == 0:
            return

        self.Cond['TrtTypes.Pass='] = []
        for s in sel:
            self.Cond['TrtTypes.Pass='].append(s.text())
        self.UpdateSearchList(Devices=False, Pass=False)

    def LstAreaChange(self):
        sel = self.LstArea.selectedItems()
        if len(sel) == 0:
            return

        self.Cond['TrtTypes.Area='] = []
        for s in sel:
            self.Cond['TrtTypes.Area='].append(s.text())
        self.UpdateSearchList(Devices=False, Area=False)

    def ButSetDataClick(self):
        Ids = self.GetTableSelectCol(self.TblTrts)

        self.TblDC.cellChanged.disconnect()
        self.FillCharactTable(self.TblDC, 'DCcharacts', self.DCFields, Ids)
        self.TblDC.cellChanged.connect(self.TblDCCellChanged)

        self.TblAC.cellChanged.disconnect()
        self.FillCharactTable(self.TblAC, 'ACcharacts', self.ACFields, Ids)
        self.TblAC.cellChanged.connect(self.TblACCellChanged)

    def FillCharactTable(self, Tbl, DBTable, Fields, Ids):
        Tbl.setSortingEnabled(False)
        Tbl.setRowCount(0)
        for ids in Ids:
            rows = self.DB.GetCharactInfo(DBTable,
                                          {'{}.Trt_id='.format(DBTable, DBTable):(ids, )},
                                          Fields.keys())
            self.FillTable(Tbl, rows, Fields, Offset=Tbl.rowCount())

    def FillTable(self, Tbl, Rows, Fields, Offset=0):

        Tbl.setSortingEnabled(False)

        Tbl.setRowCount(len(Rows)+Offset)
        for ir, r in enumerate(Rows):
            for f, c in Fields.items():
                item = QtWidgets.QTableWidgetItem()

                if c[3] == 0:
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)

                if c[2] == 1:
                    item.setData(Qt.DisplayRole, str(r[f]))
                else:
                    item.setData(Qt.DisplayRole, r[f])
                Tbl.setItem(ir+Offset, c[1], item)

        Tbl.setSortingEnabled(True)

    def GetTableSelectCol(self, Table, Col=0, String=False):
        ids = []
        # Get selected
        Selection = Table.selectedRanges()
        for sel in Selection:
            for i in range(sel.rowCount()):
                if String:
                    ids.append(Table.item(sel.topRow()+i, Col).text())
                else:
                    ids.append(int(Table.item(sel.topRow()+i, Col).text()))
        return ids

    def ButDeleteDCClick(self):
        self.DeleteCharacts(self.TblDC, 'DCcharacts')

    def ButDeleteACClick(self):
        self.DeleteCharacts(self.TblAC, 'ACcharacts')

    def DeleteCharacts(self, Tbl, DBTable):
        ids = self.GetTableSelectCol(Tbl)

        Msg = "Do you want to delete following Ids {} \n from {}".format(' , '.join([str(i) for i in ids]), DBTable)
        buttonReply = QMessageBox.question(self,
                                           'WARNING !!!!!!',
                                           Msg,
                                           QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if buttonReply == QMessageBox.Yes:       
            self.DB.DeleteCharact(DBTable, ids)
            self.UpdateSearchList(Devices=False, Types=False)
            self.ButSetDataClick()

    def TblACCellChanged(self, row, column):

        if self.ChkUpdateAC.isChecked():
            Val = self.TblAC.item(row, column).text()
            ids = int(self.TblAC.item(row, 0).text())
            self.DB.UpdateRow('ACcharacts',
                              {self.ACUpdateFields[column]: Val},
                              ('idACcharacts=', ids))
            self.DB.db.commit()

    def TblDCCellChanged(self, row, column):

        if self.ChkUpdateDC.isChecked():
            Val = self.TblDC.item(row, column).text()
            ids = int(self.TblDC.item(row, 0).text())
            self.DB.UpdateRow('DCcharacts',
                              {self.DCUpdateFields[column]: Val},
                              ('idDCcharacts=', ids))
            self.DB.db.commit()

    def TblTrtsCellChanged(self, row, column):

        if self.ChkUpdateTrts.isChecked():
            Val = self.TblTrts.item(row, column).text()
            ids = int(self.TblTrts.item(row,
                                        self.TrtsUpdateFields[column][3]).text())

#            print ('table', self.TrtsUpdateFields[column][0], 
#                  'Field', self.TrtsUpdateFields[column][1],
#                  'Val', Val,
#                  'Codition',self.TrtsUpdateFields[column][2],
#                  'Codition val', ids)

            self.DB.UpdateRow(Table=self.TrtsUpdateFields[column][0],
                              Fields={self.TrtsUpdateFields[column][1]: Val},
                              Condition=(self.TrtsUpdateFields[column][2], ids))
            self.DB.db.commit()

    def GetDataFromDb(self, AC=False, DC=False):
        if AC:
            ids = self.GetTableSelectCol(self.TblAC)
            Trts = self.GetTableSelectCol(self.TblAC,
                                          Col=self.ACFields['Trts.Name'][1],
                                          String=True)
            self.DataAC = self.DB.GetCharactFromId('ACcharacts', ids, Trts)

        if DC:
            ids = self.GetTableSelectCol(self.TblDC)
            Trts = self.GetTableSelectCol(self.TblDC,
                                          Col=self.DCFields['Trts.Name'][1],
                                          String=True)
            self.DataDC = self.DB.GetCharactFromId('DCcharacts', ids, Trts)

    def ButExportDCClick(self):
        self.GetDataFromDb(DC=True)

        fileName, _ = QFileDialog.getSaveFileName(self,
                                                  "Export Data", "",
                                                  "Pickle Files (*.pkl);;All Files (*)")
        ExptData = {}
        ExptData['DataDC'] = (self.DataDC)
        pickle.dump(ExptData,
                    open(fileName, 'wb'))

    def ButExportACClick(self):
        self.GetDataFromDb(AC=True)
        fileName, _ = QFileDialog.getSaveFileName(self,
                                                  "Export Data", "",
                                                  "Pickle Files (*.pkl);;All Files (*)")
        ExptData = {}
        ExptData['DataAC'] = (self.DataAC)
        pickle.dump(ExptData,
                    open(fileName, 'wb'))

    def ButViewACClick(self):
        self.GetDataFromDb(AC=True)
        Plot = PyFETpl.PyFETPlot()
        Plot.AddAxes(self.ViewAxsAC)
        Plot.PlotDataSet(self.DataAC, self.DataAC.keys(), PltIsOK=True)
        Plot.AddLegend(Axn='GMPoly')

    def ButAnalyzeACClick(self):
        self.GetDataFromDb(AC=True)
        self.DataExp = AppDataExp(self.DataAC,
                                  self.ChkCalcIrms.isChecked())
        self.DataExp.show()

    def ButViewDCClick(self):
        self.GetDataFromDb(DC=True)
        Plot = PyFETpl.PyFETPlot()
        Plot.AddAxes(self.ViewAxsDC)
        Plot.PlotDataSet(self.DataDC, self.DataDC.keys(), PltIsOK=True)
        Plot.AddLegend()

    def ButAnalyzeDCClick(self):
        self.GetDataFromDb(DC=True)
        self.DataExp = AppDataExp(self.DataDC, IsDC=True)
        self.DataExp.show()


class AppDataExp(QtWidgets.QMainWindow):

    def __init__(self, ACData, CalcIrmsNok=False, IsDC=False):
        QtWidgets.QMainWindow.__init__(self)
        uipath = os.path.join(os.path.dirname(__file__), 'GuiDataExplorer.ui')
        uic.loadUi(uipath, self)

        self.setWindowTitle('Data Explorer')
        self.Data = ACData

        self.LstTrt.itemSelectionChanged.connect(self.LstTrtChange)
        self.LstCy.itemSelectionChanged.connect(self.LstCyChange)

        self.LstVds.itemSelectionChanged.connect(self.LstVdsChange)
        self.LstVgs.itemSelectionChanged.connect(self.LstVgsChange)

        self.PltVsVgs.clicked.connect(self.PltVsVgsClick)

        self.ButPltVsX.clicked.connect(self.ButPltVsXClick)

        self.PlotFreq = None

        Vds = []
#        Vgs=[]
        for TrtN, Trt in sorted(self.Data.items()):
            self.LstTrt.addItem(TrtN)
            for cy in list(Trt.values()):
                for vd in cy['Vds']:
                    Vds.append(vd)
#                for vg in cy['Vgs']:Vgs.append(vg)
                if CalcIrmsNok:
                    print ('Calc Irms for', TrtN)
                    PyFETData.CalcNoiseIrmsCh(cy, IsOkFilt=False)

                if not cy['IsOK']:
                    continue
                if 'Irms' not in cy and not IsDC:
                    print ('Calc Irms for', TrtN)
                    PyFETData.CalcNoiseIrmsCh(cy)

        for vd in sorted(set(Vds)):
            self.CmbVds.addItem(str(vd))

#        self.SpinVgs.setMinimum(min(Vgs))
#        self.SpinVgs.setMaximum(max(Vgs))

    def LstTrtChange(self):
        sel = self.LstTrt.selectedItems()
        self.LstCy.clear()
        for Cy in self.Data[sel[0].text()]:
            self.LstCy.addItem(Cy)

        self.LstCy.setCurrentRow(0,QItemSelectionModel.Select)

    def LstCyChange(self):
        scy = self.LstCy.selectedItems()
        if len(scy) == 0:
            return

        strt = self.LstTrt.selectedItems()

        trt = strt[0].text()
        cy = scy[0].text()

        self.LstVds.clear()
        for iVd, Vd in enumerate(self.Data[trt][cy]['Vds']):
            self.LstVds.addItem('{} : {:.4f}'.format(iVd, Vd))
        self.LstVds.setCurrentRow(0,QItemSelectionModel.Select)

        self.LstVgs.clear()
        for iVg, Vg in enumerate(self.Data[trt][cy]['Vgs']):
            self.LstVgs.addItem('{} : {:.4f}'.format(iVg, Vg))
        self.LstVgs.setCurrentRow(0,QItemSelectionModel.Select)

        self.GenTrtInfo(trt, cy)

    def GenTrtInfo(self, trt, cy):

        Dat = self.Data[trt][cy]

        Info = 'Trt Name : ' + Dat['Name']
        for k, c in Dat['TrtTypes'].items():
            Info = '{}\n{} : {}'.format(Info, k, c)

        self.LblInfo.setText(Info)

    def ButPltVsXClick(self):

        self.PlotXX = PyFETpl.PyFETPlotParam()

        Axs = []
        for ck in self.GrpYY.findChildren(QtWidgets.QCheckBox):
            if ck.isChecked():
                Axs.append(ck.text())

        if len(Axs) == 0:
            return

        Bias = (None, self.SpinVgs.value())
        for Rad in self.GrpXVar.findChildren(QtWidgets.QRadioButton):
            if Rad.isChecked():
                xVar = Rad.text()

        PltUd0 = self.ChkVd0YY.isChecked()

        self.PlotXX.AddAxes(Axs, xVar)

        Trts = []
        for TrtN in self.LstTrt.selectedItems():
            Trts.append(TrtN.text())

        self.PlotXX.PlotDataSet(self.Data,
                                Trts,
                                xVar,
                                Bias,
                                PltUd0=PltUd0)

    def PltVsVgsClick(self):
        self.PlotVgs = PyFETpl.PyFETPlot()       
        
        Axs = []        
        for ck in self.GrpVgs.findChildren(QtWidgets.QCheckBox):    
            if ck.isChecked(): Axs.append(ck.text())
            
        self.PlotVgs.AddAxes(Axs)
        
        PltUd0 = self.ChkVd0.isChecked() 
        PltIsOk = self.ChkIsOK.isChecked() 
        
        for Rad in self.GrpVsVgsColor.findChildren(QtWidgets.QRadioButton):    
            if Rad.isChecked(): cPar = Rad.text()
        
        Trts=[]
        for TrtN in self.LstTrt.selectedItems():
            Trts.append(TrtN.text())

        self.PlotVgs.PlotDataSet(self.Data,Trts, 
                                 PltUd0=PltUd0, 
                                 PltIsOK=PltIsOk, 
                                 ColorOn=cPar)

        self.PlotVgs.AddLegend()
    
    def LstVdsChange(self):
        self.UpdatePltVsFreq()
    def LstVgsChange(self):
        self.UpdatePltVsFreq()    
    def UpdatePltVsFreq(self):
 
        if not self.ChkPltVsFreq.isChecked():return

        scy = self.LstCy.selectedItems()         
        strt = self.LstTrt.selectedItems()
        if len(scy)==0 or len(strt)==0:return
        trt = strt[0].text()
        cy = scy[0].text()   

        ivgs = [] 
        sel = self.LstVgs.selectedItems()
        for s in sel:
            st = s.text()
            ivgs.append(int(st.split(':')[0]))
            
        ivds = [] 
        sel = self.LstVds.selectedItems()
        for s in sel:
            st = s.text()
            ivds.append(int(st.split(':')[0]))

        if len(ivds)==0 or len(ivgs)==0:
            return

        if not self.PlotFreq:
            self.CreateNewPlotFreq()
        if not plt.fignum_exists(self.PlotFreq.Fig.number):
            del self.PlotFreq
            self.CreateNewPlotFreq()

        self.PlotFreq.ClearAxes()
        self.PlotFreq.Plot(self.Data[trt][cy], iVds=ivds, iVgs=ivgs,
                           ColorOnVgs=True)
        self.PlotFreq.Fig.canvas.draw()

    def CreateNewPlotFreq(self):
        print ('New plot')
        self.PlotFreq = PyFETpl.PyFETPlot()
        Axs = []
        for ck in self.GrpFreq.findChildren(QtWidgets.QCheckBox):
            if ck.isChecked():
                Axs.append(ck.text())
        self.PlotFreq.AddAxes(Axs)


def main():
    import argparse
    import pkg_resources

    # Add version option
    __version__ = pkg_resources.require("PyGFETdb")[0].version
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(
                            version=__version__))
    parser.parse_args()

    app = QtWidgets.QApplication(sys.argv)
    w = DBViewApp()
    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()



