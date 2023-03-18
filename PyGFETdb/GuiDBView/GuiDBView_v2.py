
import argparse
import multiprocessing
import sys
import os
from multiprocessing import Pool
from qtpy.QtWidgets import QMessageBox
from qtpy import QtWidgets, uic
from qtpy.QtCore import QSortFilterProxyModel
from PyGFETdb.DBCore2 import PyFETdb, Data2Pandas
import pandas as pd
from PyGFETdb.GuiDBView import UpdateDialogs
from PyGFETdb import DBInterface, __version__
import copy
from PyGFETdb.GuiDBView.GuiDBViewDataExplorer import DataExplorer
from PyGFETdb.GuiDBView.GuiHelpers import PandasModel, AddCycle, AddColRow


class DBViewApp(QtWidgets.QMainWindow):
    WafersColumns = ('Name',
                     'idWafers',
                     'Substrate',
                     'Run',
                     'Masks',
                     'Comments')

    DevicesTableCols = ('Wafers.Name',
                        'Wafers.Run',
                        'Wafers.Masks',
                        'Wafers.Graphene',
                        'Wafers.Comments',
                        'Devices.Name',
                        'Devices.State',
                        'Devices.ExpOK',
                        'Wafers.Comments',
                        'Devices.idDevices',
                        'Wafers.idWafers',
                        )

    TrtsTableCols = ('Trts.idTrts',
                     'Trts.Name',
                     'VTrts.DCMeas',
                     'VTrts.ACMeas',
                     'VTrts.GMeas',
                     'TrtTypes.Name',
                     'TrtTypes.Length',
                     'TrtTypes.Width',
                     'TrtTypes.Contact',
                     'TrtTypes.Pass',
                     'TrtTypes.Area',
                     'Wafers.Comments',
                     'Devices.Comments',
                     'Wafers.Comments',
                     'Wafers.Run',
                     'Devices.idDevices',
                     'Wafers.idWafers',
                     'TrtTypes.idTrtTypes',
                     'Trts.Comments',
                     'Devices.State',
                     'Devices.ExpOK',
                     )

    DCCharTableCols = ('DCcharacts.idDCcharacts',
                       'Trts.Name',
                       'DCcharacts.MeasDate',
                       'DCcharacts.IsOK',
                       'DCcharacts.IsCmp',
                       'DCcharacts.Ph',
                       'DCcharacts.IonStrength',
                       'DCcharacts.FuncStep',
                       'DCcharacts.AnalyteCon',
                       'DCcharacts.Comments',
                       'DCcharacts.UpdateDate',
                       'DCcharacts.FileName',
                       )

    ACCharTableCols = ('ACcharacts.idACcharacts',
                       'Trts.Name',
                       'ACcharacts.MeasDate',
                       'ACcharacts.IsOK',
                       'ACcharacts.IsCmp',
                       'ACcharacts.Ph',
                       'ACcharacts.IonStrength',
                       'ACcharacts.FuncStep',
                       'ACcharacts.AnalyteCon',
                       'ACcharacts.Comments',
                       'ACcharacts.UpdateDate',
                       'ACcharacts.DC_id',
                       'ACcharacts.FileName')

    ClassQueriesDC = copy.deepcopy(DBInterface.ClassQueriesDC)
    pdAttrDC = copy.deepcopy(DBInterface.pdAttrDC)
    ClassQueries = copy.deepcopy(DBInterface.ClassQueries)
    pdAttr = copy.deepcopy(DBInterface.pdAttr)

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)

        uipath = os.path.join(os.path.dirname(__file__), 'GuiDBView_v2.ui')
        uic.loadUi(uipath, self)

        self.setWindowTitle('PyFETdb Viewer v' + __version__)

        keypath = os.path.join(os.path.dirname(__file__), 'key.key')

        self.DB = PyFETdb(open(keypath, 'rb').read())
        self.DB._DEBUG = False

        self.InitLists()

        users = self.DB.MultiSelect(Table='Users',
                                    Conditions={'Users.idUsers >': (0,), },
                                    Output=('Name',))
        self.LstUsers.addItem('None')
        self.LstUsers.addItems(pd.DataFrame(users)['Name'])

        self.ButSetData.clicked.connect(self.ButSetData_Click)
        self.ButViewDC.clicked.connect(self.ButViewDC_Click)
        self.ButViewAC.clicked.connect(self.ButViewAC_Click)
        self.ButEditWafer.clicked.connect(self.ButEditWafer_Click)
        self.ButEditDevice.clicked.connect(self.ButEditDevice_Click)
        self.ButEditTrts.clicked.connect(self.ButEditTrts_Click)
        self.ButEditTrtsType.clicked.connect(self.ButEditTrtsType_Click)
        self.ButEditDC.clicked.connect(self.ButEditDC_Click)
        self.ButEditAC.clicked.connect(self.ButEditAC_Click)
        self.ButEditElecParsDC.clicked.connect(self.ButEditEleDC_Click)
        self.ButEditElecParsAC.clicked.connect(self.ButEditEleAC_Click)

        self.LstWafers.itemSelectionChanged.connect(self.LstWafersChange)
        self.LstDevices.itemSelectionChanged.connect(self.LstDevicesChange)
        self.LstSubstrates.itemSelectionChanged.connect(self.LstChange)
        self.LstMasks.itemSelectionChanged.connect(self.LstChange)
        self.LstRuns.itemSelectionChanged.connect(self.LstChange)
        self.LstUsers.itemSelectionChanged.connect(self.InitLists)

    def InitLists(self):
        self.LstSubstrates.clear()
        self.LstMasks.clear()
        self.LstRuns.clear()
        self.LstWafers.clear()

        sel = self.LstUsers.selectedItems()
        ucond = {'Users.Name =': []}
        for s in sel:
            if s.text() == 'None':
                break
            ucond['Users.Name ='].append(s.text())
        if len(ucond['Users.Name =']):
            data = self.DB.GetCharactInfo(Table='DCcharacts',
                                          Conditions=ucond,
                                          Output=('Wafers.Name',))
            wcond = {'Wafers.Name =': list(pd.DataFrame(data)['Wafers_Name'])}
        else:
            wcond = {'Wafers.idWafers >': (0,), }

        data = self.DB.MultiSelect(Table='Wafers',
                                   Conditions=wcond,
                                   Output=self.WafersColumns)
        self.WafersDF = pd.DataFrame(data)

        self.LstSubstrates.addItem('None')
        self.LstSubstrates.addItems(sorted(list(self.WafersDF.Substrate.unique())))
        self.LstMasks.addItem('None')
        self.LstMasks.addItems(sorted(list(self.WafersDF.Masks.unique())))
        self.LstRuns.addItem('None')
        self.LstRuns.addItems(sorted(list(self.WafersDF.Run.unique())))
        self.LstWafers.addItems(sorted(list(self.WafersDF.Name.unique())))

    def LstChange(self):
        df = self.WafersDF.copy()
        sel = self.LstSubstrates.selectedItems()
        qs = []
        for s in sel:
            if s.text() == 'None':
                break
            qs.append("Substrate == '{}'".format(s.text()))
        if len(qs):
            df.query(' | '.join(qs), inplace=True)

        sel = self.LstMasks.selectedItems()
        qs = []
        for s in sel:
            if s.text() == 'None':
                break
            qs.append("Masks == '{}'".format(s.text()))
        if len(qs):
            df.query(' | '.join(qs), inplace=True)

        sel = self.LstRuns.selectedItems()
        qs = []
        for s in sel:
            if s.text() == 'None':
                break
            qs.append("Run == '{}'".format(s.text()))
        if len(qs):
            df.query(' | '.join(qs), inplace=True)

        self.LstWafers.clear()
        self.LstWafers.addItems(df.Name.unique())

    def LstWafersChange(self):
        sel = self.LstWafers.selectedItems()
        if len(sel) == 0:
            return

        Conditions = {'Wafers.Name=': []}
        for s in sel:
            Conditions['Wafers.Name='].append(s.text())
        data = self.DB.GetTrtsInfo(Conditions=Conditions,
                                   Output=('Devices.Name',))
        data = pd.DataFrame(data)
        self.LstDevices.clear()
        self.LstDevices.addItems(data['Devices_Name'].unique())

    def LstDevicesChange(self):
        sel = self.LstDevices.selectedItems()
        if len(sel) == 0:
            return
        Conditions = {'Devices.Name=': []}
        for s in sel:
            Conditions['Devices.Name='].append(s.text())

        self.dfDevs = pd.DataFrame(self.DB.GetTrtsInfo(Conditions=Conditions,
                                                       Output=self.DevicesTableCols))

        self.modelDevs = PandasModel(self.dfDevs.copy())
        self.proxyDevs = QSortFilterProxyModel()
        self.proxyDevs.setSourceModel(self.modelDevs)
        self.TblDevices.setModel(self.proxyDevs)
        self.TblDevices.show()

        self.dfTrts = pd.DataFrame(self.DB.GetTrtsInfo(Conditions=Conditions,
                                                       Output=self.TrtsTableCols))

        self.modelTrts = PandasModel(self.dfTrts.copy())
        self.proxyTrts = QSortFilterProxyModel()
        self.proxyTrts.setSourceModel(self.modelTrts)
        self.TblTrts.setModel(self.proxyTrts)
        self.TblTrts.show()

    def ButSetData_Click(self):
        Sel = self.TblTrts.selectedIndexes()
        rows = set([self.proxyTrts.mapToSource(s).row() for s in Sel])
        idtrts = list(self.dfTrts.loc[list(rows)]['Trts_idTrts'])

        data = self.DB.GetCharactInfo(Table='DCcharacts',
                                      Conditions={'Trts.idTrts =': idtrts},
                                      Output=self.DCCharTableCols)
        self.dfDCchars = pd.DataFrame(data)

        self.modelDCchars = PandasModel(self.dfDCchars.copy())
        self.proxyDCChars = QSortFilterProxyModel()
        self.proxyDCChars.setSourceModel(self.modelDCchars)
        self.TblDC.setModel(self.proxyDCChars)
        self.TblDC.show()

        self.dfACchars = pd.DataFrame(self.DB.GetCharactInfo(Table='ACcharacts',
                                                             Conditions={'Trts.idTrts =': idtrts},
                                                             Output=self.ACCharTableCols))

        self.modelACchars = PandasModel(self.dfACchars.copy())
        self.proxyACChars = QSortFilterProxyModel()
        self.proxyACChars.setSourceModel(self.modelACchars)
        self.TblAC.setModel(self.proxyACChars)
        self.TblAC.show()

    def GetData(self, Table, idchars):
        Args = [(Table, Id, True) for Id in idchars]
        print('Get records ', len(Args))
        with Pool(8) as p:
            Data = p.starmap(self.DB.GetCharactFromId, Args)

        return Data


    def ButViewDC_Click(self):
        Sel = self.TblDC.selectedIndexes()
        rows = set([self.proxyDCChars.mapToSource(s).row() for s in Sel])
        idchars = list(self.dfDCchars.loc[list(rows)]['DCcharacts_idDCcharacts'])

        Data = self.GetData(Table='DCcharacts', idchars=idchars)
        dfRaw = Data2Pandas(Data)
        AddCycle(dfRaw)

        self.DataExp = DataExplorer(dfRaw,
                                    ClassQueries=self.ClassQueriesDC,
                                    pdAttr=self.pdAttrDC)
        self.DataExp.show()

    def ButViewAC_Click(self):
        Sel = self.TblAC.selectedIndexes()
        rows = set([self.proxyACChars.mapToSource(s).row() for s in Sel])
        idchars = list(self.dfACchars.loc[list(rows)]['ACcharacts_idACcharacts'])

        Data = self.GetData(Table='ACcharacts', idchars=idchars)
        dfRaw = Data2Pandas(Data)
        AddCycle(dfRaw)

        self.DataExp = DataExplorer(dfRaw,
                                    ClassQueries=self.ClassQueries,
                                    pdAttr=self.pdAttr)
        self.DataExp.show()

    def ButEditWafer_Click(self):
        sel = self.LstWafers.selectedItems()
        if len(sel) == 0:
            return

        Columns = UpdateDialogs.EditWaferColumns.copy()
        cond = {'Wafers.Name = ': []}
        for s in sel:
            cond['Wafers.Name = '].append(s.text())

        self.ShowUpdateDialog(Columns=Columns,
                              Table='Wafers',
                              Records=self.DB.MultiSelect(Table='Wafers',
                                                          Conditions=cond,
                                                          Output=Columns.keys()))

    def ButEditDevice_Click(self):
        sel = self.LstDevices.selectedItems()
        if len(sel) == 0:
            return

        Columns = UpdateDialogs.EditDeviceColumns.copy()
        cond = {'Devices.Name = ': []}
        for s in sel:
            cond['Devices.Name = '].append(s.text())

        self.ShowUpdateDialog(Columns=Columns,
                              Table='Devices',
                              Records=self.DB.MultiSelect(Table='Devices',
                                                          Conditions=cond,
                                                          Output=Columns.keys()))

    def ShowUpdateDialog(self, Table, Columns, Records):
        records = []
        for r in Records:
            c = copy.deepcopy(Columns)
            for n, par in c.items():
                try:
                    if par['type'] == 'int':
                        par['value'] = int(r[n])
                    elif par['type'] == 'float':
                        par['value'] = float(r[n])
                    else:
                        par['value'] = r[n]
                except ValueError:
                    par['value'] = 0
            records.append(c)

        self.UpdateWind = UpdateDialogs.TableEditorWindow(Table=Table,
                                                          Records=records,
                                                          DB=self.DB)
        self.UpdateWind.show()

    def ButEditTrts_Click(self):
        Sel = self.TblTrts.selectedIndexes()
        rows = set([self.proxyTrts.mapToSource(s).row() for s in Sel])
        idtrts = list(self.dfTrts.loc[list(rows)]['Trts_idTrts'])

        Columns = UpdateDialogs.EditTrtColumns.copy()
        cond = {'Trts.idTrts = ': []}
        for s in idtrts:
            cond['Trts.idTrts = '].append(s)

        self.ShowUpdateDialog(Columns=Columns,
                              Table='Trts',
                              Records=self.DB.MultiSelect(Table='Trts',
                                                          Conditions=cond,
                                                          Output=Columns.keys()))

    def ButEditTrtsType_Click(self):
        Msg = "WARNING!!!! This action can affect other Transistors \n\n do you want to continue??"
        buttonReply = QMessageBox.question(self,
                                           'WARNING !!!!!!',
                                           Msg,
                                           QMessageBox.Yes | QMessageBox.No,
                                           QMessageBox.No)
        if buttonReply == QMessageBox.No:
            return

        Sel = self.TblTrts.selectedIndexes()
        rows = set([self.proxyTrts.mapToSource(s).row() for s in Sel])
        types = list(self.dfTrts.loc[list(rows)]['TrtTypes_Name'])

        Columns = UpdateDialogs.EditTrtTypesColumns.copy()
        cond = {'TrtTypes.Name = ': []}
        for s in set(types):
            cond['TrtTypes.Name = '].append(s)

        self.ShowUpdateDialog(Columns=Columns,
                              Table='TrtTypes',
                              Records=self.DB.MultiSelect(Table='TrtTypes',
                                                          Conditions=cond,
                                                          Output=Columns.keys()))

    def ButEditDC_Click(self):
        Sel = self.TblDC.selectedIndexes()
        rows = set([self.proxyDCChars.mapToSource(s).row() for s in Sel])
        idchars = list(self.dfDCchars.loc[list(rows)]['DCcharacts_idDCcharacts'])

        Columns = UpdateDialogs.EditDCCharColumns.copy()
        cond = {'DCcharacts.idDCcharacts = ': []}
        for s in set(idchars):
            cond['DCcharacts.idDCcharacts = '].append(s)

        self.ShowUpdateDialog(Columns=Columns,
                              Table='DCcharacts',
                              Records=self.DB.MultiSelect(Table='DCcharacts',
                                                          Conditions=cond,
                                                          Output=Columns.keys()))

    def ButEditAC_Click(self):
        Sel = self.TblAC.selectedIndexes()
        rows = set([self.proxyACChars.mapToSource(s).row() for s in Sel])
        idchars = list(self.dfACchars.loc[list(rows)]['ACcharacts_idACcharacts'])

        Columns = UpdateDialogs.EditACCharColumns.copy()
        cond = {'ACcharacts.idACcharacts = ': []}
        for s in set(idchars):
            cond['ACcharacts.idACcharacts = '].append(s)

        self.ShowUpdateDialog(Columns=Columns,
                              Table='ACcharacts',
                              Records=self.DB.MultiSelect(Table='ACcharacts',
                                                          Conditions=cond,
                                                          Output=Columns.keys()))

    def ButEditEleDC_Click(self):
        self.UpdateWind = UpdateDialogs.ElectricalParamsEditor(ClassQueries=self.ClassQueriesDC,
                                                               pdAttr=self.pdAttrDC)
        self.UpdateWind.exec_()
        self.ClassQueriesDC = copy.deepcopy(self.UpdateWind.ClassQueries)
        self.pdAttrDC = copy.deepcopy(self.UpdateWind.pdAttr)

    def ButEditEleAC_Click(self):
        self.UpdateWind = UpdateDialogs.ElectricalParamsEditor(ClassQueries=self.ClassQueries,
                                                               pdAttr=self.pdAttr)
        self.UpdateWind.exec_()
        self.ClassQueries = copy.deepcopy(self.UpdateWind.ClassQueries)
        self.pdAttr = copy.deepcopy(self.UpdateWind.pdAttr)


def main():
    # Add version option
    # __version__ = pkg_resources.require("PyGFETdb")[0].version
    parser = argparse.ArgumentParser()
    # parser.add_argument('--version', action='version',
    #                     version='%(prog)s {version}'.format(
    #                         version=__version__))
    parser.parse_args()

    app = QtWidgets.QApplication(sys.argv)
    w = DBViewApp()
    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    if sys.platform.startswith('win'):
        # On Windows calling this function is necessary.
        multiprocessing.freeze_support()

    main()
