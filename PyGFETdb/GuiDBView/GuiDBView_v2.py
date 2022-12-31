import argparse
import importlib
import multiprocessing
import sys
import os
import warnings
from multiprocessing import Pool

import numpy as np

from matplotlib import pyplot as plt
import matplotlib as mpl

from qtpy.QtWidgets import QHeaderView, QMessageBox
from qtpy.QtWidgets import QFileDialog, QAction, QInputDialog
from qtpy import QtWidgets, uic
from qtpy.QtCore import Qt, QSortFilterProxyModel, QItemSelectionModel
from PyGFETdb.DBCore2 import PyFETdb, Data2Pandas
from qtpy.QtCore import QAbstractTableModel, QModelIndex
import pandas as pd
import seaborn as sns
import math
from PyGFETdb.GuiDBView import UpdateDialogs
from PyGFETdb import DBInterface, __version__
import copy
from scipy import stats


class PandasModel(QAbstractTableModel):
    """A model to interface a Qt view with pandas dataframe """

    def __init__(self, dataframe: pd.DataFrame, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self._dataframe = dataframe

    def rowCount(self, parent=QModelIndex()) -> int:
        """ Override method from QAbstractTableModel

        Return row count of the pandas DataFrame
        """
        if parent == QModelIndex():
            return len(self._dataframe)
        return 0

    def columnCount(self, parent=QModelIndex()) -> int:
        """Override method from QAbstractTableModel

        Return column count of the pandas DataFrame
        """
        if parent == QModelIndex():
            return len(self._dataframe.columns)
        return 0

    def data(self, index: QModelIndex, role=Qt.ItemDataRole):
        """Override method from QAbstractTableModel

        Return data cell from the pandas DataFrame
        """
        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            data = self._dataframe.iloc[index.row(), index.column()]
            if type(data) == str:
                return data
            st = str(data)
            return st
        return None

    def headerData(self, section: int, orientation: Qt.Orientation, role: Qt.ItemDataRole):
        """Override method from QAbstractTableModel

        Return dataframe index as vertical header data and columns as horizontal header data.
        """
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._dataframe.columns[section])
            if orientation == Qt.Vertical:
                return str(self._dataframe.index[section])
        return None


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


LogPars = ('Irms', 'Vrms', 'NoA', 'NoC', 'IrmsNorm', 'VrmsNorm', 'NoANorm', 'NoCNorm')


def FormatAxis(PlotPars, dfAttrs):
    # figure generation
    nRows = math.ceil(np.sqrt(len(PlotPars)))
    nCols = math.ceil(len(PlotPars) / nRows)
    fig, Axs = plt.subplots(nrows=nRows, ncols=nCols)
    if nRows == 1 and nCols == 1:
        Axs = [Axs, ]
    else:
        Axs = Axs.flatten()

    # Format axis
    AxsDict = {}
    for ip, par in enumerate(PlotPars):
        ax = Axs[ip]
        AxsDict[par] = ax

        slabel = '{}[{}]'.format(par, dfAttrs['ColUnits'][par])
        ax.set_ylabel(slabel)

        # # select Vgs vector
        if (par in dfAttrs['ArrayColsNorm']) or par.endswith('Norm'):
            ax.set_xlabel('Vgs - CNP [V]')
        else:
            ax.set_xlabel('Vgs [V]')

        if (par in LogPars) or ('rms' in par):
            ax.set_yscale('log')
        elif par == 'PSD':
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlabel('Frequency [Hz]')
        elif par == 'BodeMag':
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlabel('Frequency [Hz]')
        elif par == 'BodePhase':
            ax.set_xscale('log')
            ax.set_xlabel('Frequency [Hz]')

    return fig, AxsDict


class DataExplorer(QtWidgets.QMainWindow):
    BoxPlotFunctions = {'Boxplot': sns.boxplot,
                        'violinplot': sns.violinplot,
                        'swarmplot': sns.swarmplot,
                        'boxenplot': sns.boxenplot,
                        'barplot': sns.barplot,
                        'stripplot': sns.stripplot,
                        }

    def __init__(self, dfRaw, ClassQueries, pdAttr):
        QtWidgets.QMainWindow.__init__(self)

        uipath = os.path.join(os.path.dirname(__file__), 'GuiDataExplorer_v2.ui')
        uic.loadUi(uipath, self)

        self.setWindowTitle('PyFETdb DataExplorer')

        self.dfDat = DBInterface.CalcElectricalParams(dfRaw,
                                                      ClassQueries,
                                                      pdAttr)

        self.modelData = PandasModel(self.dfDat.copy())
        self.proxyData = QSortFilterProxyModel()
        self.proxyData.setSourceModel(self.modelData)
        self.TblData.setModel(self.proxyData)
        self.TblData.show()

        self.LstBoxYPars.addItems(self.dfDat.attrs['ScalarCols'])

        includetypes = ('float', 'category', 'bool', 'datetime64[ns]')
        catCols = list(dfRaw.select_dtypes(include=includetypes).columns)
        self.CmbBoxX.addItems(catCols)
        self.CmbBoxX.setCurrentText('Device')
        self.CmbBoxHue.addItems(catCols)
        self.CmbBoxHue.setCurrentText('Device')
        self.CmbBoxType.addItems(self.BoxPlotFunctions.keys())

        self.LstLinesPars.addItems(self.dfDat.attrs['ArrayCols'] + ['PSD', 'BodeMag', 'BodePhase'])
        self.dfDat.attrs['ColUnits'].update({'PSD': 'A**2/Hz',
                                             'BodeMag': 'S',
                                             'BodePhase': 'Deg',
                                             })
        self.CmbLinesGroup.addItems(catCols)
        self.CmbLinesGroup.setCurrentText('Device')

        self.ButBoxPlot.clicked.connect(self.ButBoxPlot_Click)
        self.ButPlotVect.clicked.connect(self.ButPlotVect_Click)
        self.ButExportPkl.clicked.connect(self.ButExportPkl_Click)
        self.ButExportCSV.clicked.connect(self.ButExportCSV_Click)

    def GetSelection(self):
        Sel = self.TblData.selectedIndexes()
        rows = set([self.proxyData.mapToSource(s).row() for s in Sel])
        dSel = self.dfDat.loc[list(rows)]

        if not self.ChckQueries.isChecked():
            return dSel

        try:
            for q in self.TxtSelQueries.toPlainText().split('\n'):
                dSel.query(q, inplace=True)
        except:
            print("Error in query execution")

        return dSel

    def ButBoxPlot_Click(self):
        Sel = self.LstBoxYPars.selectedItems()
        if len(Sel) == 0:
            print('Select Y var')
            return
        PlotPars = [s.text() for s in Sel]

        dSel = self.GetSelection()

        nRows = math.ceil(np.sqrt(len(PlotPars)))
        nCols = math.ceil(len(PlotPars) / nRows)
        fig, Axs = plt.subplots(nrows=nRows, ncols=nCols)
        if nRows == 1 and nCols == 1:
            Axs = [Axs, ]
        else:
            Axs = Axs.flatten()

        PltFunct = self.BoxPlotFunctions[self.CmbBoxType.currentText()]

        for ic, p in enumerate(PlotPars):
            PltFunct(x=self.CmbBoxX.currentText(),
                     y=p,
                     hue=self.CmbBoxHue.currentText(),
                     data=dSel,
                     ax=Axs[ic])
            if p in LogPars:
                Axs[ic].set_yscale('log')

            plt.setp(Axs[ic].get_legend().get_texts(),
                     fontsize='xx-small')
            plt.setp(Axs[ic].get_legend().get_title(),
                     fontsize='xx-small')
            plt.setp(Axs[ic].get_xticklabels(),
                     rotation=45,
                     ha="right",
                     fontsize='xx-small',
                     rotation_mode="anchor")

        plt.show()

    def ButPlotVect_Click(self):
        Sel = self.LstLinesPars.selectedItems()
        if len(Sel) == 0:
            print('Select Y var')
            return
        PlotPars = [s.text() for s in Sel]

        dSel = self.GetSelection()
        GroupBy = self.CmbLinesGroup.currentText()
        dgroups = dSel.groupby(GroupBy, observed=True)

        # Colors iteration
        Norm = mpl.colors.Normalize(vmin=0, vmax=len(dgroups))
        cmap = mpl.cm.ScalarMappable(norm=Norm, cmap='jet')

        fig, AxsDict = FormatAxis(PlotPars, self.dfDat.attrs)

        for ic, gn in enumerate(dgroups.groups):
            gg = dgroups.get_group(gn)
            Col = cmap.to_rgba(ic)
            for parn in PlotPars:
                if parn in self.dfDat.attrs['ArrayColsNorm']:
                    xVar = self.dfDat.attrs['VgsNorm']
                else:
                    xVar = self.dfDat.attrs['Vgs']

                Vals = np.array([])
                ax = AxsDict[parn]
                for index, row in gg.iterrows():
                    if parn == 'PSD':
                        xVar = row.CharCl.GetFpsd()
                        if xVar is None:
                            continue
                        v = row.CharCl.GetPSD().transpose()
                        ax.loglog(xVar, row.CharCl.GetFitPSD(), 'k--', alpha=1)
                    elif parn == 'BodeMag':
                        xVar = row.CharCl.GetFgm()
                        if xVar is None:
                            continue
                        v = row.CharCl.GetGmMag().transpose()
                    elif parn == 'BodePhase':
                        xVar = row.CharCl.GetFgm()
                        if xVar is None:
                            continue
                        v = row.CharCl.GetGmPh().transpose()
                    else:
                        v = row[parn]

                    if type(v) == float:
                        continue
                    try:
                        if self.CheckLines.isChecked():
                            ax.plot(xVar, v.transpose(),
                                    color=Col,
                                    alpha=0.8,
                                    label=gn)
                        Vals = np.vstack((Vals, v)) if Vals.size else v
                    except:
                        pass

                if Vals.size == 0:
                    continue

                Vals = Vals.magnitude.transpose()
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    mean = np.nanmedian(Vals, axis=1)

                if self.CheckLinesMean.isChecked():
                    ax.plot(xVar, mean, '-.', color=Col, lw=1.5, label=gn)

                if self.CheckLinesStd.isChecked():
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        std = stats.tstd(Vals, axis=1)  # changed to mad for robustness
                    ax.fill_between(xVar, mean + std, mean - std, color=Col, alpha=0.2)

                handles, labels = ax.get_legend_handles_labels()
                by_label = dict(zip(labels, handles))
                ax.legend(by_label.values(), by_label.keys(), fontsize='xx-small')

        plt.show()

    def ButExportPkl_Click(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self, "Save PKl file", "", "All Files (*);; (*.pkl)", options=options)
        if fileName:
            dSel = self.GetSelection()
            dSel.to_pickle(fileName)

    def ButExportCSV_Click(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self, "Save CSV file", "", "All Files (*);; (*.csv)", options=options)
        if fileName:
            dSel = self.GetSelection()
            dSel.to_csv(fileName)


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
