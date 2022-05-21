import argparse
import sys
import os

import numpy as np
from matplotlib import pyplot as plt
from qtpy.QtWidgets import QHeaderView, QMessageBox
from qtpy.QtWidgets import QFileDialog, QAction, QInputDialog
from qtpy import QtWidgets, uic
from qtpy.QtCore import Qt, QSortFilterProxyModel, QItemSelectionModel
from PyGFETdb.DBCore2 import PyFETdb, Data2Pandas
from PyGFETdb import DBInterface
from qtpy.QtCore import QAbstractTableModel, QModelIndex
import pandas as pd
import seaborn as sns
import math
import matplotlib as mpl
from PyGFETdb.GuiDBView import UpdateDialogs


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

    # def sort(self, Ncol, order):
    #     """Sort table by given column number.
    #     """
    #     try:
    #         self.layoutAboutToBeChanged.emit()
    #         self._dataframe = self._dataframe.sort_values(self._dataframe.columns[Ncol],
    #                                                       ascending=not order)
    #         self.layoutChanged.emit()
    #     except Exception as e:
    #         print(e)


class DBViewApp(QtWidgets.QMainWindow):
    WafersColumns = ('Name',
                     'idWafers',
                     'Substrate',
                     'Run',
                     'Masks',
                     'Comments')

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

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)

        uipath = os.path.join(os.path.dirname(__file__), 'GuiDBView_v2.ui')
        uic.loadUi(uipath, self)

        self.setWindowTitle('PyFETdb Viewer')
        self.DB = PyFETdb(open('key.key', 'rb').read())
        self.DB._DEBUG = False

        data = self.DB.MultiSelect(Table='Wafers',
                                   Conditions={'Wafers.idWafers >': (0,), },
                                   Output=self.WafersColumns)

        self.WafersDF = pd.DataFrame(data)

        self.LstSubstrates.addItem('None')
        self.LstSubstrates.addItems(sorted(list(self.WafersDF.Substrate.unique())))
        self.LstMasks.addItem('None')
        self.LstMasks.addItems(sorted(list(self.WafersDF.Masks.unique())))
        self.LstRuns.addItem('None')
        self.LstRuns.addItems(sorted(list(self.WafersDF.Run.unique())))
        self.LstWafers.addItems(sorted(list(self.WafersDF.Name.unique())))

        self.ButSetData.clicked.connect(self.ButSetDataClick)
        self.ButViewDC.clicked.connect(self.ButViewDCClick)
        self.ButViewAC.clicked.connect(self.ButViewACClick)
        self.ButEditWafer.clicked.connect(self.ButEditWaferlick)

        self.LstWafers.itemSelectionChanged.connect(self.LstWafersChange)
        self.LstDevices.itemSelectionChanged.connect(self.LstDevicesChange)
        self.LstSubstrates.itemSelectionChanged.connect(self.LstChange)
        self.LstMasks.itemSelectionChanged.connect(self.LstChange)
        self.LstRuns.itemSelectionChanged.connect(self.LstChange)

    def LstChange(self):
        df = self.WafersDF.copy()

        sel = self.LstSubstrates.selectedItems()
        for s in sel:
            if s.text() == 'None':
                break
            df.query("Substrate == '{}'".format(s.text()), inplace=True)

        sel = self.LstMasks.selectedItems()
        MasksNone = False
        for s in sel:
            if s.text() == 'None':
                break
            df.query("Masks == '{}'".format(s.text()), inplace=True)

        sel = self.LstRuns.selectedItems()
        for s in sel:
            if s.text() == 'None':
                break
            df.query("Run == '{}'".format(s.text()), inplace=True)

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

        self.dfTrts = pd.DataFrame(self.DB.GetTrtsInfo(Conditions=Conditions,
                                                       Output=self.TrtsTableCols))

        self.modelTrts = PandasModel(self.dfTrts.copy())
        self.proxyTrts = QSortFilterProxyModel()
        self.proxyTrts.setSourceModel(self.modelTrts)
        self.TblTrts.setModel(self.proxyTrts)
        self.TblTrts.show()

    def ButSetDataClick(self):
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

    def ButViewDCClick(self):
        Sel = self.TblDC.selectedIndexes()
        rows = set([self.proxyDCChars.mapToSource(s).row() for s in Sel])
        idchars = list(self.dfDCchars.loc[list(rows)]['DCcharacts_idDCcharacts'])

        Data = []
        for Id in idchars:
            Data.append(self.DB.GetCharactFromId(Table='DCcharacts', Id=Id))
        dfRaw = Data2Pandas(Data)

        self.DataExp = DataExplorer(dfRaw, Paramters='DC')
        self.DataExp.show()

    def ButViewACClick(self):
        Sel = self.TblAC.selectedIndexes()
        rows = set([self.proxyACChars.mapToSource(s).row() for s in Sel])
        idchars = list(self.dfACchars.loc[list(rows)]['ACcharacts_idACcharacts'])

        Data = []
        for Id in idchars:
            Data.append(self.DB.GetCharactFromId(Table='ACcharacts', Id=Id))
        dfRaw = Data2Pandas(Data)

        self.DataExp = DataExplorer(dfRaw, Paramters='AC')
        self.DataExp.show()

    def ButEditWaferlick(self):
        sel = self.LstWafers.selectedItems()
        if len(sel) == 0:
            return

        Columns = UpdateDialogs.EditWaferColumns.copy()

        WafName = sel[0].text()
        res = self.DB.MultiSelect(Table='Wafers',
                                  Conditions={'Wafers.Name =': (WafName,), },
                                  Output=Columns.keys())

        for n, par in Columns.items():
            par['value'] = res[0][n]

        self.UpdateWind = UpdateDialogs.TableEditorWindow('Wafers', Columns, self.DB)
        self.UpdateWind.show()


LogPars = ('Irms', 'Vrms', 'NoA', 'NoC')


def FormatAxis(PlotPars, dfAttrs):
    # figure generation
    nRows = math.ceil(np.sqrt(len(PlotPars)))
    nCols = math.ceil(len(PlotPars) / nRows)
    fig, Axs = plt.subplots(nrows=nRows, ncols=nCols)
    Axs = Axs.flatten()

    # Format axis
    AxsDict = {}
    for ip, par in enumerate(PlotPars):
        ax = Axs[ip]
        AxsDict[par] = ax

        slabel = '{}[{}]'.format(par, dfAttrs['ColUnits'][par])
        ax.set_ylabel(slabel)

        if par in LogPars:
            ax.set_yscale('log')

        # # select Vgs vector
        if par.endswith('Norm'):
            ax.set_xlabel('Vgs - CNP [V]')
        else:
            ax.set_xlabel('Vgs [V]')

    return fig, AxsDict


class DataExplorer(QtWidgets.QMainWindow):
    def __init__(self, dfRaw, Paramters='DC'):
        QtWidgets.QMainWindow.__init__(self)
        uipath = os.path.join(os.path.dirname(__file__), 'GuiDataExplorer_v2.ui')
        uic.loadUi(uipath, self)

        if Paramters == 'DC':
            self.dfDat = DBInterface.CalcElectricalParams(dfRaw,
                                                          DBInterface.ClassQueriesDC,
                                                          DBInterface.pdAttrDC)
        elif Paramters == 'AC':
            self.dfDat = DBInterface.CalcElectricalParams(dfRaw,
                                                          DBInterface.ClassQueries,
                                                          DBInterface.pdAttr)

        self.modelData = PandasModel(self.dfDat.copy())
        self.proxyData = QSortFilterProxyModel()
        self.proxyData.setSourceModel(self.modelData)
        self.TblData.setModel(self.proxyData)
        self.TblData.show()

        catCols = list(dfRaw.select_dtypes(include=('category', 'bool', 'datetime64[ns]')).columns)
        self.CmbBoxX.addItems(catCols)
        self.LstBoxYPars.addItems(self.dfDat.attrs['ScalarCols'])
        self.CmbBoxHue.addItems(catCols)
        self.CmbBoxType.addItems(['boxplot', 'swarmplot'])

        self.LstLinesPars.addItems(self.dfDat.attrs['ArrayCols'])
        self.CmbLinesGroup.addItems(catCols)

        self.ButBoxPlot.clicked.connect(self.ButBoxPlotClick)
        self.ButPlotVect.clicked.connect(self.ButPlotVectClick)

    def ButBoxPlotClick(self):
        PlotPars = self.LstBoxYPars.selectedItems()
        if len(PlotPars) == 0:
            print('Select Y var')
            return

        Sel = self.TblData.selectedIndexes()
        rows = set([self.proxyData.mapToSource(s).row() for s in Sel])
        dSel = self.dfDat.loc[list(rows)]

        nRows = math.ceil(np.sqrt(len(PlotPars)))
        nCols = math.ceil(len(PlotPars) / nRows)
        fig, Axs = plt.subplots(nrows=nRows, ncols=nCols)
        Axs = Axs.flatten()

        for ic, p in enumerate(PlotPars):
            sns.boxplot(x=self.CmbBoxX.currentText(),
                        y=p.text(),
                        hue=self.CmbBoxHue.currentText(),
                        data=dSel,
                        ax=Axs[ic])

    def ButPlotVectClick(self):
        Sel = self.LstLinesPars.selectedItems()
        if len(Sel) == 0:
            print('Select Y var')
            return
        PlotPars = [s.text() for s in Sel]

        Sel = self.TblData.selectedIndexes()
        rows = set([self.proxyData.mapToSource(s).row() for s in Sel])
        dSel = self.dfDat.loc[list(rows)]

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
                if parn.endswith('Norm'):
                    xVar = self.dfDat.attrs['VgsNorm']
                else:
                    xVar = self.dfDat.attrs['Vgs']

                Vals = np.array([])
                ax = AxsDict[parn]
                for index, row in gg.iterrows():
                    v = row[parn]
                    if type(v) == float:
                        continue
                    try:
                        ax.plot(xVar, v.transpose(), color=Col, alpha=0.8)
                        Vals = np.vstack((Vals, v)) if Vals.size else v
                    except:
                        pass

                # Vals = Vals.magnitude.transpose()
                # mean = np.nanmedian(Vals, axis=1)
                # std = stats.tstd(Vals, axis=1)  # changed to mad for robustness
                #
                # ax.plot(xVar, mean, color='k', lw=1.5, label=gn)
                # if not dic['Ylog']:
                #     ax.fill_between(xVar, mean + std, mean - std, color='k', alpha=0.2)

            # AxsDict['GM'].legend(fontsize='xx-small')
            # fig.tight_layout()
            # PDF.savefig(fig)
            # plt.close(fig)


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
    main()
