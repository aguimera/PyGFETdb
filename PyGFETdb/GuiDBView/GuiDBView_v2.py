
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
                     'Trts.Comments',
                     'Devices.Comments',
                     'Wafers.Comments',
                     'Devices.State',
                     'Devices.ExpOK',
                     'Devices.idDevices',
                     'Wafers.idWafers',
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
                       'DCcharacts.FileName')

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

        self.FillList(Lst=self.LstWafers,
                      Field='Name',
                      Data=self.DB.MultiSelect(Table='Wafers',
                                               Conditions={'Wafers.idWafers >': (0,)},
                                               Output=('Name',)))

        self.ConnectLst()

        self.ButSetData.clicked.connect(self.ButSetDataClick)
        self.ButViewDC.clicked.connect(self.ButViewDCClick)
        self.ButViewAC.clicked.connect(self.ButViewACClick)

    def ConnectLst(self):
        self.LstWafers.itemSelectionChanged.connect(self.LstWafersChange)
        self.LstDevices.itemSelectionChanged.connect(self.LstDevicesChange)

    def FillList(self, Lst, Field, Data):
        Vals = []
        for r in Data:
            Vals.append(r[Field])
        Lst.clear()

        for d in sorted(set(Vals)):
            Lst.addItem(str(d))

    def LstWafersChange(self):
        sel = self.LstWafers.selectedItems()
        if len(sel) == 0:
            return

        Conditions = {'Wafers.Name=': []}
        for s in sel:
            Conditions['Wafers.Name='].append(s.text())

        self.FillList(Lst=self.LstDevices,
                      Field='Devices_Name',
                      Data=self.DB.GetTrtsInfo(Conditions=Conditions,
                                               Output=('Devices.Name',)))

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

        self.dfDCchars = pd.DataFrame(self.DB.GetCharactInfo(Table='DCcharacts',
                                                             Conditions={'Trts.idTrts =': idtrts},
                                                             Output=self.DCCharTableCols))


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

        # if dic['Ylog']:
        #     ax.set_yscale('log')
        #
        # # select Vgs vector
        # if par.endswith('Norm'):
        #     dic['xVar'] = dfAttrs['VgsNorm']
        #     ax.set_xlabel('Vgs - CNP [V]')
        # else:
        #     ax.set_xlabel('Vgs [V]')
        #     dic['xVar'] = dfAttrs['Vgs']

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
        dgroups = dSel.groupby(GroupBy)

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
