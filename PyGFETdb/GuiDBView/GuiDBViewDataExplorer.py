

import os
import warnings

import numpy as np

from matplotlib import pyplot as plt
import matplotlib as mpl
from qtpy.QtWidgets import QFileDialog
from qtpy import QtWidgets, uic
from qtpy.QtCore import QSortFilterProxyModel
import seaborn as sns
import math
from PyGFETdb import DBInterface
from scipy import stats

from PyGFETdb.GuiDBView.GuiHelpers import FormatAxis, GenPSDBodeReport, PandasModel, LogPars


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
        self.ButRepPSDBode.clicked.connect(self.RepPSDBode_Click)

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

    def RepPSDBode_Click(self):
        dSel = self.GetSelection()
        GenPSDBodeReport(dfData=dSel)

