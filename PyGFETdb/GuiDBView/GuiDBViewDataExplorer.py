import os


from matplotlib import pyplot as plt
from qtpy.QtWidgets import QFileDialog
from qtpy import QtWidgets, uic
from qtpy.QtCore import QSortFilterProxyModel
import seaborn as sns
from PyGFETdb import DBInterface
from PyGFETdb.GuiDBView.GuiHelpers import GenPSDBodeReport, PandasModel, GenScalarFigures, GenVectorFigures, GenDeviceReportGui
import tempfile


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

        includetypes = ('float', 'category', 'bool', 'datetime64[ns]', 'int')
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
        self.ButRepDevice.clicked.connect(self.ButRepDevice_Click)

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

        GenScalarFigures(Data=dSel,
                         PlotPars=PlotPars,
                         Xvar=self.CmbBoxX.currentText(),
                         Huevar=self.CmbBoxHue.currentText(),
                         PltFunct=self.BoxPlotFunctions[self.CmbBoxType.currentText()],
                         )
        plt.show()

    def ButPlotVect_Click(self):
        Sel = self.LstLinesPars.selectedItems()
        if len(Sel) == 0:
            print('Select Y var')
            return
        PlotPars = [s.text() for s in Sel]

        dSel = self.GetSelection()

        GenVectorFigures(data=dSel,
                         GroupBy=self.CmbLinesGroup.currentText(),
                         PlotPars=PlotPars,
                         Lines=self.CheckLines.isChecked(),
                         Mean=self.CheckLinesMean.isChecked(),
                         Std=self.CheckLinesStd.isChecked(),
                         )
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

    def ButRepDevice_Click(self):
        dSel = self.GetSelection()
        Sel = self.LstBoxYPars.selectedItems()
        if len(Sel) == 0:
            print('Select Y var')
            return
        ScalarPars = [s.text() for s in Sel]
        Sel = self.LstLinesPars.selectedItems()
        if len(Sel) == 0:
            print('Select Y var')
            return
        VectorPars = [s.text() for s in Sel]

        FileName = tempfile.NamedTemporaryFile(suffix='.pdf', delete=False)

        GenDeviceReportGui(Data=dSel,
                           FileOut=FileName.name,
                           ScalarPars=ScalarPars,
                           ScalarPltFunct=self.BoxPlotFunctions[self.CmbBoxType.currentText()],
                           VectorPars=VectorPars)

        # if platform.system() == 'Darwin':  # macOS
        #     subprocess.call(('open', FileName.name))
        # elif platform.system() == 'Windows':  # Windows
        #     os.startfile(FileName.name)
        # else:  # linux variants
        #     subprocess.call(('xdg-open', FileName.name))
