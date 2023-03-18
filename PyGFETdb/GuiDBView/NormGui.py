import os

from matplotlib import pyplot as plt
from qtpy.QtWidgets import QFileDialog
from qtpy import QtWidgets, uic
from qtpy.QtCore import QSortFilterProxyModel
from PyGFETdb import DBInterface
from PyGFETdb.GuiDBView.GuiHelpers import GenPSDBodeReport, PandasModel, GenScalarFigures, GenVectorFigures, \
    GenDeviceReportGui
import tempfile


class NormGui(QtWidgets.QMainWindow):
    def __init__(self, dfData):
        QtWidgets.QMainWindow.__init__(self)

        self.SelCol = None
        uipath = os.path.join(os.path.dirname(__file__), 'GuiNormalization.ui')
        uic.loadUi(uipath, self)

        self.dfData = dfData

        self.setWindowTitle('Normalization Dialog')

        self.ButGen.clicked.connect(self.ButGen_Click)

        includetypes = ('float', 'category', 'bool', 'datetime64[ns]', 'int')
        catCols = list(self.dfData.select_dtypes(include=includetypes).columns)
        self.LstColumn.addItems(catCols)

        self.LstColumn.itemSelectionChanged.connect(self.LstColumn_Changed)
        self.LstValue.itemSelectionChanged.connect(self.LstValue_Changed)

    def LstColumn_Changed(self):
        Sel = self.LstColumn.selectedItems()
        if len(Sel) == 0:
            print('Select Y var')
            return
        SelCol = [s.text() for s in Sel][0]

        self.SelCol = SelCol
        self.LstValue.clear()
        self.LstValue.addItems([str(s) for s in list(self.dfData[SelCol].unique())])

    def LstValue_Changed(self):
        Sel = self.LstValue.selectedItems()
        if len(Sel) == 0:
            print('Select Y var')
            return
        SelVal = [s.text() for s in Sel][0]

        if self.dfData[self.SelCol].dtype == int or self.dfData[self.SelCol].dtype == float:
            q = "{} == {}".format(self.SelCol, SelVal)
        else:
            q = "{} == '{}'".format(self.SelCol, SelVal)
        self.dRef = self.dfData.query(q)

        self.modelData = PandasModel(self.dRef.copy())
        self.proxyData = QSortFilterProxyModel()
        self.proxyData.setSourceModel(self.modelData)
        self.TblRefValues.setModel(self.proxyData)
        self.TblRefValues.show()

    def ButGen_Click(self):
        pass