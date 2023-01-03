import glob
import math
import os
import platform
import subprocess
import tempfile

import matplotlib as mpl
import numpy as np
import pandas as pd
from PyPDF2 import PdfMerger
from PyQt5.QtCore import QAbstractTableModel, QModelIndex
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from qtpy.QtCore import Qt
from tqdm.contrib.concurrent import process_map


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


def GenPSDBodeFigure(Args):
    Cl = Args['Cl']
    tmpDir = Args['tmpDir']
    dfAttrs = Args['dfAttrs']
    Vgs = Cl.GetVgs()
    Norm = mpl.colors.Normalize(vmin=0, vmax=len(Vgs))
    cmap = mpl.cm.ScalarMappable(norm=Norm, cmap='jet')

    fig, AxsDict = FormatAxis(('Ids', 'BodeMag', 'PSD', 'BodePhase'), dfAttrs)
    fig.set_size_inches(11, 8)
    fig.suptitle(Args['title'])

    ax = AxsDict['Ids']
    ax.plot(Vgs, Cl.GetIds().flatten())
    fpsd = Cl.GetFpsd()
    if fpsd is not None:
        psd = Cl.GetPSD()
        ax = AxsDict['PSD']
        for iv, vg in enumerate(Vgs):
            ax.loglog(fpsd, psd[:, iv],
                      color=cmap.to_rgba(iv),
                      label='Vgs = {:0.2}'.format(vg))
        ax.legend()

    Fgm = Cl.GetFgm()
    if Fgm is not None:
        GmMag = Cl.GetGmMag()
        ax = AxsDict['BodeMag']
        for iv, vg in enumerate(Vgs):
            ax.loglog(Fgm, GmMag[:, iv],
                      color=cmap.to_rgba(iv),
                      label='Vgs = {:0.2}'.format(vg))
            Gm = np.abs(Cl.GetGM(Vgs=vg).flatten())
            ax.plot(Fgm[0], Gm, color=cmap.to_rgba(iv), marker='*')
        ax.legend()
        GmPh = Cl.GetGmPh()
        ax = AxsDict['BodePhase']
        for iv, vg in enumerate(Vgs):
            ax.semilogx(Fgm, GmPh[:, iv],
                        color=cmap.to_rgba(iv),
                        label='Vgs = {:0.2}'.format(vg))
        ax.legend()
    fig.tight_layout()
    file = tempfile.NamedTemporaryFile(dir=tmpDir.name, suffix='.pdf', delete=False)
    fig.savefig(file, format='pdf')
    plt.close(fig)


def GenPSDBodeReport(dfData):
    dfData.attrs['ColUnits'].update({'PSD': 'A**2/Hz',
                                     'BodeMag': 'S',
                                     'BodePhase': 'Deg',
                                     })

    dgroups = dfData.groupby('Name', observed=True)

    tmpDir = tempfile.TemporaryDirectory()
    print(tmpDir.name)

    plt.ioff()

    Args = []
    for ic, (gn, gg) in enumerate(dgroups):
        for index, row in gg.iterrows():
            Args.append({'tmpDir': tmpDir,
                         'Cl': row.CharCl,
                         'title': '{} \n {}'.format(gn, row.Date),
                         'dfAttrs': dfData.attrs})

    print('Generating Report')
    process_map(GenPSDBodeFigure, Args)

    plt.ion()

    merger = PdfMerger()
    FileName = tempfile.NamedTemporaryFile(suffix='.pdf', delete=False)

    for pdf in glob.glob(tmpDir.name + '/*.pdf'):
        merger.append(pdf)

    print(FileName.name)

    merger.write(FileName.name)
    merger.close()
    tmpDir.cleanup()

    if platform.system() == 'Darwin':  # macOS
        subprocess.call(('open', FileName.name))
    elif platform.system() == 'Windows':  # Windows
        os.startfile(FileName.name)
    else:  # linux variants
        subprocess.call(('xdg-open', FileName.name))



def PDFForce():
    PDF = PdfPages('test.pdf')
    PDF.close()