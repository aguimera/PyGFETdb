import glob
import math
import os
import platform
import subprocess
import tempfile
import warnings
from multiprocessing import Pool

import matplotlib as mpl
import numpy as np
import pandas as pd
from PyPDF2 import PdfMerger
from PyQt5.QtCore import QAbstractTableModel, QModelIndex
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from qtpy.QtCore import Qt
from scipy import stats
from tqdm.contrib.concurrent import process_map
import seaborn as sns


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


def GenFigure(PlotPars):
    nRows = math.ceil(np.sqrt(len(PlotPars)))
    nCols = math.ceil(len(PlotPars) / nRows)
    fig, Axs = plt.subplots(nrows=nRows, ncols=nCols)
    if nRows == 1 and nCols == 1:
        Axs = [Axs, ]
    else:
        Axs = Axs.flatten()

    return fig, Axs


def FormatAxis(PlotPars, dfAttrs):
    # figure generation
    fig, Axs = GenFigure(PlotPars)

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


def GenScalarFigures(Data, PlotPars, Xvar, Huevar, PltFunct):
    """
        BoxPlotFunctions = {'Boxplot': sns.boxplot,
                        'violinplot': sns.violinplot,
                        'swarmplot': sns.swarmplot,
                        'boxenplot': sns.boxenplot,
                        'barplot': sns.barplot,
                        'stripplot': sns.stripplot,
                        }
    """
    fig, Axs = GenFigure(PlotPars)

    for ic, p in enumerate(PlotPars):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            PltFunct(x=Xvar,
                     y=p,
                     hue=Huevar,
                     data=Data,
                     ax=Axs[ic])
        if p in LogPars:
            Axs[ic].set_yscale('log')

        leg = Axs[ic].get_legend()
        if leg is not None:
            plt.setp(leg.get_texts(),
                     fontsize='xx-small')
            plt.setp(leg.get_title(),
                     fontsize='xx-small')
        plt.setp(Axs[ic].get_xticklabels(),
                 rotation=45,
                 ha="right",
                 fontsize='xx-small',
                 rotation_mode="anchor")
        plt.setp(Axs[ic].get_yticklabels(),
                 fontsize='xx-small')

    return fig


def GenVectorFigures(data, GroupBy, PlotPars, Lines=True, Mean=False, Std=False):
    dgroups = data.groupby(GroupBy, observed=True)

    # Colors iteration
    Norm = mpl.colors.Normalize(vmin=0, vmax=len(dgroups))
    cmap = mpl.cm.ScalarMappable(norm=Norm, cmap='jet')

    fig, AxsDict = FormatAxis(PlotPars, data.attrs)

    for ic, gn in enumerate(dgroups.groups):
        gg = dgroups.get_group(gn)
        Col = cmap.to_rgba(ic)
        for parn in PlotPars:
            if parn in data.attrs['ArrayColsNorm']:
                xVar = data.attrs['VgsNorm']
            else:
                xVar = data.attrs['Vgs']

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
                    if Lines:
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

            if Mean:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    mean = np.nanmedian(Vals, axis=1)
                ax.plot(xVar, mean, '-.', color=Col, lw=1.5, label=gn)

            if Std:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    std = stats.tstd(Vals, axis=1)  # changed to mad for robustness
                ax.fill_between(xVar, mean + std, mean - std, color=Col, alpha=0.2)

            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys(), fontsize='xx-small')

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


def AddColRow(dfData):
    Col = []
    Row = []
    nCol = []
    nRow = []
    for ir, d in dfData.iterrows():
        n = d['Name'].split('-')
        Col.append(n[-1])
        Row.append(n[-2])
        nCol.append(int(n[-1][-2:]))
        nRow.append(int(n[-2][-2:]))

    dfData['Col'] = Col
    dfData['Row'] = Row
    dfData['nCol'] = nCol
    dfData['nRow'] = nRow


def AddCycle(dfData):
    gTrts = dfData.groupby('Name')
    dfData.insert(2, 'Cycle', np.zeros(len(dfData), dtype=int))
    for gn, gt in gTrts:
        for cy, (ir, d) in enumerate(gt.iterrows()):
            dfData.at[ir, 'Cycle'] = cy


def GenDeviceFigures(Data, Name, ScalarPars, ScalarPltFunct, VectorPars, **kwargs):
    Figs = []
    Fig = GenScalarFigures(Data=Data,
                           PlotPars=ScalarPars,
                           Xvar='Cycle',
                           Huevar=None,
                           PltFunct=ScalarPltFunct)
    Fig.set_size_inches(11, 8)
    Fig.suptitle(Name)
    Fig.tight_layout()
    Figs.append(Fig)
    Fig, axdict = GenVectorFigures(data=Data,
                                   GroupBy='Cycle',
                                   PlotPars=VectorPars)

    Fig.set_size_inches(11, 8)
    Fig.suptitle(Name)
    Fig.tight_layout()
    plt.close(Fig)

    if 'Col' not in Data.keys():
        return Figs

    gCys = Data.groupby('Cycle')
    for gcn, g in gCys:
        Fig, Axs = GenFigure(ScalarPars)
        for ic, p in enumerate(ScalarPars):
            sns.heatmap(g.pivot('Col', 'Row', p),
                        ax=Axs[ic])
            Axs[ic].set_title(p)

        Fig.suptitle(str(Name) + ' Cycle ' + str(gcn))
        Fig.set_size_inches(11, 8)
        Fig.tight_layout()
        Figs.append(Fig)

    return Name, Figs


def GenDeviceReport(Data, FileOut, ScalarPars, ScalarPltFunct, VectorPars):
    gDevs = Data.groupby('Device')

    PDF = PdfPages(FileOut)
    plt.ioff()
    Args = []
    for ig, (gn, gd) in enumerate(gDevs):
        print('Generating {} - {} of {}'.format(gn, ig, len(gDevs)))
        # Figs = GenDeviceFigures(gd, gn, ScalarPars, ScalarPltFunct, VectorPars)
        Args.append((gd, gn, ScalarPars, ScalarPltFunct, VectorPars))

    with Pool() as p:
        Result = p.starmap_async(GenDeviceFigures, Args)
        for ic, (Name, Figs) in enumerate(Result.get()):
            print('Get Figures from {}, {} of {}'.format(Name, ic, len(Args)))
            for f in Figs:
                PDF.savefig(f)
                plt.close(f)

    plt.ion()
    PDF.close()


def GenDeviceReportBrigde(Args):
    Name, Figs = GenDeviceFigures(**Args)
    for ic, f in enumerate(Figs):
        file = tempfile.NamedTemporaryFile(dir=Args['tmpDir'].name,
                                           prefix=Name + '{0:02}'.format(ic),
                                           suffix='.pdf',
                                           delete=False)
        f.savefig(file, format='pdf')
        plt.close(f)


def GenDeviceReportGui(Data, FileOut, ScalarPars, ScalarPltFunct, VectorPars):
    gDevs = Data.groupby('Device')

    tmpDir = tempfile.TemporaryDirectory()
    print(tmpDir.name)

    plt.ioff()

    Args = []
    for ig, (gn, gd) in enumerate(gDevs):
        print('Generating {} - {} of {}'.format(gn, ig, len(gDevs)))
        # Figs = GenDeviceFigures(gd, gn, ScalarPars, ScalarPltFunct, VectorPars)
        Args.append({'Data': gd,
                     'Name': gn,
                     'ScalarPars': ScalarPars,
                     'ScalarPltFunct': ScalarPltFunct,
                     'VectorPars': VectorPars,
                     'tmpDir': tmpDir})

    print('Generating Report')
    process_map(GenDeviceReportBrigde, Args)

    merger = PdfMerger()
    files = glob.glob(tmpDir.name + '/*.pdf')
    for pdf in sorted(files):
        merger.append(pdf)
    merger.write(FileOut)
    merger.close()
    tmpDir.cleanup()

    plt.ion()
