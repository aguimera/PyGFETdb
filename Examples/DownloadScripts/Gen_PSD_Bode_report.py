import glob
import subprocess

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from PyGFETdb.GuiDBView.GuiHelpers import FormatAxis
from PyPDF2 import PdfMerger
import tempfile
from tqdm.contrib.concurrent import process_map


def AddBodePSDFigure(Args):
    Cl = Args['Cl']
    tmpDir = Args['tmpDir']
    Vgs = Cl.GetVgs()
    Norm = mpl.colors.Normalize(vmin=0, vmax=len(Vgs))
    cmap = mpl.cm.ScalarMappable(norm=Norm, cmap='jet')

    fig, AxsDict = FormatAxis(('Ids', 'BodeMag', 'PSD', 'BodePhase'), dfAttrs)
    fig.set_size_inches(11, 8)
    fig.suptitle('{} \n {}'.format(gn, row.Date))

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


if __name__ == '__main__':

    mpl.use("Qt5Agg")
    plt.close('all')

    FileName = 'Report.pdf'
    dfData = pd.read_pickle('Data.pkl')
    dfData.attrs['ColUnits'].update({'PSD': 'A**2/Hz',
                                     'BodeMag': 'S',
                                     'BodePhase': 'Deg',
                                     })
    dfAttrs = dfData.attrs

    dgroups = dfData.groupby('Name', observed=True)

    tmpDir = tempfile.TemporaryDirectory()
    print(tmpDir.name)

    plt.ioff()

    Args = []
    for ic, (gn, gg) in enumerate(dgroups):
        # if ic > 20:
        #     break
        for index, row in gg.iterrows():
            Args.append({'tmpDir': tmpDir,
                         'Cl': row.CharCl})

    print('Generating Report')
    process_map(AddBodePSDFigure, Args)

    plt.ion()

    merger = PdfMerger()

    for pdf in glob.glob(tmpDir.name + '/*.pdf'):
        merger.append(pdf)

    merger.write(FileName)
    merger.close()
    tmpDir.cleanup()
    subprocess.call(['xdg-open', FileName])
