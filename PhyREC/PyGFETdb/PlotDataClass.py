#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:31:02 2017

@author: aguimera
"""

import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as mpcolors
import matplotlib.cm as cmx
#import matplotlib.markers.MarkerStyle.filled_markers as pltMarks
from itertools import cycle
import numpy as np
import math
# import PyGFETdb.NoiseModel as noise
from scipy.interpolate import interp1d
import sys


class MyCycle():
    def __init__(self, iterable):
        self.iterable = iterable
        self.cy = cycle(self.iterable)

    def next(self):
        return next(self.cy)

    def reset(self):
        del self.cy
        self.cy = cycle(self.iterable)


class PyFETPlotBase:
    marks = MyCycle(('', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', 's', 'p',
                     '*', 'D'))
    lines = MyCycle(('-', '--', '-.', ':'))
    ColorSet = cmx.jet
    mark = marks.next()
    line = lines.next()

    def setNColors(self, ncolors):
        cmap = cmx.ScalarMappable(mpcolors.Normalize(vmin=0, vmax=ncolors),
                                  self.ColorSet)
        col = []
        for i in range(ncolors):
            col.append(cmap.to_rgba(i))
        self.colors = MyCycle(col)

    def NextColor(self):
        self.color = self.colors.next()

    def NextMark(self):
        self.mark = self.marks.next()

    def NextLine(self):
        self.line = self.lines.next()
    
    def CreateFigure(self, Size=(9, 6)):
        self.Fig = plt.figure(figsize=Size)

    def FigExists(self):
        return plt.fignum_exists(self.Fig.number)

    def AddAxes(self, AxNames, Xvar=None):
        self.Axs = {}

        if len(AxNames) == 1:
            r = 1
        else:
            r = 2

        gs = gridspec.GridSpec(r, int(math.ceil(float(len(AxNames))/2)))

        for ia, ax in enumerate(AxNames):
#            if self.AxsProp[ax][3]:continue
            self.Axs[ax] = self.Fig.add_subplot(gs[ia])

        self.SetAxesLabels()
        self.SetAxesXLabels(Xvar)
        self.Fig.tight_layout()
        plt.show()

    def ClearAxes(self):
        for ax in self.Axs.values():
            while ax.lines:
                ax.lines[0].remove()

    def SetAxesXLabels(self, Xvar=None):
#        print 'empty'
        a = None

    def SetAxesLabels(self, fontsize='medium', labelsize=5,
                      scilimits=(-2, 2), RdsLim=(100, 15e3)):

        Axs = self.Axs

        if 'Ids' in Axs:
            Axs['Ids'].set_ylabel('Ids[A]', fontsize=fontsize)
            Axs['Ids'].set_xlabel('Vgs[V]', fontsize=fontsize)
            Axs['Ids'].grid()
            Axs['Ids'].tick_params(axis='both', which='both',
                                   labelsize=labelsize)
            Axs['Ids'].ticklabel_format(axis='y', style='sci',
                                        scilimits=scilimits)

        if 'IdsPoly' in Axs:
            Axs['IdsPoly'].set_ylabel('Ids[A]', fontsize=fontsize)
            Axs['IdsPoly'].set_xlabel('Vgs[V]', fontsize=fontsize)
            Axs['IdsPoly'].grid()
            Axs['IdsPoly'].tick_params(axis='both', which='both',
                                       labelsize=labelsize)
            Axs['IdsPoly'].ticklabel_format(axis='y', style='sci',
                                            scilimits=scilimits)

        if 'Ig' in Axs:
            Axs['Ig'].set_ylabel('Ig[A]', fontsize=fontsize)
            Axs['Ig'].tick_params(axis='both', which='both',
                                  labelsize=labelsize)
            Axs['Ig'].ticklabel_format(axis='y', style='sci',
                                       scilimits=scilimits)

        if 'Rds' in Axs:
            Axs['Rds'].set_ylabel(r'$RDS [\Omega]$', fontsize=fontsize)
            Axs['Rds'].set_xlabel('Vgs[V]', fontsize=fontsize)
            Axs['Rds'].grid()
            Axs['Rds'].tick_params(axis='both', which='both',
                                   labelsize=labelsize)
            Axs['Rds'].ticklabel_format(axis='y', style='sci',
                                        scilimits=scilimits)
            if RdsLim:
                Axs['Rds'].set_ylim(RdsLim[0], RdsLim[1])

        if 'GMPoly' in Axs:
            Axs['GMPoly'].set_ylabel('Gm[S]', fontsize=fontsize)
            Axs['GMPoly'].set_xlabel('Vgs[V]', fontsize=fontsize)
            Axs['GMPoly'].grid()
            Axs['GMPoly'].tick_params(axis='both', which='both',
                                      labelsize=labelsize)
            Axs['GMPoly'].ticklabel_format(axis='y', style='sci',
                                           scilimits=scilimits)

        if 'Gm' in Axs:
            Axs['Gm'].set_ylabel('Gm[S]', fontsize=fontsize)
            Axs['Gm'].set_xlabel('Vgs[V]', fontsize=fontsize)
            Axs['Gm'].grid()
            Axs['Gm'].tick_params(axis='both', which='both',
                                  labelsize=labelsize)
            Axs['Gm'].ticklabel_format(axis='y', style='sci',
                                       scilimits=scilimits)

        if 'FEMn' in Axs:
            Axs['FEMn'].set_ylabel('n[1/cm2]', fontsize=fontsize)
            Axs['FEMn'].set_xlabel('Vgs[V]', fontsize=fontsize)
            Axs['FEMn'].grid()
            Axs['FEMn'].tick_params(axis='both', which='both',
                                    labelsize=labelsize)
            Axs['FEMn'].ticklabel_format(axis='y', style='sci',
                                         scilimits=scilimits)

        if 'FEMmu' in Axs:
            Axs['FEMmu'].set_ylabel('mu[cm2/Vs]', fontsize=fontsize)
            Axs['FEMmu'].set_xlabel('Vgs[V]', fontsize=fontsize)
            Axs['FEMmu'].grid()
            Axs['FEMmu'].tick_params(axis='both', which='both',
                                     labelsize=labelsize)
            Axs['FEMmu'].ticklabel_format(axis='y', style='sci',
                                          scilimits=scilimits)

        if 'FEMmuGm' in Axs:
            Axs['FEMmuGm'].set_ylabel('mu[cm2/Vs]', fontsize=fontsize)
            Axs['FEMmuGm'].set_xlabel('Vgs[V]', fontsize=fontsize)
            Axs['FEMmuGm'].grid()
            Axs['FEMmuGm'].tick_params(axis='both', which='both',
                                       labelsize=labelsize)
            Axs['FEMmuGm'].ticklabel_format(axis='y', style='sci',
                                            scilimits=scilimits)

        if 'GM' in Axs:
            Axs['GM'].set_ylabel('Gm[S]', fontsize=fontsize)
            Axs['GM'].set_xlabel('Vgs[V]', fontsize=fontsize)
            Axs['GM'].grid()
            Axs['GM'].tick_params(axis='both', which='both',
                                  labelsize=labelsize)
            Axs['GM'].ticklabel_format(axis='y', style='sci',
                                       scilimits=scilimits)

        if 'GmMag' in Axs:
            Axs['GmMag'].set_ylabel('Gm[S]', fontsize=fontsize)
            Axs['GmMag'].set_xlabel('Frequency [Hz]', fontsize=fontsize)
            Axs['GmMag'].grid()
            Axs['GmMag'].tick_params(axis='both', which='both',
                                     labelsize=labelsize)
            Axs['GmMag'].ticklabel_format(axis='y', style='sci',
                                          scilimits=scilimits)

        if 'GmPh' in Axs:
            Axs['GmPh'].set_ylabel('Phase[0]', fontsize=fontsize)
            Axs['GmPh'].set_xlabel('Frequency [Hz]', fontsize=fontsize)
            Axs['GmPh'].grid()
            Axs['GmPh'].tick_params(axis='both', which='both',
                                    labelsize=labelsize)
            Axs['GmPh'].ticklabel_format(axis='y', style='sci',
                                         scilimits=scilimits)

        if 'PSD' in Axs:
            Axs['PSD'].set_ylabel('PSD [A2/Hz]', fontsize=fontsize)
            Axs['PSD'].set_xlabel('Frequency [Hz]', fontsize=fontsize)
            Axs['PSD'].grid()
            Axs['PSD'].tick_params(axis='both', which='both',
                                   labelsize=labelsize)
            Axs['PSD'].ticklabel_format(axis='y', style='sci',
                                        scilimits=scilimits)

        if 'NoA' in Axs:
            Axs['NoA'].set_ylabel('a [A^2]', fontsize=fontsize)
            Axs['NoA'].set_xlabel('Vgs [V]', fontsize=fontsize)
            Axs['NoA'].grid()
            Axs['NoA'].tick_params(axis='both', which='both',
                                   labelsize=labelsize)
            Axs['NoA'].ticklabel_format(axis='y', style='sci',
                                        scilimits=scilimits)

        if 'NoB' in Axs:
            Axs['NoB'].set_ylabel('b []', fontsize=fontsize)
            Axs['NoB'].set_xlabel('Vgs [V]', fontsize=fontsize)
            Axs['NoB'].grid()
            Axs['NoB'].tick_params(axis='both', which='both',
                                   labelsize=labelsize)
            Axs['NoB'].ticklabel_format(axis='y', style='sci',
                                        scilimits=scilimits)

        if 'Irms' in Axs:
            Axs['Irms'].set_ylabel('Irms [Arms]', fontsize=fontsize)
            Axs['Irms'].set_xlabel('Vgs [V]', fontsize=fontsize)
            Axs['Irms'].grid()
            Axs['Irms'].tick_params(axis='both', which='both',
                                    labelsize=labelsize)
            Axs['Irms'].ticklabel_format(axis='y', style='sci',
                                         scilimits=scilimits)

        if 'IrmsIds' in Axs:
            Axs['IrmsIds'].set_ylabel('Irms [Arms]', fontsize=fontsize)
            Axs['IrmsIds'].set_xlabel('Ids [A]', fontsize=fontsize)
            Axs['IrmsIds'].grid()
            Axs['IrmsIds'].tick_params(axis='both', which='both',
                                       labelsize=labelsize)
            Axs['IrmsIds'].ticklabel_format(axis='y', style='sci',
                                            scilimits=scilimits)

        if 'Vrms' in Axs:
            Axs['Vrms'].set_ylabel('Vrms [Vrms]', fontsize=fontsize)
            Axs['Vrms'].set_xlabel('Vgs [V]', fontsize=fontsize)
            Axs['Vrms'].grid()
            Axs['Vrms'].tick_params(axis='both', which='both',
                                    labelsize=labelsize)
            Axs['Vrms'].ticklabel_format(axis='y', style='sci',
                                         scilimits=scilimits)

        if 'FitErrA' in Axs:
            Axs['FitErrA'].set_ylabel('a fit error', fontsize=fontsize)
            Axs['FitErrA'].set_xlabel('Vgs [V]', fontsize=fontsize)
            Axs['FitErrA'].grid()
            Axs['FitErrA'].tick_params(axis='both', which='both',
                                       labelsize=labelsize)
            Axs['FitErrA'].ticklabel_format(axis='y', style='sci',
                                            scilimits=scilimits)

        if 'FitErrB' in Axs:
            Axs['FitErrB'].set_ylabel('b fit error', fontsize=fontsize)
            Axs['FitErrB'].set_xlabel('Vgs [V]', fontsize=fontsize)
            Axs['FitErrB'].grid()
            Axs['FitErrB'].tick_params(axis='both', which='both',
                                       labelsize=labelsize)
            Axs['FitErrB'].ticklabel_format(axis='y', style='sci',
                                            scilimits=scilimits)

        if 'GMax' in Axs:
            Axs['GMax'].set_ylabel('Max. transconductance [S]',
                                   fontsize=fontsize)
            Axs['GMax'].set_xlabel('Cicle', fontsize=fontsize)
            Axs['GMax'].grid()
            Axs['GMax'].tick_params(axis='both', which='both',
                                    labelsize=labelsize)
            Axs['GMax'].ticklabel_format(axis='y', style='sci',
                                         scilimits=scilimits)

        if 'Ud0' in Axs:
            Axs['Ud0'].set_ylabel('Vd [V]', fontsize=fontsize)
            Axs['Ud0'].grid()
            Axs['Ud0'].tick_params(axis='both', which='both',
                                   labelsize=labelsize)
            Axs['Ud0'].ticklabel_format(axis='y', style='sci',
                                        scilimits=scilimits)

        if 'Ud' in Axs:
            Axs['Ud'].set_ylabel('Dirac Point [V]',  fontsize=fontsize)
            Axs['Ud'].set_xlabel('Cicle', fontsize=fontsize)
            Axs['Ud'].grid()
            Axs['Ud'].tick_params(axis='both', which='both',
                                  labelsize=labelsize)
            Axs['Ud'].ticklabel_format(axis='y', style='sci',
                                       scilimits=scilimits)

        if 'Imin' in Axs:
            Axs['Imin'].set_ylabel('Min. Current [A]', fontsize=fontsize)
            Axs['Imin'].set_xlabel('Cicle', fontsize=fontsize)
            Axs['Imin'].grid()
            Axs['Imin'].tick_params(axis='both', which='both',
                                    labelsize=labelsize)
            Axs['Imin'].ticklabel_format(axis='y', style='sci',
                                         scilimits=scilimits)

    def AddLegend(self, Axn=None, fontsize='xx-small'):

        if Axn:
            Ax = self.Axs[Axn]
        else:
            Ax = self.Axs[list(self.Axs.keys())[0]]

        Ax.legend(fontsize=fontsize,
                  ncol=3, framealpha=0.2, loc=0)


class PyFETPlotParam(PyFETPlotBase):

    xVarProp = {'Length': ('TrtTypes', 'Length', 1),
                'Width': ('TrtTypes', 'Width', 1),
                'Pass': ('TrtTypes', 'Pass', 0),
                'W/L': (None, None, 1),
                'Area': ('TrtTypes', 'Area', 1),
                'Date': (None, 'DateTime', 0)}

                ### (Ispoly, NormVar,logscale)
    AxsProp = {'Vrms': (0, None, 1),
               'Irms': (0, None, 1),
               'NoA': (0, None, 1),
               'FitErrA': (0, None, 1),
               'FitErrB': (0, None, 1),
               'NoB': (0, None, 0),
               'GMPoly': (1, 'Vds', 0),
               'IdsPoly': (1, None, 0),
               'GmMax': (0, None, 0),
               'Ud0': (0, None, 0),
               'Rds': (0, None, 0),
               }

    def __init__(self):
        self.CreateFigure()
                ### (Name, X variable, Ispoly, OtherAxs)

    def PlotDataSet(self, Data, Trts, xVar, Bias, PltUd0=False):

        self.setNColors(len(Trts))
        for TrtN in sorted(Trts):
            self.NextColor()
            for cyn, cy in Data[TrtN].items():
                try:
                    self.Plot(cy, xVar, Bias, PltUd0)
                except:  # catch *all* exceptions
                    print (TrtN, cyn, sys.exc_info()[0])

    def Plot(self, Data, xVar, Bias, PltUd0):

#        Vds = Bias[0]
        BVgs = Bias[1]

        # Set Vds iterator
        if Bias[0]:
            vdind = Bias[0] ### TODO find correct Vds
        else:
            vdind = range(len(Data['Vds']))
        SLvds = ['Vd{}'.format(v) for v in vdind]

        self.marks.reset()
        self.NextMark() ### TODO WARNING will fail at 0 ''
        for svds, ivd in zip(SLvds, vdind):
            self.NextMark()

            for axn, ax in self.Axs.items():
                if PltUd0:
                    Vgs = BVgs + Data['Ud0'][ivd]
                else:
                    Vgs = BVgs

                # Calc X value
                if self.xVarProp[xVar][1]:
                    if self.xVarProp[xVar][0]:
                        valx = Data[self.xVarProp[xVar][0]][self.xVarProp[xVar][1]]
                    else:
                        valx = Data[self.xVarProp[xVar][1]]
                elif xVar == 'W/L':
                    valx = Data['TrtTypes']['Width']/Data['TrtTypes']['Length']

                # Calc Y value
                if axn == 'GmMax':
                    v = np.polyval(Data['GMPoly'][:, ivd], Data['Vgs'])
                    valy = np.max(np.abs(v))/Data['Vds'][ivd]
                elif axn == 'Rds':
                    v = np.polyval(Data['IdsPoly'][:, ivd], Vgs)
                    valy = Data['Vds'][ivd]/v
                elif axn == 'Vrms':
                    if 'GMPoly' in Data:
                        gm = np.polyval(Data['GMPoly'][:, ivd], Data['Vgs'])
                        valy = interp1d(Data['Vgs'],
                                        Data['Irms'][:, ivd]/np.abs(gm))(Vgs)
                    else:
                        continue
                elif axn == 'Ud0':
                    valy = Data['Ud0'][ivd]
                elif self.AxsProp[axn][0]:          # Polynom
                    valy = np.polyval(Data[axn][:, ivd], Vgs)
                else:
                    valy = interp1d(Data['Vgs'], Data[axn][:, ivd])(Vgs)

                # Plot data
                ax.plot(valx, valy, self.mark, color=self.color)

                if self.AxsProp[axn][2]:
                    ax.set_yscale('log')
                if self.xVarProp[xVar][2]:
                    ax.set_xscale('log')

    def SetAxesXLabels(self, Xvar, fontsize='medium', scilimits=(-2, 2)):
        for axn, ax in self.Axs.items():
            ax.set_xlabel(Xvar, fontsize=fontsize)
            ax.ticklabel_format(axis='x', style='sci', scilimits=scilimits)


class PyFETPlot(PyFETPlotBase):

                # (logY, X variable, Ispoly, OtherAxs)
    AxsProp = {'Vrms': (1, 'Vgs', 0),
               'IrmsIds': (2, 'Ids', 0),
               'Irms': (1, 'Vgs', 0),
               'NoA': (1, 'Vgs', 0),
               'FitErrA': (1, 'Vgs', 0),
               'FitErrB': (1, 'Vgs', 0),
               'NoB': (0, 'Vgs', 0),
               'GMPoly': (0, 'Vgs', 1),
               'IdsPoly': (0, 'Vgs', 1),
               'Ids': (0, 'Vgs', 1),
               'Ig': (0, 'Vgs', 0),
               'Gm': (0, 'Vgs', 1),
               'Rds': (0, 'Vgs', 1),
               'PSD': (2, 'Fpsd', 0),
               'GmMag': (2, 'Fgm', 0),
               'GmPh': (3, 'Fgm', 0)}

    ColorParams = {'Contact': ('TrtTypes', 'Contact'),
                   'Length': ('TrtTypes', 'Length'),
                   'Width': ('TrtTypes', 'Width'),
                   'Pass': ('TrtTypes', 'Pass'),
                   'W/L': (None, None),
                   'Trt': (None, 'Name'),
                   'Date': (None, 'DateTime'),
                   'Ud0': (None, 'Ud0')} #TODO fix with arrays

    def __init__(self):
        self.CreateFigure()

    def GetColorValue(self, cy, ColorOn):

        if self.ColorParams[ColorOn][1]:
            if self.ColorParams[ColorOn][0]:
                v = cy[self.ColorParams[ColorOn][0]][self.ColorParams[ColorOn][1]]
            else:
                v = cy[self.ColorParams[ColorOn][1]]
        elif ColorOn == 'W/L':
            v = cy['TrtTypes']['Width']/cy['TrtTypes']['Length']

        return v

    def PlotDataCh(self, Data, PltUd0=False, PltIsOK=False):

        self.setNColors(len(Data))
        for Trtv in Data.values():
            self.color = self.NextColor()
#            self.Plot(Trtv)
            try:
                self.Plot(Trtv, PltUd0=PltUd0, PltIsOK=PltIsOK)
            except:  # catch *all* exceptions
                print (sys.exc_info()[0])

    def PlotDataSet(self, Data, Trts, PltUd0=False, PltIsOK=False,
                    ColorOn='Trt'):

        Par = []
        for TrtN in sorted(Trts):
            for cyn, cy in Data[TrtN].items():
                Par.append(self.GetColorValue(cy, ColorOn))

        Par = sorted(set(Par))
        self.setNColors(len(Par))
        ColPar = {}
        for p in Par:
            self.NextColor()
            ColPar[p] = self.color

        for TrtN in sorted(Trts):
            self.marks.reset()

            for cyn, cy in Data[TrtN].items():
                self.color = ColPar[self.GetColorValue(cy, ColorOn)]
                self.NextMark()
                try:
                    self.Plot(cy, PltUd0=PltUd0, PltIsOK=PltIsOK)
                except:  # catch *all* exceptions
                    print (TrtN, cyn, sys.exc_info()[0])

    def Plot(self, Data, iVds=None, iVgs=None,
             PltUd0=False, PltIsOK=False, ColorOnVgs=False):

        label = Data['Name']

        # Set Vds iterator
        if iVds:
            vdind = iVds
        else:
            vdind = range(len( Data['Vds']))
        SLvds = ['Vd{}'.format(v) for v in vdind]

        self.lines.reset()
        for svds, ivd in zip(SLvds, vdind):
            self.NextLine()

            for axn, ax in self.Axs.items():
                if label == 'Gate' and axn is not 'Ig': continue
                if label is not 'Gate' and axn is 'Ig': continue

                # Calc X value
                if PltUd0:
                    Valx = Data['Vgs'] - Data['Ud0'][ivd]
                    Valxp = Data['Vgs']
                elif self.AxsProp[axn][1] == 'Ids':
                    if 'IdsPoly' in Data:
                        Valx = np.polyval(Data['IdsPoly'][:, ivd], Data['Vgs'])
                else:
                    Valx = Data[self.AxsProp[axn][1]]
                    Valxp = Valx

                # Set Vgs iterator
                if not self.AxsProp[axn][1] == 'Vgs':
                    if self.AxsProp[axn][1] == 'Ids':
                        ivgs = (0, )
                    elif iVgs:
                        ivgs = iVgs
                    else:
                        ivgs = range(len(Data['Vgs']))
                else:
                    ivgs = (0, )

                if ColorOnVgs:
                    self.setNColors(len(ivgs))

                for ivg in ivgs:
                    if ColorOnVgs:
                        self.NextColor()

                    Mark = self.line + self.mark

                    if PltIsOK:
                        if not Data['IsOK']:
                            Mark = '+'

                    # Calc Y value
                    if axn == 'Vrms':
                        if 'GMPoly' in Data:
                            gm = np.polyval(Data['GMPoly'][:, ivd],
                                            Data['Vgs'])
                            Valy = Data['Irms'][:, ivd]/np.abs(gm)
                        else:
                            continue
                    elif axn == 'Ig':
                        if 'Ig' in Data:
                            Valy = Data['Ig'][:,ivd]
                    elif axn == 'GmMag':
                        Valy = np.abs(Data['gm'][svds][ivg, :])
                    elif axn == 'GmPh':
                        Valy = np.angle(Data['gm'][svds][ivg, :])*180/np.pi
                    elif axn == 'PSD':
                        Valy = Data[axn][svds][ivg, :]
                        if 'NoA' in Data:
                            ax.loglog(Valx[1:], noise.Fnoise(Valx[1:],
                                      Data['NoA'][ivg, ivd],
                                      Data['NoB'][ivg, ivd]), '--')
                    elif axn == 'Gm':
                        if 'GMPoly' in Data:
                            Valy = np.polyval(Data['GMPoly'][:, ivd],
                                              Data['Vgs'])
                        else:
                            Valy = np.diff(Data['Ids'][:,ivd])/np.diff(Data['Vgs'])
                            Valx = Valx[1:]  # To check
                    elif axn == 'Ids':
                        Valy = Data['Ids'][:, ivd]
                        if 'IdsPoly' in Data:
                            ax.plot(Data['Vgs'],
                                    np.polyval(Data['IdsPoly'][:, ivd],
                                               Data['Vgs']), 'k-', alpha=0.3)
                    elif axn == 'Rds':
                        Valy = Data['Vds'][ivd]/Data['Ids'][:, ivd]
                    elif axn == 'IrmsIds':
                        Valy = Data['Irms'][:, ivd]
                    elif self.AxsProp[axn][2]:  # Polynom
                        Valy = np.polyval(Data[axn][:, ivd], Valxp)
                    else:
                        Valy = Data[axn][:, ivd]

                    # Plot data
                    if self.AxsProp[axn][0] == 1:
                        ax.semilogy(Valx, Valy,
                                    Mark,
                                    color=self.color,
                                    label=label)
                    elif self.AxsProp[axn][0] == 2:
                        ax.loglog(Valx, Valy,
                                  Mark,
                                  color=self.color,
                                  label=label)
                    elif self.AxsProp[axn][0] == 3:
                        ax.semilogx(Valx, Valy,
                                    Mark,
                                    color=self.color,
                                    label=label)
                    else:
                        ax.plot(Valx, Valy,
                                Mark,
                                color=self.color,
                                label=label)
