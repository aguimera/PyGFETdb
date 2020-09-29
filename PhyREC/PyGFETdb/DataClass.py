#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 17:43:50 2017

@author: aguimera
"""

import numpy as np
from scipy.interpolate import interp1d
import scipy.optimize as optim
from scipy.integrate import simps
import quantities as pq

from PhyREC.PyGFETdb import PlotDataClass
import sys

DebugPrint = True

DefaultUnits = {'Vds': pq.V,
                'Ud0': pq.V,
#                'PSD': pq.A**2/pq.Hz,
#                'Fgm': pq.S,
#                'gm': pq.S,
                'Vgs': pq.V,
#                'Fpsd': pq.Hz,
#                'Ig': pq.A,
                'Irms': pq.V,
                'Ids': pq.A,
                }


class DataCharDC(object):
    PolyOrder = 10
#    FEMn0 = 8e11  # 1/cm^2
#    FEMq = 1.602176565e-19
#    FEMRc = 300  # absolute value Ohms
#    FEMCdl = 2e-6  # F/cm^2
    IntMethod = 'linear'

    def __init__(self, Data):
        for k, v in Data.items():
            if k == 'Gate':
                if v is None:
                    if DebugPrint:
                        print ('No gate values')
                elif np.isnan(v['Ig']).any():
                    if DebugPrint:
                        print ('NaN in gate values')
                else:
                    self.__setattr__('Ig', v['Ig'])
#            if k in DefaultUnits:
##                print(k, v, DefaultUnits[k])
#                v = v * DefaultUnits[k]
            self.__setattr__(k, v)

        if 'Ud0' not in self.__dict__:
            self.CalcUd0()
           
    def _FormatOutput(self, Par, **kwargs):
        if 'Units' in kwargs:
            Par.rescale(kwargs['Units'])

        if not hasattr(Par, '__iter__'):
            return Par[None, None]
        s = Par.shape
        if len(s) == 0:
            return Par[None, None]
        if len(s) == 1:
            return Par[:, None]

        return Par.transpose()
     
    def UpdateData(self, Data):
        for k, v in Data.items():
            self.__setattr__(k, v)

    def CalcUd0(self):
        if 'IdsPoly' not in self.__dict__:
            self.CalcIdsPoly()

        vgs = np.linspace(self.Vgs[0], self.Vgs[-1], len(self.Vgs)*10000)

        self.Ud0 = np.ones(len(self.Vds))*np.NaN
        for ivd, Vds in enumerate(self.Vds):
            self.Ud0[ivd] = vgs[np.argmin(np.polyval(self.IdsPoly[:, ivd],
                                                     vgs))]

    def CalcIdsPoly(self, PolyOrder=None):
        if PolyOrder:
            self.PolyOrder = PolyOrder
        Ord = self.PolyOrder

        self.IdsPoly = np.ones((Ord+1, len(self.Vds)))*np.NaN
        for ivd, Vds in enumerate(self.Vds):
            self.IdsPoly[:, ivd] = np.polyfit(self.Vgs, self.Ids[:, ivd], Ord)

    def CalcGMPoly(self, PolyOrder=None):
        if 'IdsPoly' not in self.__dict__:
            self.CalcIdsPoly(PolyOrder=PolyOrder)

        if PolyOrder:
            self.PolyOrder = PolyOrder
        Ord = self.PolyOrder

        self.GMPoly = np.ones((Ord, len(self.Vds)))*np.NaN

        for ivd, Vds in enumerate(self.Vds):
            self.GMPoly[:, ivd] = np.polyder(self.IdsPoly[:, ivd])

    def CalcFEM(self,
                FEMn0=8e11,  # 1/cm^2
                FEMq=1.602176565e-19,
                FEMRc=300,  # absolute value Ohms
                FEMCdl=2e-6,
                FEMRcVgs=None, **kwargs):  # F/cm^2)
        # TODO interpolate IDSpoly from all vgs....
        if 'IdsPoly' not in self.__dict__:
            self.CalcIdsPoly()

        self.FEMn = np.ones((len(self.Vgs), len(self.Vds)))*np.NaN
        self.FEMmu = np.ones((len(self.Vgs), len(self.Vds)))*np.NaN
        self.FEMmuGm = np.ones((len(self.Vgs), len(self.Vds)))*np.NaN

        if FEMRcVgs is not None:
            FEMRc = np.ones(len(self.Vgs))*np.NaN 
            RcVgs = FEMRcVgs[0, :]
            RcVgsRc = FEMRcVgs[1, :]
            rcint = interp1d(RcVgs, RcVgsRc)
            Vgmeas = self.GetVgs(Ud0Norm=True)
            VgInds =  np.where((Vgmeas>np.min(RcVgs)) & (Vgmeas<np.max(RcVgs)))[0]
            FEMRc[VgInds] = rcint(Vgmeas[VgInds,0])
            print (FEMRc)


        L = self.TrtTypes['Length']
        W = self.TrtTypes['Width']

        VgUd = np.abs(self.GetVgs(Ud0Norm=True))
        Ids = self.GetIds()
        Gm = np.abs(self.GetGM())
        for ivd, Vds in enumerate(self.Vds):
            n = (FEMCdl * VgUd[:, ivd])/FEMq
            self.FEMn[:, ivd] = np.sqrt(n**2 + FEMn0**2)

            Ieff = Vds/(Vds/Ids[:, ivd] - FEMRc)
            mu = (Ieff*L)/(W*Vds*n*FEMq)
            self.FEMmu[:, ivd] = mu

            Vdseff = Vds - Ids[:, ivd]*FEMRc
            muGM = (Gm[:, ivd]*L)/(FEMCdl*Vdseff*W)
            self.FEMmuGm[:, ivd] = muGM

    def GetUd0(self, Vds=None, Vgs=None, Ud0Norm=False,
               Normalize=False, **kwargs):
        if 'Ud0' not in self.__dict__:
            self.CalcUd0()

        iVds = self.GetVdsIndexes(Vds)
        if len(iVds) == 0:
            return None

        Ud0 = np.array([])
        for ivd in iVds:
            ud0 = self.Ud0[ivd]
            if Normalize:
                ud0 = ud0-(self.Vds[ivd]/2)
            Ud0 = np.vstack((Ud0, ud0)) if Ud0.size else ud0

        return self._FormatOutput(Ud0, **kwargs)

    def GetDateTime(self, **kwargs):
        return self.DateTime

    def GetTime(self, **kwargs):
        return np.datetime64(self.DateTime)[None, None].transpose()

    def GetVds(self, **kwargs):
        return self._FormatOutput(self.Vds, **kwargs)

    def GetVgs(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        if not Ud0Norm:
            return self.Vgs

        iVds = self.GetVdsIndexes(Vds)
        if len(iVds) == 0:
            return None

        if 'Ud0' not in self.__dict__:
            self.CalcUd0()

        Vgs = np.array([])
        for ivd in iVds:
            vgs = self.Vgs - self.Ud0[ivd]
            Vgs = np.vstack((Vgs, vgs)) if Vgs.size else vgs

        return self._FormatOutput(Vgs, **kwargs)

    def GetVdsIndexes(self, Vds):
        if Vds:
            if not hasattr(Vds, '__iter__'):
                Vds = (Vds,)
            iVds = []
            for vd in Vds:
                ind = np.where(self.Vds == vd)[0]
                if len(ind) > 0:
                    iVds.append(ind[0])
                else:
                    print ('Vds = ', vd, 'Not in data')
        else:
            iVds = range(len(self.Vds))
        return iVds

    def GetIds(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        iVds = self.GetVdsIndexes(Vds)
        if len(iVds) == 0:
            return None

        if Vgs is not None:
            vgs = Vgs  # TODO Check Vgs range
        else:
            vgs = self.Vgs

        if 'IdsPoly' not in self.__dict__:
            self.CalcIdsPoly()

        Ids = np.array([])
        for ivd in iVds:
            if Ud0Norm and Vgs is not None:
                vg = vgs + self.Ud0[ivd]
            else:
                vg = vgs
            ids = np.polyval(self.IdsPoly[:, ivd], vg)
            Ids = np.vstack((Ids, ids)) if Ids.size else ids

#        Ids = Ids * DefaultUnits['Ids']
        return self._FormatOutput(Ids, **kwargs)

    def GetGM(self, Vgs=None, Vds=None, Normalize=False, Ud0Norm=False, **kwargs):
        iVds = self.GetVdsIndexes(Vds)
        if len(iVds) == 0:
            return None

        if Vgs is not None:
            vgs = Vgs  # TODO Check Vgs range
        else:
            vgs = self.Vgs

        if 'GMPoly' not in self.__dict__:
            self.CalcGMPoly()

        GM = np.array([])
        for ivd in iVds:
            if Ud0Norm and Vgs is not None:
                vg = vgs + self.Ud0[ivd]
            else:
                vg = vgs
            gm = np.polyval(self.GMPoly[:, ivd], vg)
            if Normalize:
                gm = gm/self.Vds[ivd] #*(self.TrtTypes['Length']/self.TrtTypes['Width'])/
            GM = np.vstack((GM, gm)) if GM.size else gm

        return self._FormatOutput(GM, **kwargs)

    def GetGMV(self, AbsVal=True, **kwargs):
        kwargs.update({'Normalize': True})
        if AbsVal:
            return np.abs(self.GetGM(**kwargs))
        else:
            return self.GetGM(**kwargs)
        
    def GetGMMax(self, **kwargs):
        return np.max(np.abs(self.GetGM(**kwargs)))
    
    def GetRds(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        Ids = self.GetIds(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)

        if Ids is None:
            return None

        iVds = self.GetVdsIndexes(Vds)

        Rds = np.array([])
        for iid, ivd in enumerate(iVds):
            rds = self.Vds[ivd]/Ids[:, iid]
            Rds = np.vstack((Rds, rds)) if Rds.size else rds

        return self._FormatOutput(Rds, **kwargs)

    def GetFEMn(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        if 'FEMn' not in self.__dict__:
            self.CalcFEM(**kwargs)
        return self._GetParam('FEMn', Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)

    def GetFEMmu(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        if 'FEMmu' not in self.__dict__:
            self.CalcFEM(**kwargs)
        return self._GetParam('FEMmu', Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)

    def GetFEMmuGm(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        if 'FEMmuGm' not in self.__dict__:
            self.CalcFEM(**kwargs)
        return self._GetParam('FEMmuGm', Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)

    def GetIg(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        if 'Ig' not in self.__dict__:
#            print 'No Gate data'
            return None
        return self._GetParam('Ig', Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)

    def CheckVgsRange(self, Vgs, iVds, Ud0Norm):
        if Vgs is not None:
            for ivd in iVds:
                if Ud0Norm is None or Ud0Norm is False:
                    vg = Vgs
                    VgsM = self.Vgs
                else:
                    vg = Vgs + self.Ud0[ivd]
                    VgsM = self.Vgs

                if (np.min(vg) < np.min(VgsM)) or (np.max(vg) > np.max(VgsM)):
                    print (self.Name, 'Vgs range not valid', vg, VgsM, self.Ud0)
                    return None
            return Vgs
        else:
            return self.Vgs

    def _GetParam(self, Param, Vgs=None, Vds=None,
                  Ud0Norm=False, Normalize=False, **kwargs):

        iVds = self.GetVdsIndexes(Vds)
        if len(iVds) == 0:
            return None

        vgs = self.CheckVgsRange(Vgs, iVds, Ud0Norm)
        if vgs is None:
            return None
        if len(self.Vgs) < 2:
            print ('self Vgs len error', self.Vgs)
            return None

        if Param not in self.__dict__:
            print ('Not Data')
            return None
        Par = self.__getattribute__(Param)

        PAR = np.array([])
        for ivd in iVds:
            if Ud0Norm and Vgs is not None:
                vg = vgs + self.Ud0[ivd]
            else:
                vg = vgs

            par = interp1d(self.Vgs, Par[:, ivd], kind=self.IntMethod)(vg)
            if Normalize:
                par = par/self.Vds[ivd]
            PAR = np.vstack((PAR, par)) if PAR.size else par
        
        return self._FormatOutput(PAR, **kwargs)
    

    def GetName(self, **kwargs):
        return self.Name

    def GetWL(self, **kwargs):
        return self.TrtTypes['Width']/self.TrtTypes['Length']

    def GetPass(self, **kwargs):
        return self.TrtTypes['Pass']

    def GetLength(self, **kwargs):
        return self.TrtTypes['Length']

    def GetWidth(self, **kwargs):
        return self.TrtTypes['Width']

    def GetContact(self, **kwargs):
        return self.TrtTypes['Contact']

    def GetTypeName(self, **kwargs):
        return self.TrtTypes['Name']

    def GetPh(self, **kwargs):
        return np.array(self.Info['Ph'])[None, None]

    def GetIonStrength(self, **kwargs):
        return np.array(self.Info['IonStrength'])[None, None]

    def GetFuncStep(self, **kwargs):
        return self.Info['FuncStep']

    def GetComments(self, **kwargs):
        return self.Info['Comments']

    def GetAnalyteCon(self, **kwargs):
        return np.array(self.Info['AnalyteCon'])[None, None]

    def GetGds(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        Gds = 1 / self.GetRds(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
        return Gds

    def GetGMNorm(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        GMNorm = self.GetGM(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm, Normalize=True)
        return GMNorm

    def GetUd0Vds(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        return self.GetUd0(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm, Normalize=True)


def Fnoise(f, a, b):
    return a/f**b


def LogFnoise(f, a, b):
    return b*f+a


def FitFNoise(Freq, psd):
    poptV, pcov = optim.curve_fit(Fnoise, Freq, psd)
    a = poptV[0]
    b = poptV[1]

    return a, b, np.sqrt(np.diag(pcov))


def FitLogFnoise(Freq, psd):
    poptV, pcov = optim.curve_fit(LogFnoise, np.log10(Freq),
                                  np.log10(psd))

    a = 10 ** poptV[0]
    b = - poptV[1]
    return a, b, np.sqrt(np.diag(pcov))


class DataCharAC(DataCharDC):
    FFmin = None
    FFmax = None
    NFmin = None
    NFmax = None

    def _GetFreqVgsInd(self, Vgs=None, Vds=None, Ud0Norm=False):
        iVds = self.GetVdsIndexes(Vds)
        if len(iVds) == 0:
            return None, None

        vgs = self.CheckVgsRange(Vgs, iVds, Ud0Norm)
        if vgs is None:
            return None, None

# TODO check for more than 1 vds
        if Ud0Norm is True:
            vgs = vgs - self.Ud0[iVds[0]]

        VGS = self.GetVgs(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
        vgsmeas = [min(VGS, key=lambda x:abs(x-vg)) for vg in vgs]
        VgsInd = [np.where(VGS == vg)[0][0] for vg in vgsmeas]

        SiVds = ['Vd{}'.format(i) for i in iVds]
        return SiVds, VgsInd

    def _CheckFreqIndexes(self, Freq, Fmin, Fmax):
        if Fmin is not None and Fmax is not None:
            return np.where(np.logical_and(Freq > Fmin, Freq < Fmax))
        if Fmin is not None:
            return np.where(Freq > Fmin)
        if Fmax is not None:
            return np.where(Freq < Fmax)
        return range(len(Freq))[1:]

    def FitNoise(self, Fmin, Fmax):
        nVgs = len(self.Vgs)
        nVds = len(self.Vds)

        NoA = np.ones((nVgs, nVds))*np.NaN
        NoB = np.ones((nVgs, nVds))*np.NaN
        FitErrA = np.ones((nVgs, nVds))*np.NaN
        FitErrB = np.ones((nVgs, nVds))*np.NaN

        for ivd in range(nVds):
            for ivg in range(nVgs):
                psd = self.PSD['Vd{}'.format(ivd)][ivg, :]
                Fpsd = self.Fpsd
                if np.any(np.isnan(psd)):
                    continue

                try:
                    Inds = self._CheckFreqIndexes(Fpsd, Fmin, Fmax)
                    a, b, err = FitLogFnoise(Fpsd[Inds], psd[Inds])
                    NoA[ivg, ivd] = a
                    NoB[ivg, ivd] = b
                    FitErrA[ivg, ivd] = err[0]
                    FitErrB[ivg, ivd] = err[1]
                    self.NoA = NoA
                    self.NoB = NoB
                    self.FitErrA = FitErrA
                    self.FitErrB = FitErrB
                except:
                    print ("Fitting error:", sys.exc_info()[0])

    def CalcIRMS(self, Fmin, Fmax):
        nVgs = len(self.Vgs)
        nVds = len(self.Vds)

        Irms = np.ones((nVgs, nVds))*np.NaN
        for ivd in range(nVds):
            for ivg in range(nVgs):
                psd = self.PSD['Vd{}'.format(ivd)][ivg, :]
                Fpsd = self.Fpsd
                if np.any(np.isnan(psd)):
                    continue

                Inds = self._CheckFreqIndexes(Fpsd, Fmin, Fmax)
                Irms[ivg, ivd] = np.sqrt(simps(psd[Inds], Fpsd[Inds]))
        self.Irms = Irms

    def GetPSD(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        SiVds, VgsInd = self._GetFreqVgsInd(Vgs, Vds, Ud0Norm)
        if VgsInd is None:
            return None
        return self.PSD[SiVds[0]][VgsInd, :].transpose()

    def GetGmMag(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        SiVds, VgsInd = self._GetFreqVgsInd(Vgs, Vds, Ud0Norm)
        if VgsInd is None:
            return None

        return np.abs(self.gm[SiVds[0]][VgsInd, :].transpose())

    def GetGmPh(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        SiVds, VgsInd = self._GetFreqVgsInd(Vgs, Vds, Ud0Norm)
        if VgsInd is None:
            return None

        return np.angle(self.gm[SiVds[0]][VgsInd, :].transpose(), deg=True)

    def GetFpsd(self, **kwargs):
        return self.Fpsd

    def GetFgm(self, **kwargs):
        return self.Fgm

    def _CheckRMS(self, NFmin, NFmax):
        if NFmin is not None or NFmax is not None:
            if self.NFmin != NFmin or self.NFmax != NFmax:
                print ('Calc IRMS')
                self.NFmin = NFmin
                self.NFmax = NFmax
                if self.IsOK:
                    self.CalcIRMS(Fmin=NFmin, Fmax=NFmax)

    def GetIrms(self, Vgs=None, Vds=None, Ud0Norm=False,
                NFmin=None, NFmax=None, **kwargs):
        self._CheckRMS(NFmax=NFmax, NFmin=NFmin)
        return self._GetParam('Irms', Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)

    def GetVrms(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        Irms = self.GetIrms(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm, **kwargs)
        if Irms is None:
            return None
        gm = np.abs(self.GetGM(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm))
        return Irms/gm

    def _CheckFitting(self, FFmin, FFmax):
        if FFmin is not None or FFmax is not None:
            if self.FFmin != FFmin or self.FFmax != FFmax:
                print ('Calc fitting')
                self.FFmin = FFmin
                self.FFmax = FFmax
                if self.IsOK:
                    self.FitNoise(Fmin=FFmin, Fmax=FFmax)

    def GetNoA(self, Vgs=None, Vds=None, Ud0Norm=False,
               FFmin=None, FFmax=None, **kwargs):
        self._CheckFitting(FFmin, FFmax)
        return self._GetParam('NoA', Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)

    def GetNoB(self, Vgs=None, Vds=None, Ud0Norm=False,
               FFmin=None, FFmax=None, **kwargs):
        self._CheckFitting(FFmin, FFmax)
        return self._GetParam('NoB', Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)

    def GetNoAIds2(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        NoA = self._GetParam('NoA', Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
        Ids = self.GetIds(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
        return NoA/(Ids**2)

    def GetIrmsVds(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        return self._GetParam('Irms', Vgs=Vgs, Vds=Vds,
                              Ud0Norm=Ud0Norm, Normalize=True)

    def GetIrmsIds2(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        Irms = self._GetParam('Irms', Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
        Ids = self.GetIds(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
        return Irms/(Ids**2)

    def GetIrmsIds15(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        Irms = self._GetParam('Irms', Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
        Ids = self.GetIds(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
        return Irms/(Ids**1.5)

    def GetIrmsIds(self, Vgs=None, Vds=None, Ud0Norm=False, **kwargs):
        Irms = self._GetParam('Irms', Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
        Ids = self.GetIds(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm)
        return Irms/Ids



class PyFETPlotDataClass(PlotDataClass.PyFETPlotBase):

    # (logY, logX, X variable)
    AxsProp = {'Vrms': (1, 0, 'Vgs'),
               'Irms': (1, 0, 'Vgs'),
               'NoA': (1, 0, 'Vgs'),
               'FitErrA': (1, 0, 'Vgs'),
               'FitErrB': (1, 0, 'Vgs'),
               'NoB': (0, 0, 'Vgs'),
               'GM': (0, 0, 'Vgs'),
               'Ids': (0, 0, 'Vgs'),
               'Ig': (0, 0, 'Vgs'),
               'Rds': (0, 0, 'Vgs'),
               'FEMn': (0, 0, 'Vgs'),
               'FEMmu': (1, 0, 'Vgs'),
               'FEMmuGm': (1, 0, 'Vgs'),
               'PSD': (1, 1, 'Fpsd'),
               'GmMag': (1, 1, 'Fgm'),
               'GmPh': (0, 1, 'Fgm')}

    ColorParams = {'Contact': ('TrtTypes', 'Contact'),
                   'Length': ('TrtTypes', 'Length'),
                   'Width': ('TrtTypes', 'Width'),
                   'Pass': ('TrtTypes', 'Pass'),
                   'W/L': (None, None),
                   'Trt': (None, 'Name'),
                   'Date': (None, 'DateTime'),
                   'Ud0': (None, 'Ud0'),
                   'Device': ('TrtTypes', 'Devices.Name'),
                   'Wafer': ('TrtTypes', 'Wafers.Name')}  # TODO fix with arrays

    def __init__(self, Size=(9, 6)):
        self.CreateFigure(Size=Size)

    def GetColorValue(self, Data, ColorOn):
        if self.ColorParams[ColorOn][1]:
            if self.ColorParams[ColorOn][0]:
                p = Data.__getattribute__(self.ColorParams[ColorOn][0])
                v = p[self.ColorParams[ColorOn][1]]
            else:
                v = Data.__getattribute__(self.ColorParams[ColorOn][1])
        elif ColorOn == 'W/L':
            p = Data.__getattribute__('TrtTypes')
            v = p['Width']/p['TrtTypes']['Length']
        return v

    def PlotDataCh(self, DataDict, Trts, Vgs=None, Vds=None, Ud0Norm=False,
                   PltIsOK=True, ColorOn='Trt'):

        self.setNColors(len(DataDict))
        for Trtv in DataDict.values():
            self.color = self.NextColor()
            try:
                self.Plot(DataDict, Vgs=Vgs, Vds=Vds,
                          Ud0Norm=Ud0Norm, PltIsOK=PltIsOK)
            except:  # catch *all* exceptions
                print (sys.exc_info()[0])

    def PlotDataSet(self, DataDict, Trts=None,
                    Vgs=None, Vds=None, Ud0Norm=False,
                    PltIsOK=True, ColorOn='Trt', MarkOn='Cycle', **kwargs):

        if Trts is None:
            Trts = DataDict.keys()

        Par = []
        for TrtN in sorted(Trts):
            for Dat in DataDict[TrtN]:
                Par.append(self.GetColorValue(Dat, ColorOn))

        Par = sorted(set(Par))
        self.setNColors(len(Par))
        ColPar = {}
        for p in Par:
            self.NextColor()
            ColPar[p] = self.color

        for TrtN in sorted(Trts):
            self.marks.reset()
            for Dat in DataDict[TrtN]:
                self.color = ColPar[self.GetColorValue(Dat, ColorOn)]
                if MarkOn == 'Cycle':
                    self.NextMark()                    

                try:
                    self.Plot(Dat, Vgs=Vgs, Vds=Vds,
                              Ud0Norm=Ud0Norm, PltIsOK=PltIsOK, **kwargs)
                except:  # catch *all* exceptions
                    print (TrtN, sys.exc_info()[0])

    def Plot(self, Data, Vgs=None, Vds=None,
             Ud0Norm=False, PltIsOK=True, ColorOnVgs=False, **kwargs):

        label = Data.Name

        for axn, ax in self.Axs.items():
            Mark = self.line + self.mark
            if not Data.IsOK and PltIsOK:
                Mark = '+'

            if self.AxsProp[axn][2] == 'Vgs':
                if Vgs is None:
                    Valx = Data.GetVgs(Ud0Norm=Ud0Norm)
                else:
                    Valx = Vgs
            else:
                func = Data.__getattribute__('Get' + self.AxsProp[axn][2])
                Valx = func(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm, **kwargs)

            func = Data.__getattribute__('Get' + axn)
            Valy = func(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm, **kwargs)

            if Valx is not None and Valy is not None:
                ax.plot(Valx, Valy, Mark, color=self.color, label=label)

                if axn == 'PSD':
                    a = Data.GetNoA(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm,
                                    **kwargs)
                    b = Data.GetNoB(Vgs=Vgs, Vds=Vds, Ud0Norm=Ud0Norm,
                                    **kwargs)
                    Valy = Fnoise(Valx, a, b).transpose()
                    ax.plot(Valx, Valy, Mark, '--',
                            color=self.color, alpha=0.5)

                if self.AxsProp[axn][0]:
                    ax.set_yscale('log')
                if self.AxsProp[axn][1]:
                    ax.set_xscale('log')
