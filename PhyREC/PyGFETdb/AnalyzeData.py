# -*- coding: utf-8 -*-
"""
@author: Anton Guimerà
@version: 0.1b

Revsion history

- 0.1b -- First version

TODO implement graph for selected channels

"""
import numpy as np
# import PyGFETdb.NoiseModel as noise
import sys
from scipy import interpolate

###############################################################################
#####
###############################################################################
def CalcDCparams (DevDC):   #calculates Dirac,MaxGM, Imin, Rmin,.. 
    
    
    for ich,Ch in enumerate(sorted(DevDC)):
        if Ch=='Gate': continue
        
        chDC = DevDC[Ch]
        chDC['Ud']=np.zeros(len(chDC['Vds'])) 
        chDC['Rmin']=np.zeros(len(chDC['Vds'])) 
        chDC['Imin']=np.zeros(len(chDC['Vds']))
        chDC['GMax']=np.zeros(len(chDC['Vds'])) 
        chDC['Cicle']=ich  
        
        for ivd,Vds in enumerate(chDC['Vds']):                    
            Rds = Vds/chDC['Ids'][:,ivd]
            #Calc Dirac
            Imin = np.min(chDC['Ids'][:,ivd])
            Ud = chDC['Vgs'][np.argmin(chDC['Ids'][:,ivd])]
            Rmin = np.min(Rds)
            gm = np.diff(chDC['Ids'][:,ivd])/np.diff(chDC['Vgs'])
            GMax = np.max(gm)           
            
            chDC['Ud'][ivd] = Ud
            chDC['Rmin'][ivd] = Rmin
            chDC['Imin'][ivd] = Imin
            chDC['GMax'][ivd] = GMax

###############################################################################
#####
###############################################################################
def CalcGM (DevDC, DevAC=None, Order=10):
        
    for cName,DCdat in DevDC.items():
        if cName=='Gate': continue
    
        if DevAC:
            ACdat=DevAC[cName]           
        else:
            ACdat = None
        
        CalcGMCh(DCdat,ACdat,Order)
        
def CalcGMCh (DCdat,ACdat=None,Order=10):
        
    DCdat['IdsPoly'] = np.ones((Order+1,DCdat['Ids'].shape[1]))*np.NaN
    DCdat['GMPoly'] = np.ones((Order,DCdat['Ids'].shape[1]))*np.NaN 
    DCdat['Ud0'] = np.ones((DCdat['Ids'].shape[1]))*np.NaN
         
    if ACdat: 
        ACdat['GMPoly'] = np.ones((Order,DCdat['Ids'].shape[1]))*np.NaN
        ACdat['IdsPoly'] = np.ones((Order+1,DCdat['Ids'].shape[1]))*np.NaN
        ACdat['Ud0']=np.ones((DCdat['Ids'].shape[1]))*np.NaN
#            print(ACdat['Ud0']) 
    for ivd,Vds in enumerate(DCdat['Vds']):
        IdsPoly = np.polyfit(DCdat['Vgs'],DCdat['Ids'][:,ivd],Order)
        GMPoly = np.polyder(IdsPoly)
        
        DCdat['IdsPoly'][:,ivd] = IdsPoly
        DCdat['GMPoly'][:,ivd] = GMPoly
        
        IntVgs=np.linspace(DCdat['Vgs'][0],DCdat['Vgs'][-1],len(DCdat['Vgs'])*10000)
        DCdat['Ud0'][ivd]=IntVgs[np.argmin(np.polyval(IdsPoly,IntVgs))]
        
        if ACdat: 
            ACdat['GMPoly'][:,ivd] = GMPoly                   
            ACdat['IdsPoly'][:,ivd] = IdsPoly
            ACdat['Ud0'][ivd]=DCdat['Ud0'][ivd]

###############################################################################
#####
###############################################################################
def CheckIsOK (DevDC, DevAC=None, RdsRange = [400,10e3]):
    
    RdsMin = RdsRange[0]
    RdsMax = RdsRange[1]
    
    for ich,Ch in enumerate(sorted(DevDC)):
        if Ch=='Gate': continue
    
        chDC = DevDC[Ch]
        if DevAC:
            chAC=DevAC[Ch]
    
        for ivd,Vds in enumerate(chDC['Vds']):
            Rds = Vds/chDC['Ids'][:,ivd]
            
            if np.any([Rds<RdsMin,Rds>RdsMax]):
                print ('{} -- NOK -- {} {}'.format(Ch,np.min(Rds),np.max(Rds)))
                chDC['IsOK']=False
                if DevAC: chAC['IsOK']=False  
            else:
                chDC['IsOK']=True
                if DevAC: chAC['IsOK']=True
            
###############################################################################
#####
###############################################################################
def InterpolatePSD (DevACVals, Points=100):
                
    for ich,Ch in enumerate(sorted(DevACVals)):
        ch = DevACVals[Ch]
            
        Flin = ch['Fpsd'][1:]
        Flog = np.logspace(np.log10(Flin[0]),
                           np.log10(Flin[-1]),Points)    
        
        Flog = np.round(Flog,9)
        Flin = np.round(Flin,9)
        ch['Fpsd'] = Flog

        for Vds in ch['PSD']:
            PSDlin = ch['PSD'][Vds][:,1:]
            psd = interpolate.interp1d(Flin,PSDlin)
            ch['PSD'][Vds] = psd(Flog)            
            
###############################################################################
#####
###############################################################################
def FitACNoise(Dev, Fmin=None, Fmax=None, IsOkFilt=True):

    for ChName, ChDat in Dev.items():
        if ChName=='Gate': continue            
        if IsOkFilt:
            if not ChDat['IsOK']:continue

        nVgs = len(ChDat['Vgs'])
        nVds = len(ChDat['Vds'])

        if not 'NoA' in ChDat:
            ChDat['NoA'] = np.ones((nVgs,nVds))*np.NaN

        if not 'NoB' in ChDat:
            ChDat['NoB'] = np.ones((nVgs,nVds))*np.NaN

        if not 'FitErrA' in ChDat:
            ChDat['FitErrA'] = np.ones((nVgs,nVds))*np.NaN

        if not 'FitErrB' in ChDat:
            ChDat['FitErrB'] = np.ones((nVgs,nVds))*np.NaN

       
        for ivd in range(nVds):
            for ivg in range(nVgs):
                psd = ChDat['PSD']['Vd{}'.format(ivd)][ivg,:]
                Fpsd = ChDat['Fpsd']
                
                if np.any(np.isnan(psd)):continue
            
                try:
                    a,b, err = noise.FitNoise(Fpsd,psd,Fmin=Fmin,Fmax=Fmax)
                
                    ChDat['NoA'][ivg,ivd] = a
                    ChDat['NoB'][ivg,ivd] = b
                    ChDat['FitErrA'][ivg,ivd] = err[0]
                    ChDat['FitErrB'][ivg,ivd] = err[1]
                    
                except:
                    print ("Unexpected error:", sys.exc_info()[0])
                    print ('Channel Name ', ChName, 'Vd ', ivd, 'Vg ',ivg)

###############################################################################
#####
###############################################################################
def CalcNoiseIrms(Dev,Fmin=None,Fmax=None,IsOkFilt=True):
    
    for ChName, ChDat in Dev.items():
        if ChName=='Gate': continue            
        CalcNoiseIrmsCh(ChDat,Fmin,Fmax,IsOkFilt)    
    
def CalcNoiseIrmsCh(ChDat,Fmin=None,Fmax=None,IsOkFilt=True):
    if IsOkFilt:
        if not ChDat['IsOK']:return
  
    nVgs = len(ChDat['Vgs'])
    nVds = len(ChDat['Vds'])             
    
    ChDat['Irms'] = np.ones((nVgs,nVds))*np.NaN
     
    Fpsd = ChDat['Fpsd']        
         
    for ivd in range(nVds):
       for ivg in range(nVgs):            
           PSD = ChDat['PSD']['Vd{}'.format(ivd)][ivg,:]               
           ChDat['Irms'][ivg,ivd] = noise.PSDintegral(Fpsd,PSD,
                                                    Fmin=Fmin,Fmax=Fmax)    
     
    

