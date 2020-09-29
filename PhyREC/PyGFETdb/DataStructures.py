# -*- coding: utf-8 -*-
"""
@author: Anton Guimerà
@version: 0.1b

Revsion history

- 0.1b -- First version

TODO implement graph for selected channels

"""
import numpy as np
import matplotlib.pyplot as plt
import datetime
#from PlotData import PlotDC, PlotAC
import deepdish as dd
import os.path
#from DaqInterface import *

##############################################################################
#### Measure DC and AC characterization
##############################################################################


#def MeasDCAC(ACchar, DCchar, nVgs, nVds, ChNames, nACvgs,
#             Axes=None, DevDCVals=None, DevACVals=None):
#    
#    if not DevDCVals:    
#        DevDCVals = InitDCRecord(nVds=nVds, nVgs=nVgs, ChNames=ChNames)
#    
#    if not DevACVals:            
#        Fpsd = np.fft.rfftfreq(ACchar.PSDnFFT,1/ACchar.PSDFs)
#        Fgm = ACchar.SigFreqs
#        DevACVals = InitACRecord(nVds=nVds, nVgs=nVgs[nACvgs], nFgm=Fgm, nFpsd=Fpsd, ChNames=ChNames)
#    
#    ##############################################################################
#    ##### Main loop
#    ##############################################################################
#    iACvg = -1
#    for Vgs,ivg in zip(nVgs,range(len(nVgs))):
#        if np.in1d(ivg,nACvgs):
#            iACvg +=1
#        for Vds,ivd in zip(nVds,range(len(nVds))):
#            ####################################################
#            ######### Acquire DC values
#            nIds,Ig,data = DCchar.GetBiasCurrent(Vds, Vgs)
#    
#            ######### Fill data record
#            DevDCVals['Gate']['Ig'][ivg,ivd] = Ig
#    
#            for Ch,Ids in zip(sorted(DevDCVals),nIds):
#                DevDCVals[Ch]['Ids'][ivg,ivd]=Ids
#    
#            ######### Plot data record
#            PlotDC(DevDCVals,Axes,legend=True)
#            plt.pause(0.1)
#    
#            if np.in1d(ivg,nACvgs):
#                ####################################################
#                ######### Acquire AC values
#                psd,Fpsd,data = ACchar.GetPSD()
#                Gm, Fgm, OutSig = ACchar.GetBode()
#                
#        
#                ######### Fill data record
#                DevACVals[Ch]['Fpsd']=Fpsd
#                DevACVals[Ch]['Fgm']=Fgm
#        
#                for Ch in DevACVals:
#                    DevACVals[Ch]['PSD']['Vd{}'.format(ivd)][iACvg,:]=psd[:,ChNames[Ch]]
#                    DevACVals[Ch]['gm']['Vd{}'.format(ivd)][iACvg,:]=Gm[:,ChNames[Ch]]
#        
#                ######### Plot data record
#                PlotAC(DevACVals,Axes,iVds=ivd,iVgs=iACvg,legend=True)
#                plt.pause(0.1)
#
#    return DevDCVals, DevACVals
    
##############################################################################
#### Measure DCcharacterization
##############################################################################
#def MeasDC(DCchar, nVgs, nVds, ChNames, AxsDC=None, DevDCVals=None):
#    
#    if not DevDCVals:
#        DevDCVals = InitDCRecord(nVds=nVds, nVgs=nVgs, ChNames=ChNames)
#    
#    for Vgs,ivg in zip(nVgs,range(len(nVgs))):
#        for Vds,ivd in zip(nVds,range(len(nVds))):
#            ####################################################
#            ######### Acquire DC values
#            nIds,Ig,data = DCchar.GetBiasCurrent(Vds, Vgs)
#    
#            ######### Fill data record
#            DevDCVals['Gate']['Ig'][ivg,ivd] = Ig
#    
#            for Ch,Ids in zip(sorted(DevDCVals),nIds):
#                DevDCVals[Ch]['Ids'][ivg,ivd]=Ids
#    
#            ######### Plot data record
#            if AxsDC:
#                PlotDC(DevDCVals,AxsDC,legend=True)
#                plt.pause(0.1)
#
#    return DevDCVals

###############################################################################
#####
###############################################################################
def InitDCRecord(nVds, nVgs, ChNames, Gate=True):

    Time = datetime.datetime.now()
    DevDCVals={}
    for Ch in ChNames:
        DCVals={'Ids':np.ones((len(nVgs),len(nVds)))*np.NaN,
                'Vds':nVds,
                'Vgs':nVgs,
                'ChName':Ch,
                'Name':Ch,
                'DateTime':Time}
        DevDCVals[Ch]=DCVals

    if Gate:
        GateDCVals = {'Ig':np.ones((len(nVgs),len(nVds)))*np.NaN,
                    'Vds':nVds,
                    'Vgs':nVgs,
                    'ChName':'Gate',
                    'Name':'Gate',
                    'DateTime':Time}
        DevDCVals['Gate']=GateDCVals

    return DevDCVals

###############################################################################
#####
###############################################################################
def InitACRecord(nVds, nVgs, nFgm, nFpsd, ChNames):

    Time = datetime.datetime.now()
    DevACVals={}
    for Ch in ChNames:
        noise = {}
        gm = {}
        for i in range(nVds.size):
            noise['Vd{}'.format(i)] = np.ones((len(nVgs),nFpsd.size))*np.NaN
            gm['Vd{}'.format(i)] = np.ones((len(nVgs),nFgm.size))*np.NaN*np.complex(1)

        ACVals={'PSD':noise,
                'gm':gm,
                'Vgs':nVgs,
                'Vds':nVds,
                'Fpsd':nFpsd,
                'Fgm':nFgm,
                'ChName':Ch,
                'Name':Ch,
                'DateTime':Time}
        DevACVals[Ch]=ACVals

    return DevACVals

###############################################################################
#####
###############################################################################
def LoadDataFromFile (FileName): # TODO check the dictionary order DC, AC
    DataIn = dd.io.load(FileName)

    if type(DataIn) == dict:
        DevDCVals = DataIn
        DevACVals = None

    if type(DataIn) == tuple:
        DevDCVals = DataIn[0]
        DevACVals = None
        if len(DataIn) > 1:
            DevACVals = DataIn[1]

#    if 'Gate' not in DevDCVals.keys():
#        if 'Ids' not in  DevDCVals[DevDCVals.keys()[0]]:
#            print('Waring Data DC AC changed')
#            tmp1 = DevDCVals.copy()          
#            tmp2 = DevACVals.copy()          
#            DevACVals = tmp1         
#            DevDCVals = tmp2
        
    if DevACVals:
        for ch in DevACVals:
            if 'DateTime' not in DevACVals[ch]:
                print(ch, ' Waring Date not detected')
                DevACVals[ch]['DateTime'] = datetime.datetime.fromtimestamp(os.path.getmtime(FileName)) 
                
    if DevDCVals:
        for ch in DevDCVals:            
            if 'DateTime' not in DevDCVals[ch]:
                print(ch, ' Waring Date not detected')
                DevDCVals[ch]['DateTime'] = datetime.datetime.fromtimestamp(os.path.getmtime(FileName))  
                   
    return DevDCVals,DevACVals
    
###############################################################################
#####
###############################################################################
