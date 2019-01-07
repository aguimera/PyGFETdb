# -*- coding: utf-8 -*-
"""
@author: Anton Guimerà
@version: 0.1b

Revsion history

- 0.1b -- First version

"""

import numpy as np
import scipy.optimize as optim
from scipy.integrate import simps

###############################################################################
#### !/f noise functions
###############################################################################
def Fnoise (f,a,b):
    '''
    return a/f^b
    '''
    return a/f**b

def FnoiseTh (f,a,b,c):
    '''
    return a/f^b+c
    '''
    return a/f**b+c

def LogFnoise (f,a,b):
    '''
    return b*f+a
    '''
    return b*f+a

###############################################################################
#### Fitting functions
###############################################################################
def FitNoise(Freq, psd, Fmin=None, Fmax=None):   ### Funtions redefinition    
#    return FitFNoise (Freq, psd, Fmin=Fmin, Fmax=Fmax)
    return FitLogFnoise (Freq, psd, Fmin=Fmin, Fmax=Fmax)

####### Linear fitting
def FitFNoise(Freq, psd, Fmin=None, Fmax=None):
    '''
    return a, b, pcov
    '''
    Inds = CalcFreqIndexes(Freq, Fmin, Fmax)    
    poptV, pcov = optim.curve_fit(Fnoise, Freq[Inds], psd[Inds])
    a = poptV[0]
    b = poptV[1]

    return a, b, np.sqrt(np.diag(pcov))

####### Log fitting
def FitLogFnoise(Freq, psd, Fmin=None, Fmax=None):
    
    Inds = CalcFreqIndexes(Freq, Fmin, Fmax)
    poptV, pcov = optim.curve_fit(LogFnoise, np.log10(Freq[Inds]), 
                                  np.log10(psd[Inds]))

    a = 10 ** poptV[0] 
    b = - poptV[1]
    
    return a, b, np.sqrt(np.diag(pcov))

###############################################################################
#### !/f noise function
###############################################################################
def PSDintegral(Freq,psd, Fmin=1, Fmax=5e3):
    '''
    return Irms
    '''
    Inds = CalcFreqIndexes(Freq, Fmin, Fmax)
    Irms = simps(psd[Inds],Freq[Inds])
 
    return np.sqrt(Irms)

###############################################################################
#### Calc indexes for Fmin Fmax
###############################################################################
def CalcFreqIndexes(Freq,Fmin=None,Fmax=None):    

    if Fmin:
        return np.where(Freq>Fmin)
    if Fmax:
        return np.where(Freq<Fmax)
    if Fmin and Fmax:
        return np.where(np.logical_and(Freq>Fmin, Freq<Fmax)) 
    
    return range(len(Freq))[1:]

