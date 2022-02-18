# -*- coding: utf-8 -*-
"""
Template for uploading data into DB
works with GFETCharact software

###############################################################################
##### WARNING !!!!!!!!!!!!!!!!!!!
##### Naming rules B9355O15-T3-*****.h5
#####              [0] Wafer                             
#####                       [1] Device                             
#####                           [] anything else

@author: aguimera
"""

import PyGFETdb.DBCore as Pydb
import glob
from PyGFETdb.DBInterface import LoadPickleData

MyDB = Pydb.PyFETdb()


###############################################################################
# Files to upload

FileFilter ='B12744W3-Xip6NS-*.pkl'
FileNames = glob.glob(FileFilter)

###############################################################################
# Dataflieds o complete

Fields = {}
Fields['User'] = 'SB'

OptFields = {}
OptFields['Solution'] = '10mM PBS + 137 NaCl'
OptFields['FuncStep'] = 'AcTBA'
OptFields['Ph'] = 7.21
OptFields['AnalyteCon'] = ''
OptFields['IonStrength'] = 10e-3
OptFields['Comments'] = ''

###############################################################################
# Transitors type

TrtType = {'Name': 'RW50L50P3p0CS',
           'Length': 50e-6,
           'Width': 50e-6,
           'Pass': 2.5e-9,
           'Contact': 'Flat',
           'Shape': 'Rectangular'}

###############################################################################
# uploading
GateFields = {}
for ifile, filen in enumerate(FileNames):    
    ##### Naming rules B9355O15-T3-*****.h5
    #####              [0] Wafer                             
    #####                       [1] Device                             
    #####                           [] anything else
    filen.replace('\\', '/')
    fName = filen.split('/')[-1]
    fNameP = fName.split('-')
    Fields['Wafer'] = fNameP[0]
    Fields['Device'] = '{}-{}'.format(Fields['Wafer'], fNameP[1])
    OptFields['FileName'] = fName
    Fields['Gate_id'] = None

    print ('Load {} {} of {}'.format(fName, ifile, len(FileNames)))
    print ('Wafer ', Fields['Wafer'])
    print ('Device ', Fields['Device'])

    ######## Load Data    
    DictClDC, DictClAC, GateDict = LoadPickleData(filen)
    
    ######## Insert Gate and get GateID
    if GateDict is not None:
        GateFields['User'] = Fields['User']
        GateFields['Name'] = '{}-Gate'.format(Fields['Device'])
        GateFields['FileName'] = fName
        Fields['Gate_id'] = MyDB.InsertGateCharact(GateDict, GateFields)
    else:
        Fields['Gate_id'] = None

    print ('Gate ID -->> ', Fields['Gate_id'])
    ######## Insert Data
    for chn, datDC in DictClDC.items():    
        Fields['Trt'] = '{}-{}'.format(Fields['Device'], chn)   

        ##### Transistor Type definition
        Fields['TrtType'] = TrtType['Name']
        print('TRT {} Type {}'.format(Fields['Trt'], Fields['TrtType']))

        ###### Update DataBase
        DCd = datDC[0].RawData() 
                         
        if DictClAC is not None:
            ACd = DictClAC[chn][0].RawData()
            MyDB.InsertCharact(DCVals = DCd,
                               ACVals = ACd,
                               Fields = Fields,
                               OptFields = OptFields,
                               TrtTypeFields = TrtType)
        else:
            MyDB.InsertCharact(DCVals = DCd,
                               ACVals = None,
                               Fields = Fields,
                               OptFields = OptFields,
                               TrtTypeFields = TrtType)
            

del MyDB

