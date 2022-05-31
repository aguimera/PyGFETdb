#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 20:20:54 2022

@author: aguimera
"""

import datetime
from PyGFETdb.ws_json_class import ws_json
from PyGFETdb.DataClass import DataCharAC
import pickle
import sys
import pandas as pd
import numpy as np
import traceback
import base64
from cryptography.fernet import Fernet
import importlib


def WS_Types(name):
    if name.endswith('id'):
        return 'int'
    elif name == 'Data':
        return 'bin'
    else:
        return 'str'


class PyFETdb():

    def __init__(self, connection):
        f = Fernet(connection)
        # da = f.decrypt(open('./Connection', 'rb').read())
        da = f.decrypt(importlib.resources.read_binary('PyGFETdb', 'Connection'))
        self.connection = pickle.loads(da)
        self._DEBUG = True

    def CreateQueryConditions(self, Conditions):
        Cond = []
        values = []
        for c in Conditions:
            cc = " %s OR ".join([c] * len(Conditions[c]))
            Cond.append("({} %s)".format(cc))
            for v in Conditions[c]:
                values.append(v)
        sCond = ' AND '.join(Cond)
        return sCond, tuple(values)

    def ParseQuery(self, query, values):
        sInds = []
        sInd = query.find('%s')
        while sInd > 0:
            sInds.append(sInd)
            sInd = query.find('%s', sInd + 1)

        if len(sInds) > len(values):
            print("Error parsing query")
            traceback.print_exc()

        pVals = []
        for v in values:
            if type(v) == str:
                pVals.append("'{}'".format(v))
            elif type(v) == bytes:
                pVals.append("'{}'".format(base64.b64encode(v).decode()))
            else:
                pVals.append(str(v))

        return query % tuple(pVals)

    def _FormatOutputQuery(self, Output):
        Output2 = []
        OutF = []
        for o in Output:
            OutF.append(o.replace('.', '_'))
            Output2.append('{} AS {}'.format(o, o.replace('.', '_')))
        return OutF, ' , '.join(Output2)

    def _decode(self, Res):
        for r in Res:
            for k, v in r.items():
                if k.endswith('Data'):
                    continue
                r[k] = v.decode()
        return Res

    def _execute(self, query, values=None, oper='query'):
        if oper == 'query':
            q = self.ParseQuery(query, values)
            if self._DEBUG:
                print(q)
        else:
            q = query

        Res = ws_json.call(self.connection, q, oper)
        if Res['error']:
            # print(q)
            print('Error in call ', Res['errorCode'])
            return None

        return self._decode(Res['result'])

    def MultiSelect(self, Table, Conditions, Output, Order=None):
        OutF, Out = self._FormatOutputQuery(Output)
        Cond, Values = self.CreateQueryConditions(Conditions)
        query = "SELECT {} FROM {} WHERE {}".format(Out, Table, Cond)
        if Order:
            query = '{} ORDER BY {}'.format(query, Order)

        return self._execute(query, Values)

    def GetId(self, Table, Value, Field='Name', NewVals=None):
        query = "SELECT * FROM {} WHERE {} = %s".format(Table, Field)
        Res = self._execute(query, (Value,))

        if len(Res) == 0:
            if NewVals:
                ID = self.NewRow(Table=Table, Fields=NewVals)
            else:
                ID = None
        elif len(Res) == 1:
            ID = Res[0]['id{}'.format(Table)]
        elif len(Res) > 1:
            print('Warning', query, Res)
            ID = Res[0]['id{}'.format(Table)]

        return int(ID)

    def NewRow(self, Table, Fields):
        data = {}
        data['table'] = Table
        data['fields'] = {}
        for ic, (k, v) in enumerate(Fields.items()):
            data['fields'][ic] = {'name': k,
                                  'value': v,
                                  'type': WS_Types(k)}
        Res = self._execute(data, oper='insert')
        return int(Res[0]['id'])

    def UpdateRow(self, Table, Fields, Condition):
        data = {}
        data['table'] = Table
        data['fields'] = {}
        data['where'] = Condition
        for ic, (k, v) in enumerate(Fields.items()):
            data['fields'][ic] = {'name': k,
                                  'value': v,
                                  'type': WS_Types(k)}
        Res = self._execute(data, oper='update')
        return int(Res[0]['id'])

    def GetTrtsInfo(self, Conditions, Output=None):
        if Output is None:
            Output = ('Trts.idTrts',
                      'Trts.Name',
                      'TrtTypes.idTrtTypes',
                      'TrtTypes.Name',
                      'TrtTypes.Contact',
                      'TrtTypes.Width',
                      'TrtTypes.Length')

        OutF, Out = self._FormatOutputQuery(Output)
        Cond, Values = self.CreateQueryConditions(Conditions)
        query = """ SELECT distinct {}
                    FROM VTrts, Trts
                    INNER JOIN Devices ON Devices.idDevices = Trts.Device_id
                    INNER JOIN Wafers ON Wafers.idWafers = Devices.Wafer_id
                    INNER JOIN TrtTypes ON TrtTypes.idTrtTypes = Trts.Type_id
                    WHERE VTrts.idTrts = Trts.idTrts and  {} """.format(Out,
                                                                        Cond)

        return self._execute(query, Values)

    def GetCharactInfo(self, Table, Conditions, Output):
        OutF, Out = self._FormatOutputQuery(Output)
        Cond, Values = self.CreateQueryConditions(Conditions)

        query = """ SELECT distinct {Out}
                    FROM {Table}
                    INNER JOIN Trts ON Trts.idTrts = {Table}.Trt_id
                    INNER JOIN TrtTypes ON TrtTypes.idTrtTypes = Trts.Type_id
                    INNER JOIN Devices ON Devices.idDevices = Trts.Device_id
                    INNER JOIN Wafers ON Wafers.idWafers = Devices.Wafer_id
                    INNER JOIN VTrts ON VTrts.idTrts = Trts.idTrts
                    INNER JOIN Users ON {Table}.User_id = Users.idUsers
                    WHERE {Cond}
                    ORDER BY {Table}.MeasDate; """.format(Out=Out,
                                                          Table=Table,
                                                          Cond=Cond
                                                          )

        return self._execute(query, Values)

    def GetCharactFromId(self, Table, Id, GetGate=False):
        DataF = '{}.Data'.format(Table)
        Output = ['Trts.Name',
                  'Devices.Name',
                  'Devices.Comments',
                  'Wafers.Name',
                  'Wafers.Run',
                  'Wafers.Substrate',
                  'Wafers.Masks',
                  'Wafers.Graphene',
                  'Wafers.Comments',
                  'TrtTypes.Name',
                  'TrtTypes.Width',
                  'TrtTypes.Length',
                  'TrtTypes.Pass',
                  'TrtTypes.Area',
                  'TrtTypes.Contact',
                  'TrtTypes.Shape',
                  'Users.Name',
                  DataF,
                  '{}.Ph'.format(Table),
                  '{}.Solution'.format(Table),
                  '{}.IonStrength'.format(Table),
                  '{}.IsCmp'.format(Table),
                  '{}.FuncStep'.format(Table),
                  '{}.AnalyteCon'.format(Table),
                  '{}.Comments'.format(Table)]

        if GetGate:
            if Table == 'DCcharacts':
                GidF = 'DCcharacts.Gate_id'
                Output.append(GidF)
            else:
                GidF = None
                DCidF = 'ACcharacts.DC_id'
                Output.append(DCidF)

        OutF, Out = self._FormatOutputQuery(Output)

        query = """ SELECT {Out}
                    FROM {Table}
                    INNER JOIN Trts ON Trts.idTrts = {Table}.Trt_id
                    INNER JOIN TrtTypes ON TrtTypes.idTrtTypes = Trts.Type_id
                    INNER JOIN Devices ON Devices.idDevices = Trts.Device_id
                    INNER JOIN Wafers ON Wafers.idWafers = Devices.Wafer_id
                    INNER JOIN VTrts ON VTrts.idTrts = Trts.idTrts
                    INNER JOIN Users ON {Table}.User_id = Users.idUsers
                    WHERE  {Table}.id{Table} = %s """.format(Out=Out,
                                                             Table=Table,
                                                             )
        res = self._execute(query, (Id,))
        if len(res) == 0:
            print('Id Not found')
            return None
        if len(res) > 1:
            print('Warning duplicated ID !!!!')

        Res = res[0]

        Dat = {}
        for f in OutF:
            if f.endswith('Data'):
                Dat.update(pickle.loads(Res[f], encoding='latin'))
                continue
            fp = f.split('_')
            if fp[0] == Table:
                fn = Dat.get('Info', {})
                fn[fp[1]] = Res[f]
                Dat['Info'] = fn
            elif fp[0] == 'Users':
                fn = Dat.get('Info', {})
                fn['User'+fp[1]] = Res[f]
                Dat['Info'] = fn
            else:
                fn = Dat.get(fp[0], {})
                fn[fp[1]] = Res[f]
                Dat[fp[0]] = fn

        Dat['Name'] = Res['Trts_Name']

        if GetGate:
            if Table == 'DCcharacts':                
                Dat['Gate'] = self.GetGateFromId(Res['DCcharacts_Gate_id'])
            elif Table == 'ACcharacts':
                cond = {'idDCcharacts = ': (Res['ACcharacts_DC_id'], )}
                Res = self.MultiSelect(Table='DCcharacts',
                                        Conditions=cond,
                                        Output=('Gate_id', ))

                Dat['Gate'] = self.GetGateFromId(Res[0]['Gate_id'])             

        return Dat

    def GetGateFromId(self, idg):
        Rows = self.MultiSelect(Table='Gcharacts',
                                Conditions={'idGcharacts=': (idg, )},
                                Output=('Data', ))
        if len(Rows):
            GateData = pickle.loads(Rows[0]['Data'], encoding='latin')
            Res = self.GetCharactInfo(Table='DCcharacts',
                                      Conditions={'DCcharacts.Gate_id = ': (idg, )},
                                      Output=('Trts.Name',))
            Trts = []
            for re in Res:
                Trts.append(re['Trts_Name'])

            GateData['GateTrts'] = Trts
            return GateData
        else:
            return None

    def GetData2(self, Conditions, Table, Last=True, GetGate=False):
        Output = ('id{}'.format(Table, Table, ),
                  'MeasDate'.format(Table),
                  'Trts.Name')

        Res = self.GetCharactInfo(Table=Table,
                                  Conditions=Conditions,
                                  Output=Output)

        if len(Res) == 0:
            return []

        pdSers = []
        for r in Res:
            pdSers.append(pd.Series(r))

        DTypes = {'id{}'.format(Table): int,
                  'MeasDate': 'datetime64[ns]',
                  'Trts_Name': str}
        dfr = pd.concat(pdSers, axis=1).transpose().astype(DTypes)

        if Last:
            gn = dfr.groupby('Trts_Name')
            ids = []
            for TrtName, dt in gn:
                dt.sort_values(by='MeasDate', ascending=False, inplace=True)
                ids.append(dt.iloc[0]['id{}'.format(Table)])
        else:
            ids = list(dfr['id{}'.format(Table)])

        Data = []
        for ic, Id in enumerate(ids):
            print("Downloading {} of {}".format(ic, len(ids)))
            Data.append(self.GetCharactFromId(Table=Table,
                                              Id=Id,
                                              GetGate=GetGate))

        return Data

    def InsertCharact(self, DCVals, Fields, ACVals=None, OptFields=None,
                      TrtTypeFields=None, OverWrite=True):

        Author = Fields['User']
        Wafer = Fields['Wafer']
        Device = Fields['Device']
        TrtType = Fields['TrtType']
        Trt = Fields['Trt']

        if not Wafer.replace('-', '').isalnum():
            print('ABORTED')
            print('Wafer name not valid characters')
            return
        if not Device.replace('-', '').isalnum():
            print('ABORTED')
            print('Device name not valid characters')
            return
        if not Trt.replace('-', '').isalnum():
            print('ABORTED')
            print('Trt name not valid characters')
            return
        if len(Wafer) > 12:
            print('ABORTED')
            print('Wafer name too long, 12 allowed')
            return
        if len(Device) > 20:
            print('ABORTED')
            print('Device name too long, 20 allowed')
            return
        if len(Trt) > 25:
            print('ABORTED')
            print('Trt name too long, 25 allowed')
            return

        # Set Get information in other tables
        Author_id = self.GetId(Table='Users',
                               Value=Author,
                               NewVals={'Name': Author})

        TrtType_id = self.GetId(Table='TrtTypes',
                                Value=TrtType,
                                NewVals=TrtTypeFields)

        Wafer_id = self.GetId(Table='Wafers',
                              Value=Wafer,
                              NewVals={'Name': Wafer})

        Device_id = self.GetId(Table='Devices',
                               Value=Device,
                               NewVals={'Name': Device,
                                        'Wafer_id': Wafer_id})

        Trt_id = self.GetId(Table='Trts',
                            Value=Trt,
                            NewVals={'Name': Trt,
                                     'Device_id': Device_id,
                                     'Type_id': TrtType_id})

        TimeNow = datetime.datetime.now()

        #######################################################################
        # Insert DC data
        #######################################################################
        # FILL NewData structure
        NewData = {'Trt_id': Trt_id,
                   'User_id': Author_id,
                   'Data': pickle.dumps(DCVals),
                   'MeasDate': str(DCVals['DateTime']),
                   'UpdateDate': str(TimeNow),
                   'Gate_id': Fields['Gate_id']}

        if 'IsOK' in DCVals:
            NewData['IsOK'] = DCVals['IsOK']
        if OptFields:
            NewData.update(OptFields)

        # Check if exists
        sMeasDate = DCVals['DateTime'].strftime("%Y-%m-%d %H:%M:%S")
        Rows = self.MultiSelect(Table='DCcharacts',
                                Conditions={'Trt_id=': (Trt_id, ),
                                            'MeasDate=': (sMeasDate, )},
                                Output=('Trt_id', 'idDCcharacts'))

        if len(Rows) == 0:  # New Record
            DCchatact_id = self.NewRow(Table='DCcharacts', Fields=NewData)
        else:
            print('WARNING EXISTS', Rows)
            if OverWrite:  # OverWrite
                DCchatact_id = Rows[0]['idDCcharacts']
                print('Overwrite Record id ', DCchatact_id)
                scond = 'idDCcharacts = {}'.format(DCchatact_id)
                self.UpdateRow(Table='DCcharacts',
                               Fields=NewData,
                               Condition=scond)

        NewData.pop('Gate_id')
        #######################################################################
        # Insert DC data
        #######################################################################
        if ACVals:
            # FILL NewData structure
            NewData = {'Trt_id': Trt_id,
                       'User_id': Author_id,
                       'UpdateDate': str(TimeNow),
                       'DC_id': DCchatact_id,
                       'Data': pickle.dumps(ACVals),
                       'MeasDate': str(ACVals['DateTime'])}
            if 'IsOK' in ACVals:
                NewData['IsOK'] = ACVals['IsOK']
            if OptFields:
                NewData.update(OptFields)

            # Check if exists
            sMeasDate = ACVals['DateTime'].strftime("%Y-%m-%d %H:%M:%S")
            Rows = self.MultiSelect(Table='ACcharacts',
                                    Conditions={'Trt_id=': (Trt_id,),
                                                'MeasDate=': (sMeasDate,)},
                                    Output=('Trt_id', 'idACcharacts'))
            if len(Rows) == 0:  # New Record
                self.NewRow(Table='ACcharacts', Fields=NewData)
            else:
                print('WARNING EXISTS', Rows)
                if OverWrite:  # OverWrite
                    ACchatact_id = Rows[0]['idACcharacts']
                    print('Overwrite Record id ', ACchatact_id)
                    scond = 'idACcharacts = {}'.format(ACchatact_id)
                    self.UpdateRow(Table='ACcharacts',
                                   Fields=NewData,
                                   Condition= scond)
                    

    def InsertGateCharact(self, GateData, Fields, OverWrite=True):
        # Set Get information in other tables
        Author = Fields['User']
        Author_id = self.GetId(Table='Users',
                               Value=Author,
                               NewVals={'Name': Author})
        Fields.pop('User')

        # Check if exists
        sMeasDate = GateData['DateTime'].strftime("%Y-%m-%d %H:%M:%S")
        Rows = self.MultiSelect(Table='Gcharacts',
                                Conditions={'Name=': (Fields['Name'],),
                                            'MeasDate=': (sMeasDate, )},
                                Output=('idGcharacts', ))

        NewData = {'Data': pickle.dumps(GateData),
                   'User_id': Author_id,
                   'MeasDate': str(GateData['DateTime']),
                   'UpdateDate': str(datetime.datetime.now())}
        NewData.update(Fields)

        if len(Rows) == 0:  # New Record
            Gate_id = self.NewRow(Table='Gcharacts', Fields=NewData)
        else:                        
            Gate_id = Rows[0]['idGcharacts']
            print ('Gate id foung {}'. format(Gate_id))
            if OverWrite:  # OverWrite                
                print ('WARNING EXISTS', Rows)
                print ('Overwrite Record id ', Gate_id)
                scond = 'idGcharacts = {}'.format(Gate_id)
                self.UpdateRow(Table='Gcharacts',
                               Fields=NewData,
                               Condition=scond)
        return Gate_id
    
    def InsertDataFrame(self, dfUpload):
        pass


def Data2Pandas(Data):
    pdSeries = []
    for ic, dd in enumerate(Data):
        print("Converting {} of {}".format(ic, len(Data)))
        pdser = {}
        d = DataCharAC(dd)
        pdser['Name'] = d.Name
        pdser['CharCl'] = d
        pdser['IsOk'] = d.IsOK
        pdser['ChName'] = d.ChName
        pdser['Date'] = d.GetDateTime()
        for k, v in d['Wafers'].items():
            if k == 'Name':
                pdser['Wafer'] = v
            else:
                pdser['Waf' + k] = v
        for k, v in d['Devices'].items():
            if k == 'Name':
                pdser['Device'] = v
            else:
                pdser['Dev' + k] = v
        for k, v in d['TrtTypes'].items():
            if k == 'Name':
                pdser['TrtType'] = v
            else:
                pdser[k] = v
        for k, v in d['Info'].items():
            if k == 'Gate_id':
                continue
            pdser[k] = v
        pdSeries.append(pd.Series(pdser))

    dfRaw = pd.concat(pdSeries, axis=1).transpose()
    DataTypes = {}
    for col in dfRaw.keys():
        if col in ('Width', 'Length', 'Pass', 'Area', 'Ph', 'IonStrength', 'AnalyteCon'):
            dfRaw[col] = pd.to_numeric(dfRaw[col], errors='coerce')
            DataTypes[col] = float
        else:
            DataTypes[col] = 'category'

    DataTypes['CharCl'] = object
    DataTypes['IsOk'] = bool
    DataTypes['Date'] = 'datetime64[ns]'

    return dfRaw.astype(DataTypes, errors='ignore')


if __name__ == '__main__':
    # f = Fernet(open('key.key', 'rb').read())
    # da = f.decrypt(open('Connection', 'rb').read())
    # connection = pickle.loads(da)
    MyDb = PyFETdb(open('key.key', 'rb').read())

# %% Test insert bin data
# OptFields = {}
# OptFields['Solution'] = 'TestWS'
# OptFields['FuncStep'] = 'TestWS'
# OptFields['Ph'] = 7.21

# by = pickle.dumps(OptFields,0)

# data = {}
# data['table'] = "DCcharacts"
# data['fields'] = {}
# data['fields'][0] = {}
# data['fields'][1] = {}
# data['fields'][2] = {}
# data['fields'][0]['name'] = "Trt_id"
# data['fields'][0]['value'] = 100
# data['fields'][0]['type'] = "int"
# data['fields'][1]['name'] = "Comments"
# data['fields'][1]['value'] = 'TestWS'
# data['fields'][1]['type'] = "str"
# data['fields'][2]['name'] = "Data"
# data['fields'][2]['value'] = by
# data['fields'][2]['type'] = "bin"

# res = ws_json.call(connection, data,"insert")
# '''
#  	Per un update, per exemple,
# data['where'] = "idDCcharacts=458586"
# res = ws_json_class.ws_json.call("pyfet", data,"update")
# '''
# try:
#  	error = res['error']
#  	if (error):
#          print("S'ha produït un error ")
#          print("\t" + str(res['errorCode']) + " : " + res['message'])
#  	else:
#          print("El resultat és: ")
#          print(res['result'])
# except:
#  	print("S'ha produït un error inesperat")
#  	sys.exit(0)


# qt = "SELECT * FROM DCcharacts WHERE Comments='TestWS'"

# Res = ws_json.call(connection, qt)


# %% Test Read
# DevicesList = ('B12744W3-Xip6NS', 'B12744W3-Xip5NS',
#                 )

# CharTable = 'DCcharacts'
# Conditions = {'Devices.name = ': DevicesList,
#               }

# GroupBase = {'Conditions': Conditions,
#               'Table': CharTable,
#               'Last': True,
#               'GetGate': True,
#               }


# res = MyDb.GetData2(**GroupBase)
# dfRaw = Data2Pandas(res)

# %% Test Read

# Fields = {}
# Fields['User'] = 'SB'

# OptFields = {}
# OptFields['Solution'] = 'TestWS'
# OptFields['FuncStep'] = 'TestWS'
# OptFields['Ph'] = 7.21
# OptFields['AnalyteCon'] = 'tt'
# OptFields['IonStrength'] = 10e-3
# OptFields['Comments'] = 'tt'

# Fields['Wafer'] = 'TestWS'
# Fields['Device'] = 'TestWS-Dev'
# Fields['Trt'] = 'TestWS-Dev'
# Fields['TrtType'] = 'RW50L50P3p0CS'
# Fields['Gate_id'] = None

# DCd = res[0]

# MyDb.InsertCharact(DCVals=DCd,
#                       ACVals=None,
#                       Fields=Fields,
#                       OptFields=OptFields,
#                       TrtTypeFields=None)
