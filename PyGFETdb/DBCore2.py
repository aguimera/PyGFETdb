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
        da = f.decrypt(open('Connection', 'rb').read())
        self.connection = pickle.loads(da)
        self._DEBUG = False

    def CreateQueryConditions(self, Conditions):
        Cond = []
        values = []
        for c in Conditions:
            cc = " %s OR ".join([c]*len(Conditions[c]))
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
            sInd = query.find('%s', sInd+1)

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
        return Res['result']

    def MultiSelect(self, Table, Conditions, FieldsOut, Order=None):
        Fields = ' , '.join(FieldsOut)
        Cond = ' %s AND '.join(Conditions.keys())
        sTable = ' JOIN '.join(Table)
        query = "SELECT {} FROM {} WHERE {} %s".format(Fields, sTable, Cond)
        if Order:
            query = '{} ORDER BY {}'.format(query, Order)

        return self._execute(query, list(Conditions.values()))

    def GetId(self, Table, Value, Field='Name', NewVals=None):
        query = "SELECT * FROM {} WHERE {} = %s".format(Table, Field)
        Res = self._execute(query, (Value,))

        if len(Res) == 0:
            if NewVals:
                ID = self.NewRow(Table=Table, Fields=NewVals)
            else:
                ID = None
        elif len(Res) == 1:
            ID = Res[0]['id{}'.format(Table)].decode()
        elif len(Res) > 1:
            print('Warning', query, Res)
            ID = Res[0]['id{}'.format(Table)].decode()

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
        return int(Res[0]['id'].decode())

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
        return int(Res[0]['id'].decode())

    def GetCharactInfo(self, Table, Conditions, Output):
        Out = ','.join(Output)
        Cond, Values = self.CreateQueryConditions(Conditions)

        query = """ SELECT distinct {}
                    FROM {}
                    INNER JOIN Trts ON Trts.idTrts = {}.Trt_id
                    INNER JOIN TrtTypes ON TrtTypes.idTrtTypes = Trts.Type_id
                    INNER JOIN Devices ON Devices.idDevices = Trts.Device_id
                    INNER JOIN Wafers ON Wafers.idWafers = Devices.Wafer_id
                    INNER JOIN VTrts ON VTrts.idTrts = Trts.idTrts
                    WHERE {}
                    ORDER BY {}.MeasDate; """.format(Out,
                                                     Table,
                                                     Table,
                                                     Cond,
                                                     Table)

        return self._execute(query, Values)

    def GetCharactFromId(self, Table, Id, GetGate=False):
        DataF = '{}.Data'.format(Table)
        Output = ['Trts.Name',
                  'Devices.Name',
                  'Devices.Comments',
                  'Wafers.Name',
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
                  DataF,
                  '{}.Ph'.format(Table),
                  '{}.Solution'.format(Table),
                  '{}.IonStrength'.format(Table),
                  '{}.IsCmp'.format(Table),
                  '{}.FuncStep'.format(Table),
                  '{}.AnalyteCon'.format(Table),
                  '{}.Comments'.format(Table)]
        Output2 = []
        OutF = []
        for o in Output:
            OutF.append(o.replace('.', '_'))
            Output2.append('{} AS {}'.format(o, o.replace('.', '_')))
        Out = ','.join(Output2)

        query = """ SELECT {}
                    FROM {}
                    INNER JOIN Trts ON Trts.idTrts = {}.Trt_id
                    INNER JOIN TrtTypes ON TrtTypes.idTrtTypes = Trts.Type_id
                    INNER JOIN Devices ON Devices.idDevices = Trts.Device_id
                    INNER JOIN Wafers ON Wafers.idWafers = Devices.Wafer_id
                    INNER JOIN VTrts ON VTrts.idTrts = Trts.idTrts
                    WHERE  {}.id{} = %s""".format(Out,
                                                  Table,
                                                  Table,
                                                  Table,
                                                  Table)
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
                fn[fp[1]] = Res[f].decode()
                Dat['Info'] = fn
            else:
                fn = Dat.get(fp[0], {})
                fn[fp[1]] = Res[f].decode()
                Dat[fp[0]] = fn

        Dat['Name'] = Res['Trts_Name'].decode()
        return Dat


        # GidF = None
        # if GetGate:
        #     if Table == 'DCcharacts':
        #         GidF = 'DCcharacts.Gate_id'
        #         Output.append(GidF)
        #     else:
        #         DCidF = 'ACcharacts.DC_id'
        #         Output.append(DCidF)
        # if GetGate:
        #     idg = None
        #     if GidF is None:
        #         Rows = self.MultiSelect(Table=('DCcharacts', ),
        #                                 Conditions={'idDCcharacts=': (re[DCidF], )},
        #                                 FieldsOut=('Gate_id', ))
        #         if len(Rows) > 0:
        #             idg = Rows[0][0]
        #     else:
        #         idg = re[GidF]

        #     if idg is not None:
        #         Data[Trt][cy]['Gate'] = self.GetGateFromId(idg)


    def GetData2(self, Conditions, Table, Last=True, GetGate=False):
        Output = ('{}.id{}'.format(Table, Table,),
                  '{}.MeasDate'.format(Table),
                  'Trts.Name AS TrtName')

        Res = self.GetCharactInfo(Table=Table,
                                     Conditions=Conditions,
                                     Output=Output)

        if len(Res) == 0:
            return []

        pdSers = []
        for r in Res:
            for k, v in r.items():
                r[k] = v.decode()
            pdSers.append(pd.Series(r))

        DTypes = {'id{}'.format(Table): int,
                  'MeasDate': 'datetime64[ns]',
                  'TrtName': str}
        dfr = pd.concat(pdSers, axis=1).transpose().astype(DTypes)

        if Last:
            gn = dfr.groupby('TrtName')
            ids = []
            for TrtName, dt in gn:
                dt.sort_values(by='MeasDate', ascending=False, inplace=True)
                ids.append(dt.iloc[0]['id{}'.format(Table)])
        else:
            ids = list(dfr['id{}'.format(Table)])

        Data = []
        for Id in ids:
            Data.append(self.GetCharactFromId(Table=Table, Id=Id))

        return Data
        # ids = []
        # Trts = []
        # for r in Res:
        #     ids.append(r['{}.id{}'.format(Table, Table)])
        #     Trts.append(r['Trts.Name'])

        # Data = self.GetCharactFromId(Table=Table,
        #                                 Ids=ids,
        #                                 Trts=set(Trts),
        #                                 Last=Last,
            #                             GetGate=GetGate)
        # return Data, list(set(Trts))

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
        Rows = self.MultiSelect(Table=('DCcharacts', ),
                                   Conditions={'Trt_id=': Trt_id,
                                               'MeasDate=': sMeasDate},
                                   FieldsOut=('Trt_id', 'idDCcharacts'))

        if len(Rows) == 0:  # New Record
            DCchatact_id = self.NewRow(Table='DCcharacts', Fields=NewData)
        else:
            print('WARNING EXISTS', Rows)
            if OverWrite:  # OverWrite
                DCchatact_id = Rows[0]['idDCcharacts'].decode()
                print('Overwrite Record id ', DCchatact_id)
                scond = 'idDCcharacts={}'.format(DCchatact_id)
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
                       'UpdateDate': TimeNow,
                       'DC_id': DCchatact_id,
                       'Data': pickle.dumps(ACVals),
                       'MeasDate': ACVals['DateTime']}
            if 'IsOK' in ACVals:
                NewData['IsOK'] = ACVals['IsOK']
            if OptFields:
                NewData.update(OptFields)

            # Check if exists
            sMeasDate = ACVals['DateTime'].strftime("%Y-%m-%d %H:%M:%S")
            Rows = self.MultiSelect(Table=('ACcharacts', ),
                                       Conditions={'Trt_id=': Trt_id,
                                                   'MeasDate=': sMeasDate},
                                       FieldsOut=('Trt_id', 'idACcharacts'))
            if len(Rows) == 0:  # New Record
                self.NewRow(Table='ACcharacts', Fields=NewData)
            else:
                print('WARNING EXISTS', Rows)
                if OverWrite:  # OverWrite
                    ACchatact_id = Rows[0]['idACcharacts'].decode()
                    print ('Overwrite Record id ', ACchatact_id)
                    self.UpdateRow(Table='ACcharacts',
                                      Fields=NewData,
                                      Condition=('idACcharacts=', ACchatact_id))


def Data2Pandas(Data):
    pdSeries = []
    
    for dd in Data:
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
            DataTypes[col] = float
        else:
            DataTypes[col] = 'category'
    
    DataTypes['CharCl'] = object
    DataTypes['IsOk'] = bool
    DataTypes['Date'] = 'datetime64[ns]'
    
    return dfRaw.astype(DataTypes)

if __name__ == '__main__':
    
    # f = Fernet(open('key.key', 'rb').read())
    # da = f.decrypt(open('Connection', 'rb').read())
    # connection = pickle.loads(da)
    MyDb = PyFETdb(open('key.key', 'rb').read())

#%% Test insert bin data
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
    
    
    #%% Test Read
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

    #%% Test Read
    
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

