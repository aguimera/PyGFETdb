# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 23:15:19 2016

@author: aguimera
"""

import datetime
import pymysql
import pickle
from PyGFETdb.DB import *
import sys

class PyFETdb():
    PrintQuery = False

    def __init__(self, host=None, user=None, passwd=None, db=None, Update=True):

        if host is None:
            self.db = pymysql.connect(host=DBhost,
                                      user=DBuser,
                                      passwd=DBpasswd,
                                      db=DBdb)
        else:
            self.db = pymysql.connect(host=host,
                                      user=user,
                                      passwd=passwd,
                                      db=db)

    def __del__(self):
        self.db.close()

    def _execute(self, query, values, LastRowID=False):
        cur = self.db.cursor()
        if self.PrintQuery:
            print (query, values)
        try:
            cur.execute(query, values)
        except:
            print ('Error')
            self.db.connect()
            cur = self.db.cursor()
            cur.execute(query, values)

        if LastRowID:
            ret = cur.lastrowid
        else:
            ret = cur.fetchall()

        cur.close()
        return ret

    def _DecodeData(self, DataF):
        if sys.version_info[0] == 2:
            return pickle.loads(DataF)
        else:
            try:
                bDict = pickle.loads(DataF, encoding='latin')
                return bDict
            except:
                Dict = {}
                for k, v in bDict.items():
                    if isinstance(v, dict):
                        dv = {}
                        for kv, vv in v.items():
                            dv[kv.decode('utf-8')] = vv
                    else:
                        dv = v
                    Dict[k.decode('utf-8')] = dv
                return Dict



    def GetId(self, Table, Value, Field='Name', NewVals=None):

        query = "SELECT * FROM {} WHERE {} = %s".format(Table, Field)
        Res = self._execute(query, (Value,))

        if len(Res) == 0:
            if NewVals:
                ID = self.NewRow(Table=Table, Fields=NewVals)
            else:
                ID = None
        if len(Res) == 1:
            ID = Res[0][0]

        if len(Res) > 1:
            print ('Warning', query, Res)
            ID = Res[0][0]

        return ID

    def NewRow(self, Table, Fields):
        colums = ' , '.join(Fields.keys())
        places = ' , '.join(['%s']*len(Fields))
        query = "INSERT INTO {}({}) VALUES ({})".format(Table,
                                                        colums,
                                                        places)
        return self._execute(query, list(Fields.values()), LastRowID=True)

    def MultiSelect(self, Table, Conditions, FieldsOut, Order=None):
        Fields = ' , '.join(FieldsOut)
        Cond = ' %s AND '.join(Conditions.keys())
        sTable = ' JOIN '.join(Table)
        query = "SELECT {} FROM {} WHERE {} %s".format(Fields, sTable, Cond)
        if Order:
            query = '{} ORDER BY {}'.format(query, Order)

        return self._execute(query, list(Conditions.values()))

    def UpdateRow(self, Table, Fields, Condition):
        colums = ' = %s , '.join(Fields.keys())
        colums = colums + ' = %s'
        query = "UPDATE {} SET {} WHERE {} %s".format(Table,
                                                      colums,
                                                      Condition[0])
        values = list(Fields.values())
        values.insert(len(values), Condition[1])
        self._execute(query, values)

    def InsertGateCharact(self, DCVals, Fields, OverWrite=True):
        # Set Get information in other tables
        Author = Fields['User']
        Author_id = self.GetId(Table='Users',
                               Value=Author,
                               NewVals={'Name': Author})
        Fields.pop('User')

        # Check if exists
        Rows = self.MultiSelect(Table=('Gcharacts', ),
                                Conditions={'Name=': Fields['Name'],
                                            'MeasDate=':DCVals['DateTime'].strftime("%Y-%m-%d %H:%M:%S")},
                                FieldsOut=('idGcharacts', ))

        NewData = {'Data': pickle.dumps(DCVals),
                   'User_id': Author_id,
                   'MeasDate': DCVals['DateTime'],
                   'UpdateDate': datetime.datetime.now()}
        NewData.update(Fields)

        if len(Rows) == 0:  # New Record
            Gate_id = self.NewRow(Table='Gcharacts', Fields=NewData)
        else:
            print ('WARNING EXISTS', Rows)
            if OverWrite:  # OverWrite
                Gate_id = Rows[0][0]
                print ('Overwrite Record id ', Gate_id)
                self.UpdateRow(Table='Gcharacts',
                               Fields=NewData,
                               Condition=('idGcharacts=', Gate_id))
        return Gate_id

    def InsertCharact(self, DCVals, Fields, ACVals=None, OptFields=None,
                      TrtTypeFields=None, OverWrite=True):

        Author = Fields['User']
        Wafer = Fields['Wafer']
        Device = Fields['Device']
        TrtType = Fields['TrtType']
        Trt = Fields['Trt']

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
                   'MeasDate': DCVals['DateTime'],
                   'UpdateDate': TimeNow,
                   'Gate_id': Fields['Gate_id']}

        if 'IsOK' in DCVals:
            NewData['IsOK'] = DCVals['IsOK']
        if OptFields:
            NewData.update(OptFields)

        # Check if exists
        Rows = self.MultiSelect(Table=('DCcharacts', ),
                                Conditions={'Trt_id=': Trt_id,
                                            'MeasDate=': DCVals['DateTime'].strftime("%Y-%m-%d %H:%M:%S")},
                                FieldsOut=('Trt_id', 'idDCcharacts'))
        if len(Rows) == 0:  # New Record
            DCchatact_id = self.NewRow(Table='DCcharacts', Fields=NewData)
        else:
            print ('WARNING EXISTS', Rows)
            if OverWrite:  # OverWrite
                DCchatact_id = Rows[0][1]
                print ('Overwrite Record id ', DCchatact_id)
                self.UpdateRow(Table='DCcharacts',
                               Fields=NewData,
                               Condition=('idDCcharacts=', DCchatact_id))

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
            Rows = self.MultiSelect(Table=('ACcharacts', ),
                                    Conditions={'Trt_id=': Trt_id,
                                                'MeasDate=': ACVals['DateTime'].strftime("%Y-%m-%d %H:%M:%S")},
                                    FieldsOut=('Trt_id', 'idACcharacts'))
            if len(Rows) == 0:  # New Record
                self.NewRow(Table='ACcharacts', Fields=NewData)
            else:
                print ('WARNING EXISTS', Rows)
                if OverWrite:  # OverWrite
                    ACchatact_id = Rows[0][1]
                    print ('Overwrite Record id ', ACchatact_id)
                    self.UpdateRow(Table='ACcharacts',
                                   Fields=NewData,
                                   Condition=('idACcharacts=', ACchatact_id))
        # Finishing
        self.db.commit()

    def CreateQueryConditions(self, Conditions):
        Cond = []
        values = []
        for c in Conditions:
            cc = ' %s OR '.join([c]*len(Conditions[c]))
            Cond.append('({} %s)'.format(cc))
            for v in Conditions[c]:
                values.append(v)
        sCond = ' AND '.join(Cond)
        return sCond, values

    def FindFillOutput(self, query, values, Output):
        res = self._execute(query, values)
        ret = []
        for ir, r in enumerate(res):
            do = {}
            for io, of in enumerate(Output):
                do[of] = r[io]
            ret.insert(len(ret), do)
        return ret

    def GetTrtsInfo(self, Conditions, Output=None):
        if Output is None:
            Output = ('Trts.idTrts',
                      'Trts.Name',
                      'TrtTypes.idTrtTypes',
                      'TrtTypes.Name',
                      'TrtTypes.Contact',
                      'TrtTypes.Width',
                      'TrtTypes.Length')

        Out = ','.join(Output)
        Cond, Values = self.CreateQueryConditions(Conditions)
        query = """ SELECT {}
                    FROM VTrts, Trts
                    INNER JOIN Devices ON Devices.idDevices = Trts.Device_id
                    INNER JOIN Wafers ON Wafers.idWafers = Devices.Wafer_id
                    INNER JOIN TrtTypes ON TrtTypes.idTrtTypes = Trts.Type_id
                    WHERE VTrts.idTrts = Trts.idTrts and  {} """.format(Out,
                                                                        Cond)

        return self.FindFillOutput(query, Values, Output)

    def GetDevicesInfo(self, Conditions, Output=None):
        if Output is None:
            Output = ('Devices.idDevices',
                      'Devices.Name')
        Out = ','.join(Output)

        Cond, Values = self.CreateQueryConditions(Conditions)
        query = """ SELECT {}
                    FROM Devices
                    INNER JOIN Wafers ON Wafers.idWafers = Devices.Wafer_id
                    WHERE {} """.format(Out, Cond)

        return self.FindFillOutput(query, Values, Output)

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

        return self.FindFillOutput(query, Values, Output)

    def GetCharactFromId(self, Table, Ids, Trts, Last=False, GetGate=False):

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
                  '{}.FuncStep'.format(Table),
                  '{}.AnalyteCon'.format(Table),
                  '{}.Comments'.format(Table)]

        GidF = None
        if GetGate:
            if Table == 'DCcharacts':
                GidF = 'DCcharacts.Gate_id'
                Output.append(GidF)
            else:
                DCidF = 'ACcharacts.DC_id'
                Output.append(DCidF)

        Cond = {}
        idf = '{}.id{}='.format(Table, Table)
        Cond[idf] = []
        for s in set(Ids):
            Cond[idf].append(s)

        Data = {}
        for Trt in Trts:
            Cond['Trts.Name='] = (Trt, )
            Res = self.GetCharactInfo(Table, Cond, Output)

            Data[Trt] = {}

            if Last:
                Re = (Res[-1],)
            else:
                Re = Res

            for cy, re in enumerate(Re):
                cy = 'Cy{0:03d}'.format(cy)                
                Data[Trt][cy] = self._DecodeData(re[DataF])
#                Data[Trt][cy]['Ph'] = pickle.loads(re[DataPh])
#                Data[Trt][cy]['FuncCond'] = pickle.loads(re[DataPh])

                if GetGate:
                    idg = None
                    if GidF is None:
                        Rows = self.MultiSelect(Table=('DCcharacts', ),
                                                Conditions={'idDCcharacts=': (re[DCidF], )},
                                                FieldsOut=('Gate_id', ))
                        if len(Rows) > 0:
                            idg = Rows[0][0]
                    else:
                        idg = re[GidF]

                    if idg is not None:
                        Data[Trt][cy]['Gate'] = self.GetGateFromId(idg)

                for f in Output:
                    if f.endswith('Data'):
                        continue
                    fp = f.split('.')
                    if fp[0] == 'Trts':
                        Data[Trt][cy][fp[1]] = re[f]
                    elif fp[0] == Table:
                        fn = Data[Trt][cy].get('Info', {})
                        fn[fp[1]] = re[f]
                        Data[Trt][cy]['Info'] = fn
                    else:
                        fn = Data[Trt][cy].get(fp[0], {})
                        fn[fp[1]] = re[f]
                        Data[Trt][cy][fp[0]] = fn

        return Data

    def GetGateFromId(self, idg):
        Rows = self.MultiSelect(Table=('Gcharacts', ),
                                Conditions={'idGcharacts=': (idg, )},
                                FieldsOut=('Data', ))
        if len(Rows) > 0:
            GateData = self._DecodeData(Rows[0][0])
            Res = self.GetCharactInfo(Table='DCcharacts',
                                      Conditions={'DCcharacts.Gate_id=': (idg, )},
                                      Output=('Trts.Name',))
            Trts = []
            for re in Res:
                Trts.append(re['Trts.Name'])

            GateData['GateTrts'] = Trts
            return GateData
        else:
            return None

    def DeleteCharact(self, Table, Ids):
        query = "DELETE FROM {} WHERE id{}=%s".format(Table, Table)
        for Id in Ids:
            self._execute(query, Id)
        self.db.commit()

    def GetData2(self, Conditions, Table, Last=True, GetGate=False):
        Output = ('{}.id{}'.format(Table, Table,),
                  'Trts.Name')

        Res = self.GetCharactInfo(Table=Table,
                                  Conditions=Conditions,
                                  Output=Output)
        ids = []
        Trts = []
        for r in Res:
            ids.append(r['{}.id{}'.format(Table, Table)])
            Trts.append(r['Trts.Name'])

        Data = self.GetCharactFromId(Table=Table,
                                     Ids=ids,
                                     Trts=set(Trts),
                                     Last=Last,
                                     GetGate=GetGate)
        return Data, list(set(Trts))

    def GetTrtCharact2(self, Table, TrtId, TrtName, Last=False):

        cond = {'{}.Trt_id='.format(Table): TrtId}

        Res = self.MultiSelect(Table=(Table,),
                               Conditions=cond,
                               FieldsOut=('{}.Data'.format(Table),
                                          '{}.MeasDate'.format(Table)),
                               Order='{}.MeasDate'.format(Table))

        if len(Res) == 0:
            return None

        Vals = {}
        if Last:
            Vals['Cy{0:03d}'.format(0)] = self._DecodeData(Res[-1][0])
            if TrtName:
                Vals['Cy{0:03d}'.format(0)]['Name'] = TrtName
        else:
            for cy, re in enumerate(Res):
                Vals['Cy{0:03d}'.format(cy)] = self._DecodeData(re[0])
                if TrtName:
                    Vals['Cy{0:03d}'.format(cy)]['Name'] = TrtName

        return Vals

    def GetData(self, Conditions, DC=True, AC=True, Last=False,
                Date=None, IsCmp=None):

        Output = ('Trts.idTrts', 'Trts.Name',
                  'TrtTypes.idTrtTypes', 'TrtTypes.Name',
                  'TrtTypes.Contact', 'TrtTypes.Width', 'TrtTypes.Length',
                  'Devices.Name', 'Wafers.Name')

        Trts = self.GetTrtsInfo(Conditions=Conditions, Output=Output)

        DataDC = {}
        DataAC = {}
        TrtsOut = []
        for T in Trts:
            if DC:
                DCVals = self.GetTrtCharact(Table='DCcharacts',
                                            TrtId=T['Trts.idTrts'],
                                            TrtName=T['Trts.Name'],
                                            Last=Last,
                                            Date=Date,
                                            IsCmp=IsCmp)
                if DCVals is not None:
                    DataDC[T['Trts.Name']] = DCVals
                    TrtsOut.append(T)
            if AC:
                ACVals = self.GetTrtCharact(Table='ACcharacts',
                                            TrtId=T['Trts.idTrts'],
                                            TrtName=T['Trts.Name'],
                                            Last=Last,
                                            Date=Date,
                                            IsCmp=IsCmp)
                if ACVals is not None:
                    TrtsOut.append(T)
                    DataAC[T['Trts.Name']] = ACVals

        return DataDC, DataAC, TrtsOut

    def GetTrtCharact(self, Table, TrtId, TrtName=None,
                      Last=False, Date=None, IsCmp=None):

        cond = {'{}.Trt_id='.format(Table): TrtId}
        if Date is not None:
            cond['{}.MeasDate>'.format(Table)] = Date[0]
            cond['{}.MeasDate<'.format(Table)] = Date[1]

        if IsCmp is not None:
            cond['{}.IsCmp='.format(Table)] = IsCmp

        Res = self.MultiSelect(Table=(Table,),
                               Conditions=cond,
                               FieldsOut=('{}.Data'.format(Table),
                                          '{}.MeasDate'.format(Table)),
                               Order='{}.MeasDate'.format(Table))

        if len(Res) == 0:
            return None

        Vals = {}
        if Last:
            Vals['Cy{0:03d}'.format(0)] = self._DecodeData(Res[-1][0])
            if TrtName:
                Vals['Cy{0:03d}'.format(0)]['Name'] = TrtName
        else:
            for cy, re in enumerate(Res):
                Vals['Cy{0:03d}'.format(cy)] = self._DecodeData(re[0])
                if TrtName:
                    Vals['Cy{0:03d}'.format(cy)]['Name'] = TrtName

        return Vals


#    
#    def GetTrts(self, Conditions, Output = None):
#        
#        if not Output:
#            Output = ('Trts.idTrts', 'Trts.Name', 
#                      'TrtTypes.idTrtTypes', 'TrtTypes.Name',  
#                      'TrtTypes.Contact', 'TrtTypes.Width', 'TrtTypes.Length')
#        Out = ','.join(Output)
#                
#        Cond = ' %s AND '.join(Conditions.keys())
#        Cond = Cond + '%s'
#        query = """ SELECT {} 
#            FROM Trts 
#            INNER JOIN Devices ON Devices.idDevices = Trts.Device_id
#            INNER JOIN Wafers ON Wafers.idWafers = Devices.Wafer_id		  	
#            INNER JOIN TrtTypes ON TrtTypes.idTrtTypes = Trts.Type_id
#            WHERE {} """.format(Out,Cond)   
#       
#        self.cur.execute(query,Conditions.values())        
#        res = self.cur.fetchall()        
#        ret = []        
#        for ir,r in enumerate(res):
#            do = {}
#            for io,of in enumerate(Output): do[of]=r[io]
#            ret.insert(len(ret),do)
#            
#        return ret
#    

    
#    def GetCharactFromID (self, Table, Ids):
#        Out = '{}.id{}, {}.Data, {}.MeasDate, Trts.Name'.format(Table,Table,Table,Table,Table)
#        Cond = ' OR '.join(['{}.id{}=%s'.format(Table,Table)]*len(Ids))
#        query = '''SELECT {} 
#                FROM {} 
#                INNER JOIN Trts ON Trts.idTrts={}.Trt_id
#                WHERE {} ORDER BY {}.MeasDate'''.format(Out,Table,Table,Cond,Table)
#        
#        self.cur.execute(query,Ids)
#        res = self.cur.fetchall()        
#        
#        Data = {}        
#        for r in res:
##            print 'Size', len(r[1])
#            Data['{}{}'.format(r[3],r[0])] = pickle.loads(r[1])               
#        
#        return Data
#        

        
        
