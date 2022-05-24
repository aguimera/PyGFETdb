
from PyQt5 import Qt
from qtpy.QtWidgets import QMessageBox
from pyqtgraph.parametertree import ParameterTree, Parameter
import pyqtgraph.parametertree.parameterTypes as pTypes

EditWaferColumns = {'idWafers': {'name': 'idWafers',
                                 'type': 'int',
                                 'readonly': True
                                 },
                    'Name': {'name': 'Name',
                             'type': 'str',
                             'readonly': True
                             },
                    'Run': {'name': 'Run',
                            'type': 'str',
                            },
                    'Substrate': {'name': 'Substrate',
                                  'type': 'str',
                                  },
                    'Masks': {'name': 'Masks',
                              'type': 'str',
                              },
                    'Graphene': {'name': 'Graphene',
                                 'type': 'str',
                                 },
                    'Comments': {'name': 'Comments',
                                 'type': 'str',
                                 },
                    }

EditDeviceColumns = {'idDevices': {'name': 'idDevices',
                                   'type': 'int',
                                   'readonly': True
                                   },
                     'Name': {'name': 'Name',
                              'type': 'str',
                              'readonly': True
                              },
                     'Wafer_id': {'name': 'Wafer_id',
                                  'type': 'int',
                                  },
                     'State': {'name': 'State',
                               'type': 'str',
                               },
                     'ExpOK': {'name': 'ExpOK',
                               'type': 'int',
                               },
                     'Comments': {'name': 'Comments',
                                  'type': 'str',
                                  },
                     }

EditTrtColumns = {'idTrts': {'name': 'idTrts',
                             'type': 'int',
                             'readonly': True
                             },
                  'Name': {'name': 'Name',
                           'type': 'str',
                           'readonly': True
                           },
                  'Device_id': {'name': 'Device_id',
                                'type': 'int',
                                },
                  'Type_id': {'name': 'Type_id',
                              'type': 'int',
                              },
                  'Comments': {'name': 'Comments',
                               'type': 'str',
                               },
                  }

EditTrtTypesColumns = {'idTrtTypes': {'name': 'idTrtTypes',
                                      'type': 'int',
                                      'readonly': True
                                      },
                       'Name': {'name': 'Name',
                                'type': 'str',
                                'readonly': True
                                },
                       'Length': {'name': 'Length',
                                  'type': 'float',
                                  'siPrefix': True,
                                  'suffix': 'm',
                                  },
                       'Width': {'name': 'Width',
                                 'type': 'float',
                                 'siPrefix': True,
                                 'suffix': 'm',
                                 },
                       'Pass': {'name': 'Pass',
                                'title': 'Passivation',
                                'type': 'float',
                                'siPrefix': True,
                                'suffix': 'm',
                                },
                       'Area': {'name': 'Area',
                                'type': 'float',
                                'siPrefix': True,
                                'suffix': 'm2',
                                },
                       'ContactOverlap': {'name': 'ContactOverlap',
                                          'type': 'float',
                                          'siPrefix': True,
                                          'suffix': 'm',
                                          },
                       'Contact': {'name': 'Contact',
                                   'type': 'str',
                                   },
                       'Shape': {'name': 'Shape',
                                 'type': 'str',
                                 },
                       'Comments': {'name': 'Comments',
                                    'type': 'str',
                                    },
                       }

EditDCCharColumns = {'idDCcharacts': {'name': 'idDCcharacts',
                                      'type': 'int',
                                      'readonly': True
                                      },
                     'UpdateDate': {'name': 'UpdateDate',
                                    'type': 'str',
                                    'readonly': True,
                                    },
                     'MeasDate': {'name': 'MeasDate',
                                  'type': 'str',
                                  'readonly': True,
                                  },

                     'User_id': {'name': 'User_id',
                                 'type': 'int',
                                 'readonly': True
                                 },
                     'Trt_id': {'name': 'Trt_id',
                                'type': 'int',
                                'readonly': True
                                },
                     'Gate_id': {'name': 'Gate_id',
                                 'type': 'int',
                                 'readonly': True
                                 },
                     'IsOk': {'name': 'IsOk',
                              'type': 'int',
                              },
                     'IsCmp': {'name': 'IsCmp',
                               'type': 'int',
                               },
                     'Ph': {'name': 'Ph',
                            'type': 'float',
                            },
                     'IonStrength': {'name': 'IonStrength',
                                     'type': 'float',
                                     },
                     'FuncStep': {'name': 'FuncStep',
                                   'type': 'str',
                                   },
                     'AnalyteCon': {'name': 'AnalyteCon',
                                    'type': 'float',
                                    },
                     'Comments': {'name': 'Comments',
                                  'type': 'str',
                                  },
                     'FileName': {'name': 'FileName',
                                  'type': 'str',
                                  'readonly': True
                                  },
                     }

EditACCharColumns = {'idACcharacts': {'name': 'idACcharacts',
                                      'type': 'int',
                                      'readonly': True
                                      },
                     'UpdateDate': {'name': 'UpdateDate',
                                    'type': 'str',
                                    'readonly': True,
                                    },
                     'MeasDate': {'name': 'MeasDate',
                                  'type': 'str',
                                  'readonly': True,
                                  },

                     'User_id': {'name': 'User_id',
                                 'type': 'int',
                                 'readonly': True
                                 },
                     'Trt_id': {'name': 'Trt_id',
                                'type': 'int',
                                'readonly': True
                                },
                     'DC_id': {'name': 'DC_id',
                               'type': 'int',
                               'readonly': True
                                 },
                     'IsOk': {'name': 'IsOk',
                              'type': 'int',
                              },
                     'IsCmp': {'name': 'IsCmp',
                               'type': 'int',
                               },
                     'Ph': {'name': 'Ph',
                            'type': 'float',
                            },
                     'IonStrength': {'name': 'IonStrength',
                                     'type': 'float',
                                     },
                     'FuncStep': {'name': 'FuncStep',
                                   'type': 'str',
                                   },
                     'AnalyteCon': {'name': 'AnalyteCon',
                                    'type': 'float',
                                    },
                     'Comments': {'name': 'Comments',
                                  'type': 'str',
                                  },
                     'FileName': {'name': 'FileName',
                                  'type': 'str',
                                  'readonly': True
                                  },
                     }

class TableEditor(pTypes.GroupParameter):
    def __init__(self, QtParent, TableColumns, DB, Table, **kwargs):
        super(TableEditor, self).__init__(**kwargs)
        self.QtParent = QtParent
        self.DB = DB
        self.Table = Table
        self.addChild({'name': 'Update',
                       'title': 'Update Record',
                       'type': 'action'}, )
        for n, col in TableColumns.items():
           self.addChild(col)

        self.param('Update').sigActivated.connect(self.on_update)

    def SetValue(self, pName, pValue):
        self.param(pName).setValue(pValue)

    def on_update(self):
        Values = {}
        sMsg = ''
        for p in self.children():
            if p.type() == 'action':
                continue
            if p.name().startswith('id'):
                scond = '{}={}'.format(p.name(), p.value())
                continue
            if p.valueIsDefault():
                continue
            Values[p.name()] = p.value()
            sMsg = "{} \n {} -> {}".format(sMsg, p.name(), p.value())

        Msg = "Do you want to update following Record {} {}".format(scond, sMsg)
        buttonReply = QMessageBox.question(self.QtParent,
                                           'WARNING !!!!!!',
                                           Msg,
                                           QMessageBox.Yes | QMessageBox.No,
                                           QMessageBox.No)

        if buttonReply == QMessageBox.Yes:
            res = self.DB.UpdateRow(Table=self.Table,
                                    Fields=Values,
                                    Condition=scond)
            print(res)


class TableEditorWindow(Qt.QWidget):
    ''' Main Window '''

    def __init__(self, Table, Records, DB):
        super(TableEditorWindow, self).__init__()
        layout = Qt.QVBoxLayout(self)

        self.setGeometry(650, 20, 400, 800)
        self.setWindowTitle('Edit {}'.format(Table, ))

        self.Records = []
        for ic, rec in enumerate(Records):
            self.Records.append(TableEditor(QtParent=self,
                                            TableColumns=rec,
                                            DB=DB,
                                            Table=Table,
                                            name='Rec{}'.format(ic)))

        if len(Records) > 1:
            self.Common = Parameter.create(name='CommonRecord',
                                           title='Common Values',
                                           type='group')
            self.Common.addChild({'name': 'Update',
                                  'title': 'Update ALL Records',
                                  'type': 'action'})

            for k, v in Records[0].items():
                if k.startswith('id'):
                    v['value'] = 0
                if k == 'Name':
                    v['value'] = ''
                self.Common.addChild(v)
            self.Records.insert(0, self.Common)
            self.Common.param('Update').sigActivated.connect(self.on_update)

        self.Parameters = Parameter.create(name='App Parameters',
                                           type='group',
                                           children=self.Records)
        self.treepar = ParameterTree()
        self.treepar.setParameters(self.Parameters, showTop=False)
        layout.addWidget(self.treepar)

    def on_update(self):
        update = False
        for p in self.Common.children():
            if p.type() == 'action':
                continue
            if p.name().startswith('id'):
                continue
            if p.valueIsDefault():
                continue
            for r in self.Records[1:]:
                update = True
                r.SetValue(p.name(), p.value())

        if update:
            for r in self.Records[1:]:
                r.on_update()
