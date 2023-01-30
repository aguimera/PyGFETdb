import pickle

import numpy as np
from PyQt5 import Qt, QtCore
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from pyqtgraph.parametertree import ParameterTree, Parameter
import pyqtgraph.parametertree.parameterTypes as pTypes
from PyGFETdb import DataClass as Datacl
import quantities as pq

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


###############################################################################
# Generic voltage sweep
###############################################################################

VoltageSweepParams = ({'name': 'Start',
                       'type': 'float',
                       'value': -0.2,
                       'default': -0.2,
                       'step': 0.01,
                       'siPrefix': True,
                       'suffix': 'V'},
                      {'name': 'Stop',
                       'type': 'float',
                       'value': 0.6,
                       'default': 0.6,
                       'step': 0.01,
                       'siPrefix': True,
                       'suffix': 'V'},
                      {'name': 'bPoints',
                       'title': 'Fixed points',
                       'type': 'bool',
                       'value': True,
                       'default': True},
                      {'name': 'Points',
                       'type': 'int',
                       'value': 100,
                       'step': 1,
                       'default': 100,
                       'suffix': 'Points',
                       'readonly': False},
                      {'name': 'bStep',
                       'title': 'Fixed step',
                       'type': 'bool',
                       'value': False,
                       'default': False},
                      {'name': 'Step',
                       'type': 'float',
                       'value': 0.01,
                       'siPrefix': True,
                       'step': 0.005,
                       'suffix': 'V',
                       'readonly': True},
                      {'name': 'Sweep',
                       'title': 'Sweep Points',
                       'type': 'text',
                       'expanded': False,
                       'readonly': True}
                      )


class VoltageSweepConfig(pTypes.GroupParameter):
    def __init__(self, **kwargs):
        pTypes.GroupParameter.__init__(self, **kwargs)
        self.addChildren(VoltageSweepParams)

        self.param('bPoints').sigValueChanged.connect(self.on_bPoints)
        self.param('bStep').sigValueChanged.connect(self.on_bStep)

        self.param('Start').sigValueChanged.connect(self.CalcSweep)
        self.param('Stop').sigValueChanged.connect(self.CalcSweep)
        self.param('Points').sigValueChanged.connect(self.CalcSweep)
        self.param('Step').sigValueChanged.connect(self.CalcSweep)

        self.CalcSweep()

    def on_bPoints(self, value):
        self.param('bStep').setValue(not value.value(),
                                     blockSignal=self.on_bStep)
        self.param('Step').setOpts(readonly=value.value())
        self.param('Points').setOpts(readonly=not value.value())
        self.CalcSweep()

    def on_bStep(self, value):
        self.param('bPoints').setValue(not value.value(),
                                       blockSignal=self.on_bPoints)
        self.param('Step').setOpts(readonly=not value.value())
        self.param('Points').setOpts(readonly=value.value())
        self.CalcSweep()

    def CalcSweep(self):
        if self.param('bPoints').value():
            sweep = np.linspace(self.param('Start').value(),
                                self.param('Stop').value(),
                                self.param('Points').value())
            if len(sweep) > 1:
                step = np.mean(np.diff(sweep))
            else:
                step = 0
            self.param('Step').setValue(step,
                                        blockSignal=self.CalcSweep)

        else:
            sweep = np.arange(self.param('Start').value(),
                              self.param('Stop').value(),
                              self.param('Step').value())
            self.param('Points').setValue(len(sweep),
                                          blockSignal=self.CalcSweep)

        self.SweepVals = sweep
        self.param('Sweep').setValue(str(sweep))


ElectParamParams = ({'name': 'Param',
                     'title': 'Parameter',
                     'type': 'list',
                     'values': Datacl.ParametersList,
                     'value': 'Ids',
                     'default': 'Ids',
                     },
                    {'name': 'Vgs',
                     'type': 'list',
                     'values': ('Vgs', 'VgsNorm', 'SinglePoint'),
                     'value': 'Vgs',
                     'default': 'Vgs',
                     'expanded': False,
                     'children': [{'name': 'Vgs',
                                   'type': 'float',
                                   'value': 0.1,
                                   'default': 0.1,
                                   'siPrefix': True,
                                   'suffix': 'V',
                                   'visible': False,
                                   }]
                     },
                    )

ExtraEPars = {'NFmin': {'name': 'NFMin',
                        'type': 'float',
                        'value': 5,
                        'default': 5,
                        'siPrefix': True,
                        'suffix': 'Hz',
                        'removable': True,
                        },
              'NFmax': {'name': 'NFmax',
                        'type': 'float',
                        'value': 5e3,
                        'default': 5e3,
                        'siPrefix': True,
                        'suffix': 'Hz',
                        'removable': True,
                        },
              'FFmin': {'name': 'FFmin',
                        'type': 'float',
                        'value': 5,
                        'default': 5,
                        'siPrefix': True,
                        'suffix': 'Hz',
                        'removable': True,
                        },
              'FFmax': {'name': 'FFmax',
                        'type': 'float',
                        'value': 7e3,
                        'default': 7e3,
                        'siPrefix': True,
                        'suffix': 'Hz',
                        'removable': True,
                        },
              'Units': {'name': 'Units',
                        'type': 'str',
                        'value': '',
                        'removable': True,
                        }
              }


class ElecParam(pTypes.GroupParameter):
    def __init__(self, VgsNormVals, VgsVals, cldatkw, **opts):
        self.VgsNormVals = VgsNormVals
        self.VgsVals = VgsVals
        opts['type'] = 'group'
        opts['addText'] = "Add arg"
        opts['addList'] = list(ExtraEPars.keys())
        pTypes.GroupParameter.__init__(self, **opts)
        self.addChildren(ElectParamParams)
        self.param('Param').setValue(cldatkw['Param'])
        self.param('Param').setDefault(cldatkw['Param'])

        self.param('Vgs').sigValueChanged.connect(self.VgsChanged)

        if cldatkw['Param'] in ('Ud0', 'IgMax'):
            self.param('Vgs').setOpts(**{'visible': False})
            self.param('Vgs').setValue('SinglePoint')
        else:
            self.addChild({'name': 'Ud0Norm',
                           'type': 'bool'})
            if 'Ud0Norm' in cldatkw:
                self.param('Ud0Norm').setValue(cldatkw['Ud0Norm'])
                self.param('Ud0Norm').setDefault(cldatkw['Ud0Norm'])
            else:
                self.param('Ud0Norm').setValue(False)
                self.param('Ud0Norm').setDefault(False)

            if cldatkw['Vgs'].size == 1:
                self.param('Vgs').setValue('SinglePoint')
                self.param('Vgs').setDefault('SinglePoint')
                self.param('Vgs').param('Vgs').setValue(cldatkw['Vgs'])
                self.param('Vgs').param('Vgs').setDefault(cldatkw['Vgs'])
            else:
                if self.param('Ud0Norm').value():
                    self.param('Vgs').setValue('VgsNorm')
                    self.param('Vgs').setDefault('VgsNorm')
                else:
                    self.param('Vgs').setValue('Vgs')
                    self.param('Vgs').setDefault('Vgs')

        for k, d in cldatkw.items():
            if k in ('Param', 'Ud0Norm', 'Vgs'):
                continue
            self.addChild({'name': k,
                           'type': 'str',
                           'value': d,
                           })

    def addNew(self, Arg):
        if Arg in [p.name() for p in self.children()]:
            return
        self.addChild(ExtraEPars[Arg])

    def VgsChanged(self):
        if self.param('Vgs').value() == 'SinglePoint':
            self.param('Vgs').param('Vgs').setOpts(**{'visible': True})
        else:
            self.param('Vgs').param('Vgs').setOpts(**{'visible': False})

    def GetArgKwargs(self):
        kwargs = {}
        for p in self.children():
            if not p.opts['visible']:
                continue
            if p.name() == 'Vgs':
                if p.value() == 'SinglePoint':
                    kwargs[p.name()] = p.param('Vgs').value() * pq.V
                elif p.value() == 'VgsNorm':
                    kwargs[p.name()] = self.VgsNormVals.SweepVals * pq.V
                elif p.value() == 'Vgs':
                    kwargs[p.name()] = self.VgsVals.SweepVals * pq.V
            else:
                kwargs[p.name()] = p.value()
        return kwargs


class ElecParamsEditor(pTypes.GroupParameter):
    def __init__(self, ClassQueries, pdAttr, **kwargs):
        pTypes.GroupParameter.__init__(self, **kwargs)
        self.addChild(VoltageSweepConfig(name='VgsVals',
                                         title='Vgs Points',
                                         expanded=False))
        self.addChild(VoltageSweepConfig(name='VgsNormVals',
                                         title='Vgs Normalized Points',
                                         expanded=False))

        self.param('VgsNormVals').param('Sweep').setOpts(**{'title': 'VgsNorm'})
        self.param('VgsNormVals').param('Start').setOpts(**{'value': np.min(pdAttr['VgsNorm']),
                                                            'default': np.min(pdAttr['VgsNorm'])})
        self.param('VgsNormVals').param('Stop').setOpts(**{'value': np.max(pdAttr['VgsNorm']),
                                                           'default': np.max(pdAttr['VgsNorm'])})
        self.param('VgsNormVals').param('Points').setOpts(**{'value': pdAttr['VgsNorm'].size,
                                                             'default': pdAttr['VgsNorm'].size})

        self.param('VgsVals').param('Sweep').setOpts(**{'title': 'Vgs'})
        self.param('VgsVals').param('Start').setOpts(**{'value': np.min(pdAttr['Vgs']),
                                                        'default': np.min(pdAttr['Vgs'])})
        self.param('VgsVals').param('Stop').setOpts(**{'value': np.max(pdAttr['Vgs']),
                                                       'default': np.max(pdAttr['Vgs'])})
        self.param('VgsVals').param('Points').setOpts(**{'value': pdAttr['Vgs'].size,
                                                         'default': pdAttr['Vgs'].size})

        EPars = []
        for pname, cldatkw in ClassQueries.items():
            EPars.append(ElecParam(name=pname,
                                   renamable=True,
                                   removable=True,
                                   expanded=False,
                                   VgsVals=self.param('VgsVals'),
                                   VgsNormVals=self.param('VgsNormVals'),
                                   cldatkw=cldatkw))

        self.addChild({'name': 'Queries',
                       'title': 'Electrical Values',
                       'type': 'group',
                       'addText': 'Add Calc',
                       'addList': ['New1', ''],
                       'children': EPars})

        self.param('Queries').sigAddNew.connect(self.on_addNew)

    def on_addNew(self, Arg):
        self.param('Queries').addChild(ElecParam(name='New',
                                                 renamable=True,
                                                 removable=True,
                                                 expanded=False,
                                                 VgsVals=self.param('VgsVals'),
                                                 VgsNormVals=self.param('VgsNormVals'),
                                                 cldatkw={'Param': 'Ids',
                                                          'Vgs': 0 * pq.V}))

    def GetQueries(self):
        Quer = {}
        ArrayCols = []
        ScalarCols = []
        ArrayColsNorm = []
        for p in self.param('Queries').children():
            Quer[p.name()] = p.GetArgKwargs()
            if p.param('Vgs').value() == 'SinglePoint':
                ScalarCols.append(p.name())
            else:
                ArrayCols.append(p.name())
                if 'Ud0Norm' in Quer[p.name()]:
                    if Quer[p.name()]['Ud0Norm']:
                        ArrayColsNorm.append(p.name())

        pdAttr = {'Vgs': self.param('VgsVals').SweepVals * pq.V,
                  'VgsNorm': self.param('VgsNormVals').SweepVals * pq.V,
                  'ScalarCols': ScalarCols,
                  'ArrayCols': ArrayCols,
                  'ArrayColsNorm': ArrayColsNorm,
                  }

        return Quer, pdAttr


class ElectricalParamsEditor(Qt.QDialog):
    ''' Main Window '''

    def __init__(self, ClassQueries, pdAttr):
        super(ElectricalParamsEditor, self).__init__()
        self.setWindowModality(QtCore.Qt.ApplicationModal)

        layout = Qt.QVBoxLayout(self)
        self.ClassQueries = ClassQueries
        self.pdAttr = pdAttr

        self.setGeometry(650, 20, 400, 800)
        self.setWindowTitle('Edit Electrical Parameters')

        self.ElecParsConf = ElecParamsEditor(ClassQueries=ClassQueries,
                                             pdAttr=pdAttr,
                                             QTparent=self,
                                             name='Elect')
        self.Parameters = Parameter.create(name='App Parameters',
                                           type='group',
                                           children=({'name': 'SaveConf',
                                                      'type': 'action'},
                                                     {'name': 'LoadConf',
                                                      'type': 'action'},
                                                     self.ElecParsConf,))
        self.Parameters.param('SaveConf').sigActivated.connect(self.on_save)
        self.Parameters.param('LoadConf').sigActivated.connect(self.on_load)

        self.treepar = ParameterTree()
        self.treepar.setParameters(self.Parameters, showTop=False)

        layout.addWidget(self.treepar)

    def on_save(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self, "Save PKl file", "", "All Files (*);; (*.pkl)", options=options)
        if fileName:
            if not fileName.endswith('.pkl'):
                fileName = fileName + '.pkl'
            ClassQueries, pdAttr = self.ElecParsConf.GetQueries()
            pickle.dump({'ClassQueries': ClassQueries,
                         'pdAttr': pdAttr,
                         }, open(fileName, 'wb'))

    def on_load(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "Load PKl file", "", "All Files (*);; (*.pkl)", options=options)
        if fileName:
            try:
                d = pickle.load(open(fileName, 'rb'))
                self.ClassQueries = d['ClassQueries']
                self.pdAttr = d['pdAttr']
                self.Parameters.removeChild(self.Parameters.param('Elect'))
                self.Parameters.addChild(ElecParamsEditor(ClassQueries=self.ClassQueries,
                                                          pdAttr=self.pdAttr,
                                                          QTparent=self,
                                                          name='Elect'))
            except:
                print('erro loading file')
                return

    def closeEvent(self, event):
        self.ClassQueries, self.pdAttr = self.ElecParsConf.GetQueries()
