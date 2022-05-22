from PyQt5 import Qt
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


class TableEditor(pTypes.GroupParameter):
    def __init__(self, TableColumns, DB, Table, **kwargs):
        super(TableEditor, self).__init__(**kwargs)
        self.DB = DB
        self.Table = Table
        self.addChild({'name': 'Update',
                       'title': 'Update Record',
                       'type': 'action'}, )
        for n, col in TableColumns.items():
            self.addChild(col)

        self.param('Update').sigActivated.connect(self.on_update)

    def on_update(self):
        Values = {}
        for p in self.children():
            if p.type() == 'action':
                continue
            if p.name().startswith('id'):
                scond = '{}={}'.format(p.name(), p.value())
                continue
            Values[p.name()] = p.value()

        print('Overwrite Record id ', scond)

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
            self.Records.append(TableEditor(TableColumns=rec,
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
        print('allupdate')