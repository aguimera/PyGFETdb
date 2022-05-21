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
    def __init__(self, TableColumns, DB, **kwargs):
        super(TableEditor, self).__init__(**kwargs)
        self.DB = DB
        self.addChild({'name': 'Update',
                       'title': 'UpdateValues',
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

        res = self.DB.UpdateRow(Table=self.opts['name'],
                                Fields=Values,
                                Condition=scond)

        print(res)


class TableEditorWindow(Qt.QWidget):
    ''' Main Window '''

    def __init__(self, Table, TableColumns, DB):
        super(TableEditorWindow, self).__init__()
        layout = Qt.QVBoxLayout(self)

        self.setGeometry(650, 20, 400, 800)
        self.setWindowTitle('Edit {}'.format(Table, ))

        self.ViewConfig = TableEditor(TableColumns, DB,
                                      name=Table)

        self.treepar = ParameterTree()
        self.treepar.setParameters(self.ViewConfig, showTop=False)
        layout.addWidget(self.treepar)
