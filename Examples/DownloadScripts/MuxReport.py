import pandas as pd
import matplotlib as mpl
from PyGFETdb.GuiDBView.GuiHelpers import GenScalarFigures, GenVectorFigures, AddColRow, GenDeviceReport, AddCycle, \
    GenDeviceReportGui
import seaborn as sns

mpl.use("Qt5Agg")

dfData = pd.read_pickle('DCData.pkl')

# %% Add columns

AddColRow(dfData)

# %% Calc cycles

AddCycle(dfData)

# %% GenCycle figure

# GenDeviceReport(Data=dfData,
#                 FileOut='ReportMux.pdf',
#                 ScalarPars=('CNP', 'GM01', 'Ids01', 'RdsCNP'),
#                 ScalarPltFunct=sns.swarmplot,
#                 VectorPars=('Ids', 'GM'),
#                 )

GenDeviceReportGui(Data=dfData,
                   FileOut='ReportMux.pdf',
                   ScalarPars=('CNP', 'GM01', 'Ids01', 'RdsCNP'),
                   ScalarPltFunct=sns.swarmplot,
                   VectorPars=('Ids', 'GM'),
                   )
