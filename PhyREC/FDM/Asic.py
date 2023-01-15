import pickle
import os
from PhyREC.FDM.DicAddress import DicAddress
import numpy as np
import PhyREC

# fileDir = os.path.dirname(os.path.realpath('__file__'))

fileDir = os.path.dirname(PhyREC.__file__)

class ASIC:

    def __init__(self,
                 A_File='S0_AS3_T_220101_V1',
                 P_File='AS3_F3_220128_C4',
                 Fs = 51.2e6/32/2.0,
                 NumOfCols = 26,
                 DEBUG = False
                 ):

        # print(DicAddress.FDM_Row)
        self.Fs = Fs
        self.NumOfCols = NumOfCols

        self.DEBUG = DEBUG
        Config_Dir = os.path.join(fileDir, 'FDM/Config/')

        self.A_File = os.path.join(Config_Dir, A_File)

        self.P_File = os.path.join(Config_Dir, P_File)

        # ##Sample CAL DICs
        self.Calibration_FDM()

    def Calibration_FDM(self):
        if self.DEBUG:
            print("CALIBRATION SAMPLE ", self.A_File)
        with open(self.A_File, "rb") as f:
            DcDat = pickle.load(f, encoding='latin1')

        self.Container = DcDat
        self.R_Gain = self.Container['RT']['Container_Res'][0, :]
        self.Gain = {}

        self.Gain[0] = self.Container['RT']['Container_Res'][1, :]
        self.Gain[1] = self.Container['RT']['Container_Res'][1, :]
        self.Gain[2] = self.Container['RT']['Container_Res'][2, :]
        self.Gain[4] = self.Container['RT']['Container_Res'][3, :]
        self.Gain[8] = self.Container['RT']['Container_Res'][4, :]

        self.Amp_I_C = {}
        self.Amp_I_C[10] = self.Container['TI']['Vp'][0, :]
        self.Amp_I_C[20] = self.Container['TI']['Vp'][1, :]
        self.Amp_I_C[50] = self.Container['TI']['Vp'][2, :]
        self.Amp_I_C[100] = self.Container['TI']['Vp'][3, :]

        self.Amp_Q_C = {}
        self.Amp_Q_C[10] = self.Container['TQ']['Vp'][0, :]
        self.Amp_Q_C[20] = self.Container['TQ']['Vp'][1, :]
        self.Amp_Q_C[50] = self.Container['TQ']['Vp'][2, :]
        self.Amp_Q_C[100] = self.Container['TQ']['Vp'][3, :]

        self.Freq = self.Container['TI']['Freq'][0, :]
        if 'FDem' in self.Container['TI'].keys():
            self.Freq_Dem = self.Container['TI']['FDem']
            self.Steps = self.Container['TI']['Steps']
        else:
            self.Freq_Dem = None
            self.Steps = None

        ######Phase Dic##########

        if self.P_File is not None:
            if self.DEBUG:
                print("Phase ", self.P_File)
            with open(self.P_File, "rb") as f:
                self.DcDat_F = pickle.load(f, encoding='latin1')

            self.Phase_Q = self.DcDat_F['Q']
            self.Phase_I = self.DcDat_F['I']
        else:
            self.Phase_Q = np.zeros(16)
            self.Phase_I = np.zeros(16)

    def RowAndByteArrayToIamp(self,Buffer_Fetch,Row,Gain,Ideal=None, VFS=2.0):

        Factor = 1
        R_Gain = 5e3

        Buffer_T = VFS * (Buffer_Fetch - (2 ** (13.0 - 1)) + 0.5) / 2 ** 13.0

        if Ideal == True:
            print("Conversion Ideal")
            Factor = Factor * 2.0 * Gain
        else:

            Factor = 2.0 * self.Gain[Gain][Row]
            R_Gain = self.R_Gain[Row]

        Buffer_T = (Buffer_T / Factor) / R_Gain

        return Buffer_T


    def _FDM_ByteArrayToIamp(self, Buffer_Fetch, Gain = 1, DIC_R=None, Ideal=None, VFS=2.0, Index = 1):

        if self.Checker(Buffer_Fetch):
            print("Error: Column periodicity or header not found")
        # Columns = Buffer_Fetch[:,0] & 0x001f
        TTL = (Buffer_Fetch[:,0] & 0x00e0) >> 5

        Buffer_T = VFS * ((Buffer_Fetch[:,Index:]) - (2 ** (13.0 - 1)) + 0.5) / 2 ** 13.0

        ##Factor 2*Cs_CF RGain
        Factor = np.ones(self.NumOfCols-1)
        R_Gain = np.ones(self.NumOfCols-1) * 5e3

        if Ideal == True:
            print("Conversion Ideal")
            Factor = Factor * 2.0 * Gain
            # print("Factor ",Factor ," R_Gain ",R_Gain)
        else:
            if DIC_R is not None:
                print("DIC_R is not None")
                for n in DIC_R:
                    if n[:3] == 'Row':
                        if DIC_R[n]['Enable'] == True:
                            print(n)
                            Pointer = int(n[4:])
                            Val = DIC_R[n]['Gain']

                            Pointer_Buffer = DicAddress.FDM_Row[Pointer] ##Needed to locate the correct Row

                            Factor[Pointer_Buffer] = 2.0 * self.Gain[Val][Pointer]
                            R_Gain[Pointer_Buffer] = self.R_Gain[Pointer]
                            # print("Factor ",Factor[Pointer] ," R_Gain ",R_Gain[Pointer])
            else:
                print("DIC_R is None")
                for n in range(32):
                    Factor[n] = 2.0 * self.Gain[Gain][n]
                    R_Gain[n] = self.R_Gain[n]

                    ## Iamp Current
        Buffer_T = (Buffer_T / Factor) / R_Gain

        return Buffer_T,TTL

    def Checker(self,Buffer):

        error = False
        if all(Buffer[:, 0] >= 40960):
            Buffer_Column = Buffer[:, 0] & 0x001f  # 5'b1 Column counter
            Pointers = np.where(np.diff(np.where(np.diff(Buffer_Column) != 1)[0]) != 32)[0]

            if len(Pointers) > 0:
                print("ERROR on received Columns")
                error = True
        else:
            print("Error on received Header")
            error = True

        return error


    def _FDM_DIC_RtoInclude(self, DIC_R):

        Row_V = np.zeros(32)
        Col_V = np.zeros(32)
        Gain_V = np.zeros(32)

        for n in DIC_R.keys():
            if n[:3] == 'Col':
                if DIC_R[n]['OE']:
                    Pointer = int(n[3:])
                    if DIC_R[n]['En_Q']:
                        Col_V[Pointer] = 1
                    # if DIC_R[n]['En_I']:
                    #     Col_V[Pointer + 16] = 1

            if n[:3] == 'Row':
                if DIC_R[n]['Enable']:
                    Pointer = int(n[3:])
                    Row_V[Pointer] = 1

                    Gain_V[Pointer] = DIC_R[n]['Gain']



        Col_Tot = np.sum(Col_V)
        Row_Tot = np.sum(Row_V)
        Vector_Length = int(Col_Tot * Row_Tot)

        Include = {}
        ColRow = {}
        ColRow_D = {}
        Cont_App = []
        Phase_Cont = []
        Rows = []
        Cols = []
        Contador = 0
        NameRowCol = []
        for i, Row in enumerate(Row_V):
            for j, Col in enumerate(Col_V):
                if Col and Row:
                    if j > 15:
                        y = int(j - 16)
                        Vp = self.Amp_I_C[DIC_R['Col ' + str(y)]['Amp_I']][y]
                        Phase_Cont.append(self.Phase_I[y])
                    else:
                        y = j
                        Vp = self.Amp_Q_C[DIC_R['Col ' + str(j)]['Amp_Q']][j]
                        Phase_Cont.append(self.Phase_Q[y])
                    # print(y)

                    Freq = self.Freq[y]
                    if self.Freq_Dem is not None:
                        Freq_D = self.Freq_Dem[y]
                    else:
                        Freq_D = Freq

                    ColRow['R'+str(i)+'C'+str(j)] = [int(i) , int(j), int(DicAddress.FDM_Row[i]), int(Contador), Vp, Freq]
                    ColRow_D['R'+str(i)+'C'+str(j)] = [int(i) ,int(j), int(DicAddress.FDM_Row[i]), int(Contador), Vp, Freq_D,Gain_V[i]]
                    Cont_App.append('R'+str(i)+'C'+str(j))
                    NameRowCol.append('R'+str(i)+'C'+str(j))
                    Rows.append(int(DicAddress.FDM_Row[i]))
                    Cols.append(int(j))

                    Contador = Contador + 1
        Include["DicMux"] = ColRow
        Include["DicMuxDem"] = ColRow_D
        Include['Include'] = Cont_App
        Include['NameRowCol'] = NameRowCol
        Include['Phase'] = Phase_Cont
        Include['CT'] = int(Col_Tot)
        Include['RT'] = int(Row_Tot)
        Include['VL'] = Vector_Length
        Include['Rows'] = np.unique(np.array(Rows))
        Include['Cols'] = np.unique(np.array(Cols))
        Include['Info'] = ['Col', 'Row', 'DicAddress', 'VpeakAmpQ', 'Freq']

        return Include


