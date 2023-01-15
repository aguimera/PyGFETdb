import numpy as np
from scipy import signal
import pickle

class Processing():

    def __init__(self, Ovf_Lim = 0.98):
        print("Process")
        self.Ovf_Lim = Ovf_Lim

    def OVF(self,Data):
        return  np.max(np.abs(Data)) >= self.Ovf_Lim

    def PSD_W(self, Data, Fs, Nseg=0.125 * 0.25 * 0.25, scaling='density'):
        FC1_Fxx, FC1_Pxx = signal.welch(Data, fs=Fs, nperseg=int(Nseg * len(Data)), scaling=scaling)

        return FC1_Fxx, FC1_Pxx

    def RemoveDC(self, sig, Type='constant'):
        st = np.array(sig)
        st = signal.detrend(st, type=Type, axis=0)
        return st  # sig.duplicate_with_new_data(signal=st*sig.units)

    def GetDC(self, sig, Type='constant'):
        st = np.array(sig)
        st = st - signal.detrend(st, type=Type, axis=0)
        st = np.mean(st)
        return st  # sig.duplicate_with_new_data(signal=st*sig.units)

    def DownSampling(self, sig, Fact, zero_phase=True):
        print(Fact)
        rs = signal.decimate(np.array(sig),
                             q=Fact,
                             zero_phase=zero_phase,
                             ftype='iir',
                             axis=0)
        return rs

    def NewTime(self, sig, Fs):
        return np.arange(0.0, len(sig) / Fs, 1.0 / Fs)


    def Filter(self, sig, Type, Order, Freqs, Fs):
        st = np.array(sig)

        # print("Fs", Fs)
        # Fs = sig.sampling_rate
        Freqs = np.array(Freqs)
        freqs = Freqs / (0.5 * Fs)

        sos = signal.butter(Order, freqs, Type, output='sos')
        st = signal.sosfiltfilt(sos, st, axis=0)

        return st

    def PeakAmp(self, sig, Tone=10, Leak=10):

        # sig = row['Psw_F']
        # Leak = 10
        # Tone = 10

        Bin = sig[0][1]

        Freq_Index = int(Tone / Bin)

        Max_peak = np.max(sig[1][Freq_Index - Leak:Freq_Index + Leak])
        try:
            print("Peak Hz ", sig[0][np.where(sig[1] == Max_peak)[0]][0])
        except:
            print("PEAK HZ Not Found", Max_peak)

        return np.sqrt(Max_peak) * 2 ** 0.5

    def GmFromPeak(self, sig, Amp):

        return sig / Amp

    def Gain(self, sig):
        return sig[0] * sig[1]

    def GainVin(self, sig):
        print("Gain", sig[1])
        # print("sig",sig)
        try:
            A = sig[0] / sig[1]
        except:
            A = 0
        return A


    def sum_One(self, sig):
        return sig + 1


    def Dem(self, sig, index, Include_C, Filt_Freq=2000.0, Order=10, Fs=800e3, Time_R=0.0):
        print("Index ", index, " shape ", np.shape(sig))

        Ts_M = 1.0 / Fs
        Time = np.arange(0.0, len(sig) * Ts_M, Ts_M)

        Col, Row, Pointer, Vds, Freq = Include_C['DicMuxDem'][index]

        Sin_C = 2.0 * np.sin(2 * np.pi * Time * Freq)
        Cos_C = 2.0 * np.cos(2 * np.pi * Time * Freq)

        # C_Dem = Cos_C + 1j*Sin_C
        C_Dem = Sin_C + 1j * Cos_C
        Dem = sig * C_Dem

        Dem_F = self.Filter(sig=Dem,
                            Type='lowpass',
                            Order=Order,
                            Freqs=(Filt_Freq,),
                            Fs=Fs
                            )
        Rec = int(Time_R * Fs)
        return Dem_F[Rec:]

    def Angle(self, sig):
        return np.angle(sig)

    def Abs(self, sig):
        return np.abs(sig)

    def RMS(self, sig):
        return np.sqrt(np.mean(sig ** 2))

    def Recorte(self, sig, Fs=800e3, Time_R=0.0):
        Rec = int(Time_R * Fs)
        return sig[Rec:]

    def ReadNyp(self,FilesIn,NumOfCols = 26):
        t = np.array([])
        for FilesInName in FilesIn:
            f = open(FilesInName, 'rb')
            print(FilesInName)
            while True:
                try:
                    b = np.load(f).flatten()
                    t = np.hstack((t, b)) if t.size else b
                except:
                    break
            f.close()
        return t.reshape((-1,NumOfCols))


    def ReadNypMany(self,FilesIn,NumOfFiles,NumOfCols = 26):
        t = np.array([])
        for FilesInName in range(NumOfFiles):
            f = open(FilesInName+'_' + str(n)+'.npy', 'rb')
            print(FilesInName)
            while True:
                try:
                    b = np.load(f).flatten()
                    t = np.hstack((t, b)) if t.size else b
                except:
                    break
            f.close()
        return t.reshape((-1,NumOfCols))

def NewDIC_R(FileName_Cond):
    with open(FileName_Cond, "rb") as f:
        Dic_Cond = pickle.load(f, encoding='latin1')
        RowsDic = Dic_Cond['ASICConf']['children']['RowsConf']['children']['Rows']['children']
        ColsDic = Dic_Cond['ASICConf']['children']['ColsConf']['children']['Cols']['children']

    NewDIC_R = {}
    for Value in RowsDic.keys():
        # print(Value)
        NewDIC_R["Row "+Value.split("R")[1]] = {'Enable': RowsDic[Value]['children']['Enable']['value'],
                                                'Gain': RowsDic[Value]['children']['Gain']['value'],
                                                'Offset Vector': [1] * 32}

    for Value in ColsDic.keys():
        # print(Value)
        NewDIC_R["Col "+Value.split("C")[1]] = {'OE': ColsDic[Value]['children']['Enable']['value'],
                                                'Amp_I': ColsDic[Value]['children']['AmpI']['value'],
                                                'Amp_Q': ColsDic[Value]['children']['AmpQ']['value'],
                                                'En_I': ColsDic[Value]['children']['EnableI']['value'],
                                                'En_Q': ColsDic[Value]['children']['EnableQ']['value']}
    return NewDIC_R

def DemAndStack(MyLock, BufferIndividual,ChunckSize):

    dem = np.array([])
    for samps in BufferIndividual.reshape((-1, ChunckSize)):
        r = MyLock.CalcData(samps[:, None])
        dem = np.hstack((dem, r)) if dem.size else r
    SigDem = np.abs(dem)
    return SigDem

def LoadConditions(Route=None, Device=None, FileName=None):

    try:
        if FileName is not None:
            FileName_Cond = FileName
        else:
            FileName_Cond = Route + Device + "_Cond.plk"

        with open(FileName_Cond, "rb") as f:
            Dic_Cond = pickle.load(f, encoding='latin1')
        Muestra = Dic_Cond['Muestra']
        DIC_R = Dic_Cond['DIC_R']
    except:
        FileName_Cond = Route + Device + ".pkl"
        DIC_R = NewDIC_R(FileName_Cond)
        Muestra = "25x5_M"

    return DIC_R,Muestra

def LoadConditionsOnTF(Route,Device):

    FileName = Route + Device + ".pkl"

    with open(FileName, "rb") as f:
        DevDCVals = pickle.load(f, encoding='latin1')
    DevDCVals = DevDCVals[0]

    Include_C = DevDCVals['IncludeJ']
    # Include_L = Include_C['Include']

    if 'Muestra' in Include_C.keys():
        print("Muestra ", Include_C['Muestra'])
        Muestra = Include_C['Muestra']
    else:
        Muestra = '25x5_M'

    Dic_SinI = {}
    Dic_DevDCVals = {}
    Include = []
    for index, val in Include_C['DicMuxDem'].items():
        # print(index,val)
        if val[0] < 16:
            Dic_SinI[index] = val
            Include.append(index)
            Dic_DevDCVals[index] = DevDCVals[index]

    Dic = {'DevDCVals':Dic_DevDCVals,
            'Include_C':Dic_SinI,
           'Include':Include,
           'Muestra':Muestra
           }

    return Dic


def ApplyProcessChain(sig, ProcessChain):
    sl = sig.copy()

    for Proc in ProcessChain:
        sl = Proc['function'](sl, **Proc['args'])

    return sl

def AmpsFromDem(MyLock,BufferIndividual,ChunckSize,Row,Gain,AS):

    Data = DemAndStack(MyLock=MyLock,
                       BufferIndividual=BufferIndividual,
                       ChunckSize=ChunckSize
                       )

    Data_Amps = AS.RowAndByteArrayToIamp(Buffer_Fetch=Data,
                                         Row=Row,
                                         Gain=Gain
                                         )

    return Data_Amps