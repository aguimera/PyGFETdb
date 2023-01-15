import numpy as np
# from PyQt5 import Qt
from scipy import signal
from fractions import Fraction

class SosFilter:
    def __init__(self, **kwargs):
        self.sos = signal.butter(output='sos', **kwargs)
        self.zi = signal.sosfilt_zi(self.sos)

    def CalcData(self, data):
        fdat, self.zi = signal.sosfilt(sos=self.sos,
                                       x=data,
                                       axis=-1,
                                       zi=self.zi)
        return fdat


class SyncCarriers:
    def __init__(self, Fs, Fc, BufferSize, GenFunct):
        self.ChunckSize = BufferSize
        self.vFc = np.ones(BufferSize)
        self.EndInd = 0
        self.Fc = Fc
        self.Fs = Fs

        Ts = 1 / Fs
        FcPoints = Fraction(Fs / Fc).numerator
        TcGen = np.arange(0, Ts * FcPoints, Ts)
        vFcSingle = GenFunct(2 * np.pi * Fc * TcGen)
        nFcBlocks = int(self.ChunckSize / FcPoints)
        self.CstBlock = np.tile(vFcSingle, nFcBlocks)
        self.CstBlockSize = self.CstBlock.size

    def GetNextBlock(self):
        StartInd = self.CstBlockSize - self.EndInd
        self.vFc[:StartInd] = self.CstBlock[self.EndInd:]
        self.EndInd = self.ChunckSize - StartInd
        self.vFc[StartInd:] = self.CstBlock[:self.EndInd]

        return self.vFc


class LockIn:
    def __init__(self, FsIn, Fc, BufferSize, Filter, DownFactor, nRow=0, Index=0):
        self.nRow = nRow
        self.Index = Index
        self.Fs = FsIn
        self.Ts = 1 / FsIn
        self.BufferSize = BufferSize
        self.timeT = np.arange(start=0,
                               stop=BufferSize * self.Ts,
                               step=self.Ts)

        self.sinT = SyncCarriers(Fs=FsIn,
                                 Fc=Fc,
                                 BufferSize=BufferSize,
                                 GenFunct=np.sin)
        self.cosT = SyncCarriers(Fs=FsIn,
                                 Fc=Fc,
                                 BufferSize=BufferSize,
                                 GenFunct=np.cos)

        self.FilterReal = SosFilter(**Filter)
        self.FilterImag = SosFilter(**Filter)

        self.DownSampInds = np.arange(start=0,
                                      stop=BufferSize,
                                      step=DownFactor)
        self.FsOut = self.Fs / DownFactor

    def CalcData(self, NewData):
        d = NewData[:, self.nRow]
        r = self.FilterReal.CalcData(d * self.cosT.GetNextBlock())[self.DownSampInds]
        i = self.FilterImag.CalcData(d * self.sinT.GetNextBlock())[self.DownSampInds]
        return r + i * 1j

