#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 11:12:23 2018

@author: aguimera
"""

#def PlotPSD(self, Time, nFFT=2**17, FMin=None, Resamp=False):
#
#    if not self.FigFFT or not plt.fignum_exists(self.FigFFT.number):
#        self.FigFFT, self.AxFFT = plt.subplots()
#
#    PSD = {}
#    for sl in self.Slots:
#        sig = sl.GetSignal(Time, Resamp=Resamp)
#        if FMin:
#            nFFT = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/FMin))+1)) 
#
#        ff, psd = signal.welch(x=sig, fs=sig.sampling_rate,
#                               window='hanning',
#                               nperseg=nFFT, scaling='density', axis=0)
#        slN = sl.SigName
#        PSD[slN] = {}
#        PSD[slN]['psd'] = psd
#        PSD[slN]['ff'] = ff
#        self.AxFFT.loglog(ff, psd, label=sl.DispName)
#
#    self.AxFFT.set_xlabel('Frequency [Hz]')
#    self.AxFFT.set_ylabel('PSD [V^2/Hz]')
#    self.AxFFT.legend()
#    return PSD

#
#    def PlotHist(self, Time, Resamp=False, Binds=250):
#        fig, Ax = plt.subplots()
#
#        for sl in self.Slots:
#            Ax.hist(sl.GetSignal(Time, Resamp=Resamp),
#                    Binds,
#                    alpha=0.5)
#            Ax.set_yscale('log')
#
#    def PlotPSD(self, Time, nFFT=2**17, FMin=None, Resamp=False):
#
#        if not self.FigFFT or not plt.fignum_exists(self.FigFFT.number):
#            self.FigFFT, self.AxFFT = plt.subplots()
#
#        PSD = {}
#        for sl in self.Slots:
#            sig = sl.GetSignal(Time, Resamp=Resamp)
#            if FMin:
#                nFFT = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/FMin))+1)) 
#
#            ff, psd = signal.welch(x=sig, fs=sig.sampling_rate,
#                                   window='hanning',
#                                   nperseg=nFFT, scaling='density', axis=0)
#            slN = sl.SigName
#            PSD[slN] = {}
#            PSD[slN]['psd'] = psd
#            PSD[slN]['ff'] = ff
#            self.AxFFT.loglog(ff, psd, label=sl.DispName)
#
#        self.AxFFT.set_xlabel('Frequency [Hz]')
#        self.AxFFT.set_ylabel('PSD [V^2/Hz]')
#        self.AxFFT.legend()
#        return PSD


#    def PlotEventAvg(self, (EventRec, EventName), Time, TimeWindow,
#                     OverLap=True, Std=False, Spect=False, Resamp=False):
#
#        ft, Axt = plt.subplots()
#
#        for sl in self.Slots:
#            avg = np.array([])
#            if Spect:
#                ft, (Ax, AxS) = plt.subplots(2, 1, sharex=True)
#            else:
#                ft, Ax = plt.subplots()
#
#            if Resamp:
#                Fs = sl.ResampleFs.magnitude
#            else:
#                Fs = sl.Signal().sampling_rate.magnitude
#
#            Ts = 1/Fs
#            nSamps = int((TimeWindow[1]-TimeWindow[0])/Ts)
#            t = np.arange(nSamps)*Ts*pq.s + TimeWindow[0]
#
#            etimes = EventRec.GetEventTimes(EventName, Time)
#            for et in etimes:
#                start = et+TimeWindow[0]
#                stop = et+TimeWindow[1]
#
#                if sl.Signal().t_start < start and sl.Signal().t_stop > stop:
#                    st = sl.GetSignal((start, stop), Resamp=Resamp)[:nSamps]
#                    try:
#                        avg = np.hstack([avg, st]) if avg.size else st
#                        if OverLap:
#                            Ax.plot(t, st, 'k-', alpha=0.1)
#                    except:
#                        print 'Error', nSamps, et, avg.shape, st.shape
#
#            print EventName, 'Counts', len(etimes)
#
#            MeanT = np.mean(avg, axis=1)
#            Ax.plot(t, MeanT, 'r-')
#            if Std:
#                StdT = np.std(avg, axis=1)
#                Ax.fill_between(t, MeanT+StdT, MeanT-StdT,
#                                facecolor='r', alpha=0.5)
#
#            ylim = Ax.get_ylim()
#            Ax.plot((0, 0), (1, -1), 'g-.', alpha=0.7)
#            Ax.set_ylim(ylim)
#            Ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
#
#            plt.title(sl.DispName)
#            if Spect:
#                nFFT = int(2**(np.around(np.log2(Fs/sl.SpecFmin))+1))
#                noverlap = int((Ts*nFFT - sl.SpecTimeRes)/Ts)
#                TWindOff = (nFFT * Ts)/8
#
#                f, tsp, Sxx = signal.spectrogram(MeanT, Fs,
#                                                 window='hanning',
#                                                 nperseg=nFFT,
#                                                 noverlap=noverlap,
#                                                 scaling='density',
#                                                 axis=0)
#
#                finds = np.where(f < sl.SpecFmax)[0][1:]
#                print Sxx.shape
#                r, c = Sxx.shape
#                S = Sxx.reshape((r, c))[finds][:]
#                pcol = AxS.pcolormesh(tsp + TimeWindow[0].magnitude + TWindOff,
#                                      f[finds],
#                                      np.log10(S),
#                                      vmin=np.log10(np.max(S))+sl.SpecMinPSD,
#                                      vmax=np.log10(np.max(S)),
#                                      cmap=sl.SpecCmap)
#                f, a = plt.subplots(1, 1)
#                f.colorbar(pcol)
#
#            Axt.plot(t, np.mean(avg, axis=1), label=sl.DispName)
#            ft.canvas.draw()
#
#        Axt.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
#        Axt.legend()
#        ft.canvas.draw()
#        plt.show()
#
#    def PlotPSDSNR(self, (EventRec, EventName), TDelay, TEval, DevDCVals,
#                   Time=None, nFFT=2**17):
#
#        etimes = EventRec.GetEventTimes(EventName, Time)
#
#        DevACVals = {}
#
#        for sl in self.Slots:
#            fig, (ax1, ax2)  = plt.subplots(1,2)
#            psd = np.array([])
#            DCVals = DevDCVals[sl.SigName]
#            for ne, et in enumerate(etimes):
#                start = et+TDelay
#                stop = et+TDelay+TEval
#
#                sig = sl.GetSignal((start, stop), Resamp=False)
#                fpsd, npsd = signal.welch(x=sig, fs=sig.sampling_rate,
#                                          window='hanning',
#                                          nperseg=nFFT, scaling='density', axis=0)                
#                psd = np.hstack([psd, npsd]) if psd.size else npsd
#
##                Flin = fpsd[1:].magnitude
##                Flog = np.logspace(np.log10(Flin[0]),
##                                   np.log10(Flin[-1]), 100)
##                Flog = np.round(Flog, 9)
##                Flin = np.round(Flin, 9)
##                intpsd = interpolate.interp1d(Flin, npsd[1:].transpose())(Flog)
##                a, b, _ = FETAna.noise.FitNoise(Flog,
##                                                intpsd, Fmin=150, Fmax=5e3)
#
#                ax1.loglog(fpsd, npsd)
##                ax2.loglog(Flog, intpsd/FETAna.noise(Flog, a, b))
#
#            psdD = {'Vd0': psd.transpose()}
#            ACVals = {'PSD': psdD,
#                      'gm': None,
#                      'Vgs': DCVals['Vgs'],
#                      'Vds': DCVals['Vds'],
#                      'Fpsd': fpsd.magnitude,
#                      'Fgm': None,
#                      'ChName': sl.SigName,
#                      'Name': sl.DispName,
#                      'GMPoly': DCVals['GMPoly'],
#                      'IdsPoly': DCVals['IdsPoly'],
#                      'Ud0': DCVals['Ud0'],
#                      'IsOK': DCVals['IsOK'],
#                      'DateTime': DCVals['DateTime']}
#            DevACVals[sl.DispName] = ACVals
#            fig.canvas.draw()
#            plt.show()
#
#        FETAna.InterpolatePSD(DevACVals)
#        FETAna.FitACNoise(DevACVals, Fmin=150, Fmax=5e3, IsOkFilt=False)
#        pltPSD = FETplt.PyFETPlot()
#        pltPSD.AddAxes(('PSD', 'NoA', 'NoB'))
#        pltPSD.AddLegend()
#        pltPSD.PlotDataCh(DevACVals)
#
#        pickle.dump(DevACVals, open('test.pkl', 'wb'))
#
#        return DevACVals
