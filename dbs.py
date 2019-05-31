import numpy as np

from swift import aswift

class DBS(object):
    def __init__(self, dt=1e-3, stim_amp=1, width=60, tstart=0):
        '''
        Initialize a DBS object

        Parameters
        ----------
        | stim_amp : float
        |     Stimulus puse amplitude (mA)
        | width : float
        |     Stimulus puse width (us)
        '''
        self._dt     = dt
        self._stim_amp    = stim_amp
        self._width  = width
        self._tstart = tstart
        self._calc_charge()
        
    def _calc_charge(self):
        self._charge = max(0, self.stim_amp * self.width * 1e-6)
        
    # Setters/Getters
    #----------------
    @property
    def stim_amp(self):
        return self._stim_amp
    @stim_amp.setter
    def stim_amp(self,x):
        self._stim_amp = x
        self._calc_charge()
        
    @property
    def width(self):
        return self._width
    @width.setter
    def width(self,x):
        self._width = x
        self._calc_charge()

    # Getters only
    #-------------
    @property
    def dt(self):
        return self._dt

    @property
    def tstart(self):
        return self._tstart
    
    @property
    def charge(self):
        return self._charge #Charge in mC

    
class cDBS(DBS):
    def __init__(self, f, dt=1e-3, stim_amp=1, width=60, tstart=0):
        '''
        Continuous DBS object.

        Parameters
        ----------
        | dt : float (s)
        |     integration timestep
        | f : float (Hz)
        |     stimulation frequency
        | stim_amp : float (mA)
        |     stimulus pulse ampltiude
        | width : float (us)
        |     stimulus pulse width
        | tstart : float (s)
        |     stimulation start time
        '''
        super(cDBS,self).__init__(dt,stim_amp,width,tstart)

        self._f = f
        
        self._steps_per_pulse    = int(round(1./self._f/self._dt))
        self._pulse_counter      = self._steps_per_pulse

        self._n_start = self.tstart / self.dt
        self._i = 0
        
    def advance(self):
        stim = False
        if self._i < self._n_start:
            self._i += 1
        else:
            if self._pulse_counter >= self._steps_per_pulse:
                self._pulse_counter = 0
                stim = True
            self._pulse_counter += 1
        return stim * self.charge

    @property
    def f(self):
        return self._f


class pDBS(DBS):
    def __init__(self, f, tau_s, tau_f, phase_thr=0, power_thr=-np.inf, ref_period=0.3, dt=1e-3, stim_amp=1, width=60, tstart=0):
        super(pDBS,self).__init__(dt,stim_amp,width,tstart)

        self._f          = f
        self._tau_s      = tau_s
        self._tau_f      = tau_f
        self._phase_thr  = phase_thr
        self._power_thr  = power_thr
        self._ref_period = ref_period

        self._aswift = aswift(tau_s = self.tau_s,
                              tau_f = self.tau_f,
                              f     = self.f,
                              fs    = 1./self.dt)

        self._X          = 0
        self._amp        = 0
        self._phase      = 0

        self._shift_phase = 0
        self._last_shift_phase = np.inf
        
        self._update()

    def _update(self):
        self._n_ref = 1/(self.f*self.dt) * self.ref_period
        self._i_ref = self._n_ref

    def advance(self,x):
        self._X = self.aswift.slide(x)
        self._amp,self._phase = np.abs(self._X), np.angle(self._X)
        self._amp /= (self.tau_s - self.tau_f) / self.dt
        with np.errstate(divide='ignore'):
            self._amp = 10*np.log10(self._amp**2)

        self._last_shift_phase = self._shift_phase
        self._shift_phase = (self._phase - self._phase_thr + np.pi) % (2*np.pi) - np.pi
        
        stim=False
        if self._last_shift_phase < 0 <= self._shift_phase:
        #if self._last_phase - self.phase_thr < 0 <= self._phase - self.phase_thr:
            if self._i_ref >= self._n_ref:
                if self._amp >= self._power_thr:
                    stim=True
                    self._i_ref = 0

        self._i_ref += 1
        
        return self.charge * stim

    # Mutable Properties
    #-------------------
    @property
    def phase_thr(self):
        return self._phase_thr
    @phase_thr.setter
    def phase_thr(self,x):
        self._phase_thr = x

    @property
    def power_thr(self):
        return self._power_thr
    @power_thr.setter
    def power_thr(self,x):
        self._power_thr = x
        
    # Immutable Properties
    #---------------------
    @property
    def f(self):
        return self._f
    @property
    def tau_s(self):
        return self._tau_s
    @property
    def tau_f(self):
        return self._tau_f
    @property
    def aswift(self):
        return self._aswift
    @property
    def amp(self):
        return self._amp
    @property
    def ref_period(self):
        return self._ref_period
    @property
    def phase(self):
        return self._phase

