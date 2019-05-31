'''
Fourier analysis module
'''

import numpy as np
from collections import deque

class aswift():
    '''
    Alpha sliding windowed infinite fourier transform (aSWIFT). The aSWIFT represents the fourier
    transform at t=n as a function of the transform at t=n-1. The aSWIFT should be used
    whenever the window overlap is large compared to the window length.

    In order to normalize, the output should be divided by (taus-tauf).

    Parameters
    ----------
    tau_s : float (seconds)
        Slow time constant
    tau_f : float (seconds)
        Fast time constant
    f : float or array_like (Hz)
        Center frequency(ies) of transform
    fs : float (Hz)
        Sampling frequency

    See Also
    --------
    sdft
    swift
    '''

    def __init__(self, tau_s, tau_f, f, fs):
        self.__tau_s = tau_s
        self.__tau_f = tau_f
        self.__f = f
        self.__fs = fs
              
        self.__slow = swift(self.tau_s,self.f,self.fs)
        self.__fast = swift(self.tau_f,self.f,self.fs)

    def slide(self, x):
        '''
        Slide forward N samples, where N is the length of x.

        Parameters
        ----------
        x : float or array_like
            Sample or array of samples

        Returns
        -------
        Xf : complex
            Numpy.array of complex numbers at frequencies f
        '''
        return self.slow.slide(x) - self.fast.slide(x)
        
    @property
    def tau_s(self):
        '''float : tau slow (s)'''
        return self.__tau_s
    @property
    def tau_f(self):
        '''float : tau fast (s)'''
        return self.__tau_f
    @property
    def f(self):
        '''
        float : center frequency (Hz)
        '''
        return self.__f
    @property
    def fs(self):
        '''float : sampling frequency (Hz)'''
        return self.__fs
    @property
    def slow(self):
        return self.__slow
    @property
    def fast(self):
        return self.__fast
    
class swift():
    '''
    Sliding windowed infinite fourier transform (SWIFT). The SWIFT represents the fourier
    transform at t=n as a function of the transform at t=n-1. The SWIFT should be used
    whenever the window overlap is large compared to the window length.

    In order to normalize, the output should be divided by tau.

    Parameters
    ----------
    tau : float or array_like (seconds)
        Exponential window time constant(s)
    f : float or array_like (Hz)
        Center frequency(ies) of transform
    fs : float (Hz)
        Sampling frequency
    
    See Also
    --------
    aswift
    sdft
    '''

    def __init__(self,tau,f,fs):
        tau,f,fs = self.__paramcheck(tau,f,fs)
        
        self.__tau  = tau
        self.__f    = f
        self.__fs   = fs
        self.__ntau = self.tau*self.fs
        
        self.Xf = np.zeros(len(self.f),dtype=complex)
        
        self.e = np.exp(2j*np.pi*self.f/self.fs)*np.exp(-1./self.ntau)

    def slide(self,x):
        '''
        Slide forward N samples, where N is the length of x.

        Parameters
        ----------
        x : float or array_like
            Sample or array of samples

        Returns
        -------
        Xf : complex
            Numpy.array of complex numbers at frequencies f
        '''
        x = np.asarray(x).reshape(-1)
        for i in range(len(x)):
            self.Xf = self.e*self.Xf + x[i]
        return self.Xf

    def __paramcheck(self,tau,f,fs):
        '''
        Ensures all parameters are the correct type and in range.
        '''
        try: tau = np.asarray(tau).reshape(-1)
        except: raise TypeError('data type not understood. tau must be array_like')
        try: f = np.asarray(f).reshape(-1)
        except: raise TypeError('data type not understood. f must be array_like')
        try: fs = float(fs)
        except: raise TypeError('data type not understood. fs must be float')

        if fs <= 0: raise ValueError('fs must be > 0')
        if np.any(tau <= 0): raise ValueError('tau must be > 0')
        if np.any(f <= 0): raise ValueError('f must be > 0')
        
        return tau,f,fs

    @property
    def tau(self):
        return self.__tau
    @property
    def ntau(self):
        return self.__ntau
    @property
    def f(self):
        return self.__f
    @property
    def fs(self):
        return self.__fs
    
