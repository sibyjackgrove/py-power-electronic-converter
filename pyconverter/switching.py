"""Class for switching signals."""

import numpy as np
import math
import cmath
import scipy
import logging

from scipy import signal
from scipy.integrate import odeint,ode

import converter_utilities

import config


class SwitchingSignals:
    """
    Converter base class.
    
    Attributes:
        count (int): Number of converter objects.
        
    """
    
    count = 0 #Object count
    
    def __init__(self,fm):
        """Creates an instance of `SwitchingSignals`.
        
        Args:
           fsw (float): Switching frequency in Hz.
           fm (float): Fundamental frequency of output waveform in Hz.
           
        Raises:
          ValueError: To be added
        
        """
        
        self.update_modulating_waveform(1.0,fm)
        self.update_carrier_waveform(1.0,config.DEFAULT_fsw)
   
    
    @property                         #Decorator used for auto updating
    def mf(self,):
        """Frequency modulation index"""
        
        return self.fsw/self.fm
    
    @property                         #Decorator used for auto updating
    def ma(self,):
        """Amplitude modulation index"""
        
        return self.Am/self.Ac
    
    def update_modulating_waveform(self,Am=1.0,fm = 50.0):
        """Update modulating amplitude and frequency."""
        
        self.Am=Am
        self.fm = fm
        self.Tm= 1/self.fm
    
    def update_carrier_waveform(self,Ac=1.0,fsw = 10000.0):
        """Update switching/carrier amplitude and frequency."""
                
        self.Ac = Ac
        self.fsw = fsw #(fsw = fcarrier)Switching and carrier frequency are same    
        self.Tsw = 1/self.fsw
    
    def show_spec(self):
        """Show waveforms."""
        
        print('fswitch:{:.2f},Tswitch:{:.5f}'.format(self.fsw,self.Tsw))
        print('fmodulating:{:.2f},Tmodulating:{:.5f}'.format(self.fm,self.Tm))
        print('ma:{:.3f},mf:{:.3f}'.format(self.ma,self.mf))
    
    def duty_cycle(self,t):
        """Create a sinusoidal PWM signal."""
        
        modulating_signal = self.sinosoid(self.Am,self.fm,t)
        
        return {'modulating':modulating_signal} 
        
    def sinosoidalPWM(self,t):
        """Create a sinusoidal PWM signal."""
        
        modulating_signal = self.sinosoid(self.Am,self.fm,t)
        carrier_signal = self.sawtooth(self.Ac,self.fsw,t)
        
        switching_signal = self.comparator(modulating_signal,carrier_signal)

        return {'switching':switching_signal,'carrier':carrier_signal,'modulating':modulating_signal} 
    
    def square_wave(self,t):
        """Create a square wave signal for ."""
        
        modulating_signal = self.sinosoid(self.Am,self.fm,t)

        #if modulating_signal>0.0:
        #    switching_signal = True
        #else:
        #    switching_signal = False
            
        switching_signal = self.comparator(modulating_signal,0.0)    

        return {'switching':switching_signal,'modulating':modulating_signal}
  
    def comparator(self,modulating_signal,carrier_signal):
        """Compare modulating signal with carrier signal."""
        
        if modulating_signal>carrier_signal:
            switching_signal = True
        else:
            switching_signal = False
        return switching_signal
    
    def square(self,amplitude,frequency,duty,t):
        """Create a sinusoid time series."""
    
        return amplitude*signal.square(2*np.pi*frequency*t,duty)
    
    def sinosoid(self,amplitude,frequency,t):
        """Create a sinusoid time series."""
        
        return amplitude*np.sin(2*np.pi*frequency*t)

    def sawtooth(self,amplitude,frequency,t):
        """Create a sinusoid time series."""
        
        return amplitude*signal.sawtooth(2*np.pi*frequency*t, width=0.5)

            
    def check_fsw(self):
        """Check if fsw is sufficient."""
        
        assert self.fsw/self.fm >= 6.0, 'Switching freq {} should be atleast 6 times greater than modulating signal frequency {}'.format(self.fsw,self.fm)
    
    
    def simulate_switching(self,tf,modulating_signal = 'sinusoid'):
        """Simulate converter operation."""
        
        Am = 1 #Amplitude of modulating signal
        fm = 10 #Frequency of modulating signal
        dt = self.dt_calc()
        
        carrier_signal_t = []
        switching_signal_t = []
        
        t_t = np.arange(0,tf,dt) #np.linspace(0, 2, 1000)
        if modulating_signal == 'sinusoid':
            modulating_signal_t = self.sinosoid(Am,fm,t_t)
        else:
            modulating_signal_t = self.square(0.5,fm,1.0,t_t)
        
        for t,modulating_signal in zip(t_t,modulating_signal_t):
            switching_signal,carrier_signal = self.sinosoidalPWM(t,modulating_signal)
            switching_signal_t.append(switching_signal)
            carrier_signal_t.append(carrier_signal)    
            print('Modulating signal:{:.3f},Carrier signal:{:.3f},Switching signal:{:.3f}'.format(modulating_signal,carrier_signal,switching_signal))
        
        plot_signal(t_t, modulating_signal_t,'modulating signal',False)
        plot_signal(t_t, carrier_signal_t,'carrier signal',False)
        plot_signal(t_t, switching_signal_t,'switching signal',True)