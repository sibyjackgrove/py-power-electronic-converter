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


class Signals:
    """
    Signals base class.
    
    Attributes:
        count (int): Number of Signals objects.
        
    """
    
    count = 0 #Object count
    
    def __init__(self):
        """Creates an instance of `SwitchingSignals`.
        
        """
        
    def square(self,amplitude,frequency,duty,t):
        """Create a sinusoid time series."""
    
        return amplitude*signal.square(2*np.pi*frequency*t,duty)
    
    def sinosoid(self,amplitude,frequency,phase,t):
        """Create a sinusoid time series."""
        
        return amplitude*np.sin(2*np.pi*frequency*t+phase)

    def sawtooth(self,amplitude,frequency,t):
        """Create a sinusoid time series."""
        
        return amplitude*signal.sawtooth(2*np.pi*frequency*t, width=0.5)
    
    def three_phase_sinosoid(self,amplitude,frequency,phase,t):
        """Create a 3 phase sinusoid time series."""
        
        phase_a = phase
        phase_b = -(2/3)*np.pi + phase
        phase_c = (2/3)*np.pi + phase
        
        va = self.sinosoid(amplitude,frequency,phase_a,t)
        vb = self.sinosoid(amplitude,frequency,phase_b,t)
        vc = self.sinosoid(amplitude,frequency,phase_c,t)
        
        return va,vb,vc


    
class SwitchingSignals(Signals):
    """
    Converter base class.
    
    Attributes:
        count (int): Number of converter objects.
        
    """
    
    count = 0 #Object count
    
    def __init__(self,fsw=config.DEFAULT_fsw,fm=config.DEFAULT_fm,phasem=0.0):
        """Creates an instance of `SwitchingSignals`.
        
        Args:
           fsw (float): Switching frequency in Hz.
           fm (float): Frequency of modulating waveform in Hz.
           
        Raises:
          ValueError: To be added
        
        """
        
        self.update_modulating_waveform(config.DEFAULT_Am,fm,phasem)
        self.update_carrier_waveform(config.DEFAULT_Asw,fsw)
   
    
    @property                         #Decorator used for auto updating
    def mf(self,):
        """Frequency modulation index"""
        
        return self.fsw/self.fm
    
    @property                         #Decorator used for auto updating
    def ma(self,):
        """Amplitude modulation index"""
        
        return self.Am/self.Ac
    
    def update_modulating_waveform(self,Am,fm,phasem):
        """Update modulating amplitude and frequency."""
        
        self.Am=Am
        self.fm = fm
        self.phasem = phasem
        self.Tm= 1/self.fm
    
    def update_carrier_waveform(self,Ac,fsw):
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
        
        modulating_signal = self.sinosoid(self.Am,self.fm,self.phasem,t)
        carrier_signal = self.sawtooth(self.Ac,self.fsw,t)
        
        switching_signal = self.comparator(modulating_signal,carrier_signal)

        return {'switching':switching_signal,'carrier':carrier_signal,'modulating':modulating_signal} 
    
    def square_wave(self,t):
        """Create a square wave signal for ."""
        
        modulating_signal = self.sinosoid(self.Am,self.fm,self.phasem,t)

        switching_signal = self.comparator(modulating_signal,0.0)    

        return {'switching':switching_signal,'modulating':modulating_signal}
  
    def comparator(self,modulating_signal,carrier_signal):
        """Compare modulating signal with carrier signal."""
        
        if modulating_signal>carrier_signal:
            switching_signal = True
        else:
            switching_signal = False
        return switching_signal    
    
            
    def check_fsw(self):
        """Check if fsw is sufficient."""
        
        assert self.fsw/self.fm >= 6.0, 'Switching freq {} should be atleast 6 times greater than modulating signal frequency {}'.format(self.fsw,self.fm)
    
    def simulate_switching(self,tf,dt,switching_signal_method,avoid_signals=['carrier_signal']):
        """Simulate converter operation."""
        
        t_t = np.arange(0.0,tf,dt)
        switching_signal_t = []
        carrier_signal_t = []
        modulating_signal_t = []
        
        for t in t_t:
            signal_dict = switching_signal_method(t)
            switching_signal_t.append(signal_dict['switching'])
            carrier_signal_t.append(signal_dict['carrier'])
            modulating_signal_t.append(signal_dict['modulating'])
        
        if 'switching_signal' not in avoid_signals:
            converter_utilities.plot_signal(t_t, switching_signal_t,'switching signal')
        if 'carrier_signal' not in avoid_signals:
            converter_utilities.plot_signal(t_t, carrier_signal_t,'carrier signal')
        if 'modulating_signal' not in avoid_signals:
            converter_utilities.plot_signal(t_t, modulating_signal_t,'modulating signal')
        
        converter_utilities.show_plot_signal()
    
    def simulate_switching1(self,tf,modulating_signal = 'sinusoid'):
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
        
        
class InverterSignals(SwitchingSignals):
    """
    SwitcInverter class.
    
    Attributes:
        ():
        
    """
    
    def __init__(self,fsw=config.DEFAULT_fsw,fm=config.DEFAULT_fm,switching_signal_type=config.DEFAULT_switching_signal):
        """Creates an instance of `InverterSignals`.
        
        Args:
           fsw (float): Switching frequency in Hz.
           fm (float): Frequency of modulating waveform in Hz.
           
        Raises:
          ValueError: To be added
        
        """
        
        super().__init__(fsw,fm)  #Initialize converter class (base class)
        self.switching_signal_type = switching_signal_type
    
    def switching_signal_single_phase(self,t):
        """Calculate switching signal."""
        
        signal = self.switching_signal_calc(t)
        
        return signal['switching']        
        
    def modulating_signal_single_phase(self,t):
        """Calculate switching signal."""
        
        signal = self.modulating_signal_calc(signals,t)
        
        return signal['modulating']       
    
    def phasor_signal_single_phase(self,signals,t):
        """Calculate switching signal."""
        
        signal = 1.0 + 0.0*1j
        
        return signal       
    
    def switching_signal_calc(self,t):
        """Calculate switching signal."""
        
        if self.switching_signal_type=='sinePWM':
            signal_dict = self.sinosoidalPWM(t)
        elif self.switching_signal_type=='square_wave':     
            signal_dict = self.square_wave(t)
        else:
            raise ValueError(f'{self.switching_signal_type} is not valid switching signal type!')
        #print(signal_dict)
        return signal_dict
    
    def modulating_signal_calc(self,t):
        """Calculate switching signal."""
        
        signal_dict = self.duty_cycle(t)
        
        return signal_dict
