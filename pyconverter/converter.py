"""Class for converter."""

import numpy as np
import math
import cmath
import scipy
import logging

from scipy import signal
from scipy.integrate import odeint,ode

from converter_utilities import plot_signal

class SwitchingSignals:
    """
    Converter base class.
    
    Attributes:
        count (int): Number of converter objects.
        
    """
    
    count = 0 #Object count
    
    def __init__(self,fsw,fm):
        """Creates an instance of `SwitchingSignals`.
        
        Args:
           fsw (float): Switching frequency in Hz.
           fm (float): Fundamental frequency of output waveform in Hz.
           
        Raises:
          ValueError: To be added
        
        """
        
        self.update_modulating_waveform(1.0,fm)
        self.update_carrier_waveform(1.0,fsw)
   
    
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
    
    def show_waveforms(self):
        """Show waveforms."""
        
        print('fswitch:{:.2f},Tswitch:{:.5f}'.format(self.fsw,self.Tsw))
        print('fmodulating:{:.2f},Tmodulating:{:.5f}'.format(self.fsw,self.Tsw))
        print('ma:{:.3f},mf:{:.3f}'.format(self.ma,self.mf))
        
        print('Time step:{:.5f}'.format(self.dt_calc()))
    
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
    
    def dt_calc(self):
        """Calcuate time step."""
        
        return 1/(10*self.fsw) #Time steps in switching simulation
    
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
    
        

class PowerElectronicConverter:
    """
    Converter base class.
    
    Attributes:
        count (int): Number of converter objects.
        
    """
    
    count = 0 #Object count
    
    def __init__(self):
        """Creates an instance of `Converter`.
        
        Args:
           fsw (float): Switching frequency in Hz.
           
        Raises:
          ValueError: If parameters corresponding to `Sinverter_rated` are not available.
        
        """
                
        PowerElectronicConverter.count = PowerElectronicConverter.count+1 #Increment count to keep track of number of converter model instances
        
        self.name = 'converter_'+str(PowerElectronicConverter.count)  #Generate a name for the instance
       
    
    def calc_primary(self,signal):
        """Calculate the primary switch."""
        
        assert isinstance(signal,bool), 'Switching signal must be boolean.'
        Sprimary = int(signal)   

        return Sprimary

    def calc_complimentary(self,signal):
        """Calculate the complimentary."""
        
        assert isinstance(signal,bool), 'Switching signal must be boolean.'
        Scomplimentary = int(not signal)       
    
        return Scomplimentary    

        
class PowerElectronicInverter(PowerElectronicConverter):
    """
    Inverter class.
    
    Attributes:
        ():
        
    """
    
    Rf = 0.01
    Lf = 2.0e-3
    
    Rload = 1.0
    
    def __init__(self,Vdc,signals,signal_type='square'):
        """Creates an instance of `Converter`.
        
        Args:
           Vdc (float): DC link voltage.
           
        Raises:
          ValueError: To be added.
        
        """
        
        super().__init__()  #Initialize converter class (base class)
        
        self.update_Vdc(Vdc)
        self.signals = signals
        self.signal_type='square'
    
    @property                         #Decorator used for auto updating
    def y(self):
        """List of initial states"""
        
        return  [self.ia, 0.0]
    
    def update_Vdc(self,Vdc):
        """Update DC link voltage."""
        
        self.Vdc = Vdc
    
    def switching_signal_calc(self,t):
        """Calculate switching signal."""
        
        if self.signal_type=='PWM':
            signal_dict = self.signals.sinosoidalPWM(t)
        elif self.signal_type=='square':     
            signal_dict = self.signals.square_wave(t)
    
        return signal_dict
        
    def half_bridge_ideal(self,Vdc,S1):
        """Simulates a bridge in inverter"""
        
        self.update_Vdc(Vdc)

        S11 = self.calc_primary(S1)
        S12 = self.calc_complimentary(S1)
        
        assert S11+S12 == 1, 'S11 and S12 switches cannot be both ON or OFF at the same time in ideal half bridge.'

        #print('S11:{},S12:{}'.format(S11,S12))
        Van = (S11 - S12)*(self.Vdc/2)
        
        #print('Van:{}'.format(Van))

        return Van
    
    def three_phase_full_bridge_ideal(Vdc,S1,S2,S3):
        """Simulates a bridge in inverter"""

        S11 = calc_primary(S1)
        S12 = calc_complimentary(S1)

        S21 = calc_primary(S2)
        S22 = calc_complimentary(S2)

        S31 = calc_primary(S3)
        S32 = calc_complimentary(S3)
        
        assert S11+S12 == 1, 'S11 and S12 switches cannot be both ON or OFF at the same time in ideal half bridge.'
        assert S21+S22 == 1, 'S21 and S22 switches cannot be both ON or OFF at the same time in ideal half bridge.'
        assert S31+S32 == 1, 'S31 and S32 switches cannot be both ON or OFF at the same time in ideal half bridge.'        

        print('S1:{},S2:{},S3:{}'.format(S11,S21,S31))

        Vno =  (self.Vdc/6)*(2*S11+2*S21+2*S31-3)

        Van = (self.Vdc/2)*(S11-S12)-Vno
        Vbn = (self.Vdc/2)*(S21-S22)-Vno
        Vcn = (self.Vdc/2)*(S31-S32)-Vno

        #Van = (2*S11 - S21 - S31)*(Vdc/3)
        #Vbn = (2*S21 - S11 - S31)*(Vdc/3)
        #Vcn = (2*S31 - S21 - S11)*(Vdc/3)
        print('Vno:{},Van+Vbn+Vcn:{}'.format(Vno,Van+Vbn+Vcn))

        print('Van:{},Vbn:{},Vcn:{}'.format(Van,Vbn,Vcn))
        print('Vab:{},Vbc:{},Vca:{}'.format(Van-Vbn,Vbn-Vcn,Vcn-Van))

        return Van,Vbn,Vcn

    def ODE_model(self,y,t):
        """ODE model of inverter branch."""
        
        self.ia,dummy = y   # unpack current values of y
        #Am = 1 #Amplitude of modulating signal
       
        Vdc = 100.0 #Get DC link voltage
        signals = self.switching_signal_calc(t)
        switching_signal = signals['switching']
        
        self.vta = self.half_bridge_ideal(Vdc,switching_signal)
        
        dia = (1/self.Lf)*(-self.Rf*self.ia + -self.Rload*self.ia + self.vta)
        
        result =  [dia,dummy]
        
        return np.array(result)
    
    def simulate_inverter(self,tf):
        """Simulate an inverter.
        
        Args:
           fm (float): Frequency of modulating signal in Hz.
           tf (float): Simulation time in seconds.
        """
        
        self.dt = self.signals.dt_calc()
        self.t_t = np.arange(0,tf,self.dt)
        print('Number of timsteps:{}'.format(len(self.t_t)))
        y0 = [0.0,0.0]
        solution,infodict = odeint(self.ODE_model,y0,self.t_t,full_output=1,printmessg=True,atol=1e-6,rtol=1e-6)
        
        self.ia_t = solution[:,0]
        self.va_t = self.ia_t*self.Rload
        
        plot_signal(self.t_t,self.ia_t,'ia',showPlot=False)
        plot_signal(self.t_t,self.va_t,'va',showPlot=True)
        
        
        return solution
    
    def show_states(self):
        """Show states."""
        
        print('Inverter states:{}'.format(self.y))