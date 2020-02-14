"""Class for converter."""

import numpy as np
import math
import cmath
import scipy
import logging

from scipy import signal
from scipy.integrate import odeint,ode

#from converter_utilities import plot_signal, plot_FFT
import converter_utilities
import config
      

class PowerElectronicConverter:
    """
    Converter base class.
    
    Attributes:
        count (int): Number of converter objects.
        
    """
    
    count = 0 #Object count
    
    def __init__(self,model_type,signal_type):
        """Creates an instance of `Converter`.
        
        Args:
           fsw (float): Switching frequency in Hz.
           
        Raises:
          ValueError: If parameters corresponding to `Sinverter_rated` are not available.
        
        """
                
        PowerElectronicConverter.count = PowerElectronicConverter.count+1 #Increment count to keep track of number of converter model instances
        
        self.name = 'converter_'+str(PowerElectronicConverter.count)  #Generate a name for the instance
        self.model_type = model_type
        self.signal_type = signal_type
        
        if self.model_type is 'switching':
            assert self.signal_type is 'square_wave' or self.signal_type is 'sinePWM', 'Switching model needs square or sine PWM as switching signal!'
        if self.model_type is 'average':
            assert self.signal_type is 'duty_cycle', 'Average model needs duty_cycle as switching signal!'
    
    def show_spec(self):
        """Print the specs."""
        
        print('Model type:{}'.format(self.model_type))
        print('Switching signal type:{}'.format(self.signal_type))
            
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
    
    def calc_average(self,m):
        """Calculate average voltage."""
        
        return Vdc
    
        
class PowerElectronicInverter(PowerElectronicConverter):
    """
    Inverter class.
    
    Attributes:
        ():
        
    """
    
    Rf = 0.01
    Lf = 2.0e-3
    
    Rload = 1.0
    
    def __init__(self,Vdc,model_type = 'average',signal_type='duty_cycle'):
        """Creates an instance of `Converter`.
        
        Args:
           Vdc (float): DC link voltage.
           
        Raises:
          ValueError: To be added.
        
        """
        
        super().__init__(model_type,signal_type)  #Initialize converter class (base class)
        
        self.update_Vdc(Vdc)
        
    
    @property                         #Decorator used for auto updating
    def y(self):
        """List of initial states"""
        
        return  [self.ia, 0.0]
    
    def update_Vdc(self,Vdc):
        """Update DC link voltage."""
        
        self.Vdc = Vdc
    
    def control_signal_calc(self,signals,t):
        """Calculate control signal."""
        
        if self.model_type is 'switching':
            signals = self.switching_signal_calc(t)
            control_signal = signals['switching']
        elif self.model_type is 'average':
            signals = self.average_signal_calc(signals,t)
            control_signal = signals['modulating']
        
        return control_signal
    
    def switching_signal_calc(self,signals,t):
        """Calculate switching signal."""
        
        if self.signal_type=='sinePWM':
            signal_dict = signals.sinosoidalPWM(t)
        elif self.signal_type=='square_wave':     
            signal_dict = signals.square_wave(t)
        #print(signal_dict)
        return signal_dict
    
    def average_signal_calc(self,signals,t):
        """Calculate switching signal."""
        
        if self.signal_type=='duty_cycle':
            signal_dict = signals.duty_cycle(t)  
        
        return signal_dict
        
    def half_bridge_switching(self,Vdc,S1):
        """Simulates a bridge in inverter"""
        
        self.update_Vdc(Vdc)

        S11 = self.calc_primary(S1)
        S12 = self.calc_complimentary(S1)
        
        assert S11+S12 == 1, 'S11 and S12 switches cannot be both ON or OFF at the same time in ideal half bridge.'

        #print('S11:{},S12:{}'.format(S11,S12))
        Van = (S11 - S12)*(self.Vdc/2)
        
        #print('Van:{}'.format(Van))

        return Van
    
    def half_bridge_average(self,Vdc,m):
        """Simulates a bridge in inverter"""
        
        self.update_Vdc(Vdc)

        assert m>=-1 and m <= 1, 'duty cycle should be between 0 and 1.'

        Van = m*(self.Vdc/2)
        
        #print('Van:{}'.format(Van))

        return Van
    
    def half_bridge_phasor(self,Vdc,m):
        """Simulates a bridge in inverter"""
        
        self.update_Vdc(Vdc)

        assert isinstance(m,complex), 'duty cycle should be complex phasor.'

        Van = m*(self.Vdc/2)
        
        #print('Van:{}'.format(Van))

        return Van
       
        #Current controller dynamics
    
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

    def get_model(self):
        """Select ODE model."""
        
        if self.model_type is 'switching' or self.model_type is 'average':
            model = self.ODE_model_EMT
        #elif self.converter.model_type is 'average':
        #    model = self.ODE_model_average
        elif self.model_type is 'dynamic_phasor':
            model = self.ODE_model_dynamicphasor
        
        return model
    
    def setup_model(self):
        """Initialize mode."""
        
        if self.model_type is 'switching' or self.model_type is 'average':
            self.ia = 0.0
        elif self.model_type is 'dynamic_phasor':
            self.iaR = 0.0
            self.iaI = 0.0        
    
    def vta_calc(self,Vdc,control_signal):
        """Calculate inverter terminal voltage."""
        
        if self.model_type is 'switching':
            vta = self.half_bridge_switching(Vdc,control_signal)
        elif self.model_type is 'average':
            vta = self.half_bridge_average(Vdc,control_signal)
       
        return vta
    
    
    def va_calc(self,t,grid,use_grid = True):
        """Calculate PCC voltage."""
        
        if use_grid:
            vpcc = grid.grid_voltage_calc(t)
        else:
            vpcc = self.Rload*self.ia
        
        return vpcc
    
    def ODE_model_switching(self,y,t):
        """ODE model of inverter branch."""
        
        self.ia,dummy = y   # unpack current values of y
        
        Vdc = 100.0 #Get DC link voltage
        switching_signal = self.control_signal_calc(t)
        
        self.vta = self.half_bridge_switching(Vdc,switching_signal)
        self.va = self.PCC_voltage_calc(self.ia,t)
        
        dia = (1/self.Lf)*(-self.Rf*self.ia -self.va + self.vta)
        
        result =  [dia,dummy]
        
        return np.array(result)
    
    def ODE_model_average(self,y,t):
        """ODE model of inverter branch."""
        
        self.ia,dummy = y   # unpack current values of y
               
        Vdc = 100.0 #Get DC link voltage
        modulating_signal = self.control_signal_calc(t)
                
        self.vta = self.half_bridge_average(Vdc,modulating_signal)
        self.va = self.PCC_voltage_calc(self.ia,t)
        
        dia = (1/self.Lf)*(-self.Rf*self.ia -self.va + self.vta)
                
        result =  [dia,dummy]
        
        return np.array(result)
    
    def ODE_model_EMT(self,y,t,signals,grid,sim):
        """ODE model of inverter branch."""
        
        self.ia,dummy = y   # unpack current values of y
        
        Vdc = 100.0 #Get DC link voltage
        control_signal = self.control_signal_calc(signals,t)
        self.vta = self.vta_calc(Vdc,control_signal)
        self.va = self.va_calc(t,grid,sim.use_grid)
                               
        dia = (1/self.Lf)*(-self.Rf*self.ia -self.va + self.vta)
                
        result =  [dia,dummy]
        
        return np.array(result)
    
    
    def ODE_model_dynamicphasor(self,y,t):
        """Dynamic phasor."""
        
        iaR,iaI = y   # unpack current values of y
        Vdc = 100.0 #Get DC link voltage
        winv = 2*math.pi*60
        
        self.ia = iaR + 1j*iaI
        self.vta = self.half_bridge_phasor(Vdc,1.0+1j*0.0)
        
        diaR = (1/self.Lf)*(-self.Rf*self.ia.real - self.Rload*self.ia.real + self.vta.real) + (winv)*self.ia.imag 
        diaI = (1/self.Lf)*(-self.Rf*self.ia.imag - self.Rload*self.ia.imag + self.vta.imag) - (winv)*self.ia.real  
        result =  [diaR,diaI]
        
        return np.array(result)    
    
    def power_calc(self,v,i):
        """Calcuate instantaneous power."""
        
        return v*i
    
    def show_states(self):
        """Show states."""
        
        print('Inverter states:{}'.format(self.y))