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
from models import InverterModels      

class PowerElectronicConverter:
    """
    Converter base class.
    
    Attributes:
        count (int): Number of converter objects.
        
    """
    
    count = 0 #Object count
    
    def __init__(self,model_type):
        """Creates an instance of `Converter`.
        
        Args:
           fsw (float): Switching frequency in Hz.
           
        Raises:
          ValueError: If parameters corresponding to `Sinverter_rated` are not available.
        
        """
                
        PowerElectronicConverter.count = PowerElectronicConverter.count+1 #Increment count to keep track of number of converter model instances
        
        self.name = 'converter_'+str(PowerElectronicConverter.count)  #Generate a name for the instance
        
        
        self.model_type = model_type
        
        """
        if self.model_type is 'switching':
            assert self.signal_type is 'square_wave' or self.signal_type is 'sinePWM', 'Switching model needs square or sine PWM as switching signal!'
        if self.model_type is 'average':
            assert self.signal_type is 'duty_cycle', 'Average model needs duty_cycle as switching signal!'
        """
        
    def check_model_type(self,model_type):
        """Check if model type is valid."""
        
        assert model_type in self.model_types, f'{model_type} is not a valid model type!'
    
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
    
        #Current controller dynamics
    
    
class PowerElectronicInverter(PowerElectronicConverter,InverterModels):
    """
    Inverter class.
    
    Attributes:
        ():
        
    """
    
    Rf = 0.01
    Lf = 1.0e-3
    
    Rload = 1.0
    
    inverter_types = ['single_phase_half_bridge','single_phase_full_bridge',
                     'three_phase_full_bridge']
    model_types = ['EMT_switching','EMT_average','dynamic_phasor']
    
    def __init__(self,Vdc,model_type = 'EMT_average',inverter_type='single_phase_half_bridge'):
        """Creates an instance of `Converter`.
        
        Args:
           Vdc (float): DC link voltage.
           
        Raises:
          ValueError: To be added.
        
        """
                
        self.check_model_type(model_type)
        
        super().__init__(model_type)  #Initialize converter class (base class)
        
        self.update_Vdc(Vdc)
        self.inverter_type =inverter_type
        
    
    @property                         #Decorator used for auto updating
    def y(self):
        """List of initial states"""
        
        return  [self.ia, 0.0]
    
    def update_Vdc(self,Vdc):
        """Update DC link voltage."""
        
        self.Vdc = Vdc
    """
    def control_signal_calc(self,signals,t):
        Calculate control signal.
        
        if self.model_type is 'EMT_switching':
            signals = self.switching_signal_calc(signals,t)
            control_signal = signals['switching']
        elif self.model_type is 'EMT_average':
            signals = self.average_signal_calc(signals,t)
            control_signal = signals['modulating']
        elif self.model_type is 'dynamicphasor':
            pass
        
        return control_signal
    """    
    def setup_model(self):
        """Initialize mode."""
        
        self.initialize_model()  
        self.vt_calc = self.select_vt_model()
        self.vpcc_calc = self.select_vpcc_model()
        self.ODE_model = self.select_ODE_model()
        #self.control_signal_calc = self.select_control_signal()
       
    
    def select_control_signal(self):
        """Select the control signal suitable for the problem."""
        
        if self.model_type is 'EMT_switching':
            if self.inverter_type == 'single_phase_half_bridge':
                control_signal = self.switching_signal_single_phase
            elif self.inverter_type == 'single_phase_full_bridge':
                raise NotImplementedError(f'{self.inverter_type} is not implemented!')
        
        elif self.model_type is 'EMT_average':
            if self.inverter_type == 'single_phase_half_bridge':
                control_signal = self.modulating_signal_single_phase
            elif self.inverter_type == 'single_phase_full_bridge':
                raise NotImplementedError(f'{self.inverter_type} is not implemented!')
        
        elif self.model_type is 'dynamic_phasor':
            if self.inverter_type == 'single_phase_half_bridge':
                control_signal = self.phasor_signal_single_phase
            elif self.inverter_type == 'single_phase_full_bridge':
                raise NotImplementedError(f'{self.inverter_type} is not implemented!')
        
        return control_signal
    
    def select_vt_model(self):
        """Get the terminal voltage model."""
        
        if self.model_type == 'EMT_switching':
            if self.inverter_type == 'single_phase_half_bridge':
                vt_model = self.single_phase_half_bridge_switching
            elif self.inverter_type == 'single_phase_full_bridge':
                vt_model = self.single_phase_full_bridge_switching
            elif self.inverter_type == 'three_phase_full_bridge':
                vt_model = self.three_phase_full_bridge_switching
            else:
                print(f'{self.inverter_type} not found for model type {self.model_type}!')
        elif self.model_type == 'EMT_average':
            if self.inverter_type == 'single_phase_half_bridge':
                vt_model = self.single_phase_half_bridge_average
            elif self.inverter_type == 'single_phase_full_bridge':
                vt_model = self.single_phase_full_bridge_average
            elif self.inverter_type == 'three_phase_full_bridge':
                raise NotImplementedError(f'{self.inverter_type} is not implemented!')
            else:
                print(f'{self.inverter_type} not found for model type {self.model_type}!')
        elif self.model_type == 'dynamicphasor':
            if self.inverter_type == 'single_phase_half_bridge':
                vt_model = self.single_phase_half_bridge_phasor
            elif self.inverter_type == 'single_phase_full_bridge':
                vt_model = self.single_phase_full_bridge_phasor
            elif self.inverter_type == 'three_phase_full_bridge':
                raise NotImplementedError(f'{self.inverter_type} is not implemented!')
            else:
                print(f'{self.inverter_type} not found for model type {self.model_type}!')    
        
        print(type(vt_model))
        return vt_model

    def select_vpcc_model(self,grid=None):
        """Get the PCC voltage model."""     
        
        if not grid:
            vpcc_model = self.v_load_model()
         
        return vpcc_model

    def select_ODE_model(self):
        """Select ODE model."""
        
        if self.model_type is 'EMT_switching' or self.model_type is 'EMT_average':
           if self.inverter_type  is 'single_phase_half_bridge' or self.inverter_type  is 'single_phase_full_bridge':
              ODE_model = self.ODE_model_single_phase_EMT
           elif self.inverter_type  is 'three_phase_full_bridge':
              raise NotImplementedError(f'{self.inverter_type} is not implemented!')
        elif self.model_type is 'dynamic_phasor':
           if self.inverter_type is 'single_phase_half_bridge' or self.inverter_type  is 'single_phase_full_bridge':
              ODE_model = self.ODE_model_single_phase_dynamicphasor
           elif self.inverter_type  is 'three_phase_full_bridge':
              raise NotImplementedError(f'{self.inverter_type} is not implemented!')
            
        return ODE_model
    
    def initialize_model(self):
        """Initialize mode."""
        
        if self.model_type is 'EMT_switching' or self.model_type is 'EMT_average':
            if self.inverter_type  is 'single_phase_half_bridge' or self.inverter_type  is 'single_phase_full_bridge':
                self.ia = 0.0
            elif self.inverter_type  is 'three_phase_full_bridge':
                raise NotImplementedError
           
        elif self.model_type is 'dynamic_phasor':
            if self.inverter_type  is 'single_phase_half_bridge' or self.inverter_type  is 'single_phase_full_bridge':
                self.iaR = 0.0
                self.iaI = 0.0
            if self.inverter_type  is 'three_phase_full_bridge':
                raise NotImplementedError           
    """
    def vta_calc(self,Vdc,control_signal):
        Calculate inverter terminal voltage.
        
        if self.model_type is 'switching':
            vta = self.half_bridge_switching(Vdc,control_signal)
        elif self.model_type is 'average':
            vta = self.half_bridge_average(Vdc,control_signal)
       
        return vta
    """
    
    def v_load_model(self):
        """Calculate voltage across load at PCC."""        
        
        return self.Rload*self.ia
    
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
    
    
    def power_calc(self,v,i):
        """Calcuate instantaneous power."""
        
        return v*i
    
    def show_states(self):
        """Show states."""
        
        print('Inverter states:{}'.format(self.y))