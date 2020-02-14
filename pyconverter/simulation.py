"""Class for simulation."""

import numpy as np
import math
import cmath
import scipy
import logging

from scipy import signal
from scipy.integrate import odeint,ode

import converter_utilities
import config
   
class Simulation:
    """
    Simulation class.
    
    Attributes:
        count (int): Number of converter objects.
        
    """
    
    count = 0 #Object count
    
    def __init__(self,signals,converter,grid=None):
        """Creates an instance of `Converter`.
        
        Args:
           fsw (float): Switching frequency in Hz.
           
        Raises:
          ValueError: If parameters corresponding to `Sinverter_rated` are not available.
        
        """
                
        Simulation.count = Simulation.count+1 #Increment count to keep track of number of converter model instances
        
        self.name = 'sim_'+str(Simulation.count)  #Generate a name for the instance
        self.converter = converter
        self.signals = signals        
        self.grid = grid       
        self.use_grid = False        
    
    @property
    def dt(self):
        """Calcuate time step."""
        
        if self.converter.model_type == 'switching':
            dt = 1/(10*self.signals.fsw) #Sampling rate in switching model simulation should be many times the switching frequency
        if self.converter.model_type == 'average':
            dt = 1/(self.signals.fsw) #Sampling rate in average model simulation can be same as the switching frequency
        if self.converter.model_type == 'dynamic_phasor':
            dt = 1/(100*self.signals.fm) #Sampling rate in average model simulation can be double the modulating frequency
        
        return dt
    
    def show_spec(self):
        """Print the specs."""
        
        self.signals.show_spec()
        self.converter.show_spec()
        print('Simulation time step:{:.5f}'.format(self.dt))    
        
    def simulate_inverter(self,tf):
        """Simulate an inverter.
        
        Args:
           tf (float): Simulation time in seconds.
        """
        
        #self.setup_simulation(tf)
        y0 = [0.0,0.0]
        self.t_t = np.arange(0,tf,self.dt)
        print('Number of timsteps:{}'.format(len(self.t_t)))
        
        solution,infodict = odeint(self.converter.ODE_model_EMT,y0,self.t_t,args=(self.signals,self.grid,self),full_output=1,printmessg=True,atol=1e-6,rtol=1e-6)
        #solution,infodict = odeint(self.ODE_model,self.y0,self.t_t,full_output=1,printmessg=True,atol=1e-6,rtol=1e-6)       
        
        if self.converter.model_type is 'switching' or self.converter.model_type is 'average':
            self.ia_t = solution[:,0]
            self.calc_timeseries()
            
            index= np.where(self.t_t==1/self.signals.fm)[0][0]
            self.iaRMS_t = converter_utilities.calc_window_RMS(self.ia_t, window_size=index)
            converter_utilities.plot_signal(self.t_t,self.ia_t,'ia',showPlot=False)
            converter_utilities.plot_signal(self.t_t,self.vta_t,'vta',showPlot=False)
            converter_utilities.plot_signal(self.t_t,self.va_t,'va',showPlot=False)
            converter_utilities.plot_signal(self.t_t[index-1:],self.iaRMS_t,'iaRMS',showPlot=True)
        elif self.converter.model_type is 'dynamic_phasor':
            self.iaR_t = solution[:,0]
            self.iaI_t = solution[:,1]
            self.ia_t = self.iaR_t+1j*self.iaI_t
            self.iaRMS_t = np.abs(self.ia_t)/math.sqrt(2)
            self.vaRMS_t = self.ia_t*self.Rload            
            converter_utilities.plot_signal(self.t_t,self.iaRMS_t,'iaRMS',showPlot=True)
        
        self.converter.show_states()
        print('iaRMS expected:{},iaRMS simulation:{}'.format((self.converter.Vdc/2)/math.sqrt(2),self.iaRMS_t[-1]))
        #converter_utilities.plot_FFT(1/self.dt,self.ia_t)
        
        return solution
    
    def setup_simulation(self,tf):
        """Setup the simulation."""
        
        self.t_t = np.arange(0,tf,self.dt)
        print('Number of timsteps:{}'.format(len(self.t_t)))
        self.y0 = [0.0,0.0]
        
        self.model = self.converter.get_model()  
        self.converter.setup_model()
    
    def ODE_model(self,y,t):
        """Select ODE model."""
        
        Vdc = 100.0 #Get DC link voltage
        control_signal = self.converter.control_signal_calc(self.signals,t)                
        
        self.converter.vta = self.converter.vta_calc(Vdc,control_signal)
        self.converter.va = self.converter.va_calc(t,self.grid,self.use_grid)
        
        result = self.model(y,t)
       
        return result
        
    
    def calc_timeseries(self):
        """Calculate time series."""
        
        self.vta_t = []
        self.va_t = []
        Vdc = 100.0 #Get DC link voltage
        
        for t,i in zip(self.t_t,self.ia_t):
            control_signal = self.converter.control_signal_calc(self.signals,t)                
            self.converter.ia = self.converter.vta_calc(Vdc,control_signal)
            self.converter.vta = self.converter.vta_calc(Vdc,control_signal)
            self.converter.va = self.converter.va_calc(t,self.grid,self.use_grid)
        
            self.vta_t.append(self.converter.vta)
            self.va_t.append(self.converter.va)
        
        self.vta_t = np.array(self.vta_t)
        self.va_t = np.array(self.va_t)
    