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


class Grid:
    """
    Grid class.
    
    Attributes:
        count (int): Number of grid objects.
        
    """
    
    count = 0 #Object count
    
    def __init__(self,Agrid,fgrid):
        """Creates an instance of `SwitchingSignals`.
        
        Args:
           fsw (float): Switching frequency in Hz.
           fm (float): Fundamental frequency of output waveform in Hz.
           
        Raises:
          ValueError: To be added
        
        """
        
        self.update_grid_waveform(Agrid,fgrid)
           
    
    @property                         #Decorator used for auto updating
    def mf(self,):
        """Frequency modulation index"""
        
        return self.fsw/self.fm
    
    def update_grid_waveform(self,Agrid=1.0,fgrid = 50.0):
        """Update grid voltage amplitude and frequency."""
        
        self.Agrid=Agrid
        self.fgrid = fgrid
        self.Tgrid= 1/self.fgrid
    
    def grid_voltage_calc(self,t):
        """Create a sinusoid time series."""
        
        return self.Agrid*np.sin(2*np.pi*self.fgrid*t)    
        
        