"""Commonly used calculations on electrical quantities."""

from __future__ import division
import numpy as np
import math
import cmath
import sys
import time
import six
import scipy.io as sio
#from numba import jit



def Uunbalance_calc(ua,ub,uc):
    """Calculate voltage/current unbalance."""
        
    uavg = (ua + ub + uc)/3
        
    return (max(ua,ub,uc) - min(ua,ub,uc))/uavg



#@jit(nopython=True)
def Urms_calc(ua,ub,uc):
    """Function to calculate rms value of scalar phasor quantities."""
    
    return math.sqrt((pow(abs(ua),2)+pow(abs(ub),2)+pow(abs(uc),2))/3.0)/math.sqrt(2)  #Pure python implementation is faster      

def Urms_calc_1phase(ua):
    """Function to calculate rms value of scalar phasor quantities for single phase."""
    
    return abs(ua)/math.sqrt(2)  #Pure python implementation is faster      

#@jit(nopython=True)    
def Ub_calc(Ua):
    """Convert phase A quantity to Phase B."""
    
    return Ua*pow(math.e,1j*(-(2/3)*math.pi))  #Shift by -120 degrees
#@jit(nopython=True)    
def Uc_calc(Ua):
    """Convert phase A quantity to Phase C."""
    
    return Ua*pow(math.e,1j*((2/3)*math.pi))  #Shift by -120 degrees

def relative_phase_calc(Uph1,Uph2,DEGREES=False):
    """Calculate relative phase between phasors between 0 to 2pi or 0 to 360 degrees."""
    
    if DEGREES:
        del_phase = math.degrees(cmath.phase(Uph1)-cmath.phase(Uph2)) % 360
        
    else:
        del_phase = math.radians(math.degrees(cmath.phase(Uph1)-cmath.phase(Uph2)) % 360)
        
    return del_phase
    

#@jit(nopython=True)    
def phasor_to_time(upha = 1+1j*0.0,uphb = -0.5-1j*0.867,uphc = -0.5+1j*0.867,w=2.0*math.pi*60.0,t=0.0):
    """Convert a,b,c quantities from phasor domain to time domain."""
    
    ra,pha = cmath.polar(upha)
    rb,phb = cmath.polar(uphb)
    rc,phc = cmath.polar(uphc)
    #ua = (ra*np.exp(1j*(w*t+pha-(math.pi/2)))).real
    #ub = (rb*np.exp(1j*(w*t+phb-(math.pi/2)))).real
    #uc = (rc*np.exp(1j*(w*t+phc-(math.pi/2)))).real
    
    ua = ra*pow(math.e,1j*(w*t+pha-(math.pi/2))).real
    ub = rb*pow(math.e,1j*(w*t+phb-(math.pi/2))).real
    uc = rc*pow(math.e,1j*(w*t+phc-(math.pi/2))).real
    return ua,ub,uc   

def phasor_to_time_1phase(uph,w,t):
    """Convert a,b,c quantities (time series) from phasor domain to time domain."""
    
    r,ph = cmath.polar(uph)
    
    return r*pow(math.e,1j*(w*t+ph-(math.pi/2))).real
#@jit(nopython=True)
def abc_to_dq0(ua,ub,uc,wt=2*math.pi):
    """Convert to d-q."""
    Us = (2/3)*(ua + ub*pow(math.e,1j*((2/3)*math.pi)) + uc*pow(math.e,1j*(-(2/3)*math.pi)))*pow(math.e,1j*(-wt))
        
    ud = Us.real
    uq = Us.imag
    u0 = (1/3)*(ua+ub+uc)
    return ud,uq,u0

def dq0_to_abc(ud,uq,u0,wt=2*math.pi):
    """Convert to abc."""
    
    ua = ud*math.cos(wt) - uq*math.sin(wt) + u0
    ub = ud*math.cos(wt-(2/3)*math.pi) - uq*math.sin(wt-(2/3)*math.pi) + u0
    uc = ud*math.cos(wt+(2/3)*math.pi) - uq*math.sin(wt+(2/3)*math.pi) + u0

    return ua,ub,uc

def alpha_beta_to_d_q(ualpha,ubeta,wt):
    """Convert alpha-beta to d-q."""
    
    Us = (ualpha + 1j*ubeta)*pow(math.e,-1j*(wt))
    #Us = (ualpha + 1j*ubeta)*pow(math.e,-1j*(wt-(math.pi/2)))
    
    ud = Us.real
    uq = Us.imag
    
    #print(ud,uq)
    #ud = ualpha*math.sin(wt) - ubeta*math.cos(wt)
    #uq = ualpha*math.cos(wt) + ubeta*math.sin(wt)
    
    #ud = ualpha*math.cos(wt) + ubeta*math.sin(wt)
    #uq = -ualpha*math.sin(wt) + ubeta*math.cos(wt)
    
    return ud,uq

def phasor_to_symmetrical(upha,uphb,uphc):
    """Convert to zero sequence."""     
    a = pow(math.e,1j*((2/3)*math.pi))
    aa = pow(math.e,1j*((4/3)*math.pi))
    u0 = (1/3)*(upha + uphb + uphc)
    u1 = (1/3)*(upha + a*uphb + (aa)*uphc)  #Positive sequence
    u2 = (1/3)*(upha + (aa)*uphb + a*uphc) #Negative sequence
    
    return u0,u1,u2

def phasor_to_zero_sequence(upha,uphb,uphc):
    """Convert to zero sequence.""" 
    
    u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
    
    return u0,u0,u0	

def phasor_to_positive_sequence(upha,uphb,uphc):
    """Convert to positive sequence.""" 
    
    u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
    
    return u1,u1*pow(math.e,1j*((4/3)*math.pi)),u1*pow(math.e,1j*((2/3)*math.pi))

def phasor_to_negative_sequence(upha,uphb,uphc):
    """Convert to negative sequence.""" 
    
    u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
    
    return u2,u2*pow(math.e,1j*((4/3)*math.pi)),u2*pow(math.e,1j*((2/3)*math.pi))    
