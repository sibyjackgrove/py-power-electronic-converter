import matplotlib.pyplot as plt
from scipy import fftpack
import numpy as np

import os
  

def calc_RMS(signals_array):
    """Calculate RMS from time series."""
        
    return np.sqrt(np.mean(np.square(signals_array)))

def calc_window_RMS(a, window_size=400):
    """Calculate RMS from time series."""    
    
    return np.sqrt(sum([a[window_size-i-1:len(a)-i]**2 for i in range(window_size-1)])/window_size)        
    
def plot_signal(time,signal,label,showPlot=False):
    """Plot signal"""
    
    plt.plot(time, signal,label=label)
    plt.legend(fontsize =18)

def show_plot_signal():
    """show plots."""
    
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')            
    plt.show()    
        
def plot_FFT(f_s,y):
    """Plot FFT."""
    
    Y = fftpack.fft(y)
    freqs = fftpack.fftfreq(len(y)) * f_s
    fig, ax = plt.subplots()
    
    ax.stem(freqs, np.abs(Y),use_line_collection=True)
    ax.set_xlabel('Frequency in Hertz [Hz]')
    ax.set_ylabel('Frequency Domain (Spectrum) Magnitude')
    plt.show()
