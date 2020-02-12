import matplotlib.pyplot as plt
import os
  
def plot_signal(time,signal,label,showPlot=False):
    """Plot signal"""
    
    plt.plot(time, signal,label=label)
    plt.legend(fontsize =18)
        
    if showPlot:
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')            
        plt.show()    