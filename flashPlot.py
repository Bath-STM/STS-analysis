#########################################################
# Explainations
# ------------------------------------------------------- 
# Example command ``
# -------------------------------------------------------
# TODO - 
#########################################################
from glob import glob
from matplotlib import axes
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import sys
import numpy as np
import cycler
from matplotlib.font_manager import FontProperties

import datetime
from matplotlib.dates import DateFormatter
#  a = datetime.timedelta(seconds=65)
# datetime.timedelta(0, 65)
# >>> str(a)


# import commonFunctions as cf

def convertSeconds(secondsArray):
    outArray = []
    for second in secondsArray:
        delta = datetime.timedelta(seconds=second)
        outArray.append(str(delta))
    return outArray

def plotFlash():
    
    flashFiles = glob('data/flash-22-05-23.txt')
    print(flashFiles)
    
    for flashF in flashFiles:
        print(flashF)
        
        time, pressure, current, temp = np.loadtxt(flashF, usecols=(0, 1, 3, 5), skiprows=2, unpack=True, delimiter='\t', dtype=float)
        
        
        hms = convertSeconds(time)
        # print(hms)
        print(hms,  pressure)
        
        fig1 = plt.figure(1, figsize=(5.75, 8))
        ax_current, ax_pressure, ax_temp = fig1.subplots(3, 1, sharex=True)
        
        # ax_current.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
        # ax_pressure.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
        # ax_temp.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
        
        ax_current.plot_date(hms, current)
        ax_current.set_ylabel('Current (A)')
        ax_current.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
        
        # plt.show()
        # sys.exit()
        
        ax_pressure.plot_date(hms, 1e12*pressure)
        ax_pressure.set_ylabel(r'Pressure (x10$^{-9}$ mbar)')
        ax_pressure.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
        
        plt.show()
        sys.exit()
        
        ax_temp.plot_date(hms, temp)
        ax_temp.set_ylabel('Temperature °C')
        ax_temp.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
        
        plt.figure(1)
        # hh_mm = DateFormatter('%H:%M')
        
        plt.xlabel('Time (h:mm:ss)')
        plt.savefig(f'figures/flashing.png')
        plt.savefig(f'figures/flashing.pdf')
        plt.close()
        
    return

def main():
    
    plotFlash()
    

if __name__=='__main__':
 	main()