from curses.panel import bottom_panel
from glob import glob
import nanonispy as nap
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import scipy
from scipy.optimize import curve_fit
import sys
import numpy as np
import cycler
from matplotlib.font_manager import FontProperties

import commonFunctions as cf

def findRamp():
    
    datFiles = glob('data/2023-05-10/IV_modA_0mV_FC_001.dat')
    
    for datF in datFiles:
        
        specObj  = nap.read.Spec(datF)
        
        print(type(specObj))
        
        biasDataFwd = cf.extract_channel(specObj, 'Bias (V)')
        biasDataBwd = cf.extract_channel(specObj, 'Bias [bwd] (V)')
        zDataFwd = cf.extract_channel(specObj, 'Z (m)')
        zDataBwd = cf.extract_channel(specObj, 'Z [bwd] (m)')
        
        dz_Fwd = np.diff(10e12*zDataFwd)
        dV_Fwd = np.diff(biasDataFwd)
        
        
        # print(dz_Fwd)
        
        plt.figure(1)
        plt.plot(biasDataFwd[1:], dz_Fwd/dV_Fwd)
        # plt.yscale('log')
        # plt.ylim(bottom=10e-2)
        gradients = np.gradient(10e12*zDataFwd)#, biasDataFwd[1]-biasDataFwd[0])
        dz_dV = gradients / dV_Fwd[0]
        print(gradients)
        print(dz_dV)
        # plt.hlines(np.mean(np.gradient(zDataFwd), biasDataFwd[1]-biasDataFwd[0]), biasDataFwd, biasDataFwd, 'k', '--')
        # plt.hlines(np.mean(dz_Fwd[25:]/dV_Fwd[25:]), biasDataFwd[26], biasDataFwd[-1], 'k', '--')
        plt.savefig('figures/dzdV_constContour.png')
        plt.close()
        
    
    
    return


def main():
   findRamp() 

if __name__=='__main__':
 	main()