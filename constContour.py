from glob import glob
import nanonispy as nap
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
import numpy as np
from matplotlib.font_manager import FontProperties

import commonFunctions as cf

def findRamp():
    
    datFiles = glob('data/2023-05-19/IV_00*.dat')
    
    # setup figure to hold the IV, dI/dV, LDOS plots
    fig1 = plt.figure(1, figsize=(7, 8))
    ax_IV, ax_zV, ax_dzdV = fig1.subplots(3, 1, sharex=True)
    
    for datF in datFiles:
        
        specObj  = nap.read.Spec(datF)
        
        biasDataFwd = cf.extract_channel(specObj, 'Bias (V)')
        biasDataBwd = cf.extract_channel(specObj, 'Bias [bwd] (V)')
        zDataFwd = cf.extract_channel(specObj, 'Z (m)')
        zDataBwd = cf.extract_channel(specObj, 'Z [bwd] (m)')
        currentDataFwd = cf.extract_channel(specObj, 'Current (A)')
        currentDataBwd = cf.extract_channel(specObj, 'Current [bwd] (A)')
        
        # plot I-V spectrum 
        ax_IV.plot(biasDataFwd, 1e12*currentDataFwd)
        ax_IV.set_ylabel('Current (pA)')
        
        # plot z-V spectrum
        ax_zV.plot(biasDataFwd, 1e10*zDataFwd)
        ax_zV.set_ylabel('z height (Å)')
        
        # plot dz/dV as function of V
        # first find dz/dV at each point, normalise here as constant dV.
        dz_dV = np.gradient(zDataFwd, biasDataFwd[1]-biasDataFwd[0])
        # then find average (ramp) value
        meanHeight = np.mean(dz_dV)
        # meanHeight = np.mean(dz_dV[25:])
        # finally plot
        ax_dzdV.plot(biasDataFwd, 1e12*dz_dV)
        ax_dzdV.set_xlabel('Bias (V)')
        ax_dzdV.set_ylabel('dz/dV (pm/V)')
        ax_dzdV.hlines(1e12*meanHeight, biasDataFwd[0], biasDataFwd[-1], 'k', '--', label=f'{round(1e12*meanHeight)} pm/V')
        # ax_dzdV.hlines(1e12*meanHeight, biasDataFwd[25], biasDataFwd[-1], 'k', '--', label=f'{round(1e12*meanHeight)} pm/V')
        ax_dzdV.legend()
        
    plt.savefig('figures/dzdV_constContour.png')
    plt.close()
        
    return


def main():
   findRamp() 

if __name__=='__main__':
 	main()