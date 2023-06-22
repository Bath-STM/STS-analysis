#########################################################
# Script to plot the constant current contour. 
# Used with I-V STS data taken with the feedback loop ON. 
########################################################
from email.policy import default
from os import sep
import sys
from glob import glob

import click
import matplotlib.pyplot as plt
import nanonispy as nap
import numpy as np
from matplotlib.font_manager import FontProperties
from scipy.optimize import curve_fit

import commonFunctions as cf


def findRamp(infile, separate, outStub):
    
    # datFiles = glob('data/2023-06-20/IV_*.dat')
    # datFiles = ['IV_003.dat', 'IV_004.dat', 'IV_005.dat']
    print(infile)
    sys.exit()
    
    # setup figure to hold the IV, dI/dV, LDOS plots
    fig1 = plt.figure(1, figsize=(7, 8))
    ax_IV, ax_zV, ax_dzdV = fig1.subplots(3, 1, sharex=True)
    
    print(datFiles)
    for datF in datFiles:
        
        expStub = datF.split('_')[-1]
        expNum = expStub.split('.')[0]
        
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
        
        if separate:
            plt.savefig(f'figures/dzdV_constContour_{outStub}_{expNum}.png')
            ax_IV.cla()
            ax_zV.cla()
            ax_dzdV.cla()
        
    if not separate:
        plt.savefig(f'figures/dzdV_constContour_{outStub}_{expNum}.png')
    
    plt.close()
        
    return

@click.command()
@click.option('--infile', '-i', help="Input .dat file(s)", type=list, prompt=True, multiple=True)
@click.option('--separate', '-s', help="Plot data on separate plots or same figure", is_flag=True, default=False)
@click.option('--outStub', '-o', help="Out file name stub", type=str, prompt=True, multiple=False)


def main(infile, separate, outStub):
   findRamp(infile, separate, outStub) 

if __name__=='__main__':
 	main()