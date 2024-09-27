#########################################################
# Plot I(z) aquisition(s) & the associated derived quantities.
# Plots average current, conductance, kappa, the effective work function (phi_eff) and dI/dz 
# z axis values may be offset by a calibration value  given in pm (defualt 0 pm)
# ------------------------------------------------------- 
# Example command `python3 izPlotter.py -iz data/2024-03-27/multi2_Inj_240327_\[2-9]_\*.txt -iz data/2024-03-27/multi2_Inj_240327_1[0-9]_\*.txt -o 240327_2-19`
# -------------------------------------------------------
# TODO - select which areas to fit in and implement these fits 
#########################################################

import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import sys
import numpy as np
from matplotlib.font_manager import FontProperties
import click

import commonFunctions as cf

def exp_function(z, I_0, kappa):
    return I_0*np.exp(-2*kappa*z)

def gaussian(x, *popt):
    A, mu, sigma = popt
    return (A*np.exp(-(x-mu)**2/(2.*sigma**2)))

def offset_gaussian(x, *popt):
    A, mu, sigma, offset = popt
    return (A*np.exp(-(x-mu)**2/(2.*sigma**2))) + offset

def line(x, *popt):
    m, c = popt
    return m*x + c

def izPlot(infilez, calibration, outstub):
    
    # Define a few constants
    hbar = 6.626e-34 / (2*np.pi)
    m_e = 9.109e-31
    e = 1.602e-19
    
    if not infilez:
        print('\n#####\nNo infile paths input. Exiting.\n#####\n')
        sys.exit(1)
    
    try:
        inFilesZ = cf.globFiles(infilez)
    except FileNotFoundError as err:
        print(err)
    
    if not inFilesZ:
        print('\n#####\nNo valid infile paths specified. Exiting.\n#####\n')
        sys.exit(1)
    
    outerHeight = [] # hold the relative tip height for each point in each experiment
    outerCurrent = [] # hold the current for each point in each experiment
    outerCond = [] # hold the conductance (in G_0) for each point in each experiment
    outerLIX = [] # hold the dIdz for each point in each experiment
    outerKappa = [] # hold the kappa for each point in each experiment
            
    for inFileZ in inFilesZ:
        try:
            dfZ = cf.extract_file(inFileZ)
        except cf.UnhandledFileError as err:
            print(err)
            continue
        except ValueError as err:
            print(err)
            continue
        
        # Check got a sufficient spectrum, not tip withdrawn
        # remove NaNs from the df and check have sufficient data points left
        if len(dfZ['current'][~np.isnan(dfZ['current'])]) < 30:
            print(f'\n### Not plotted {inFileZ}. Insufficient data points.\n')
            continue
        
        # Calibrate z data
        dfZ['height'] = dfZ['height'] + 1e-12*calibration
        ########### changed for MPhys data
        # dfZ['lIX'] = -dfZ['lIX']
        
        # Mask data to non-noisy region
        quiet = abs(dfZ['current']) > 8e-12
        
        # get the info for the averaged plots
        outerHeight.extend(dfZ[quiet]['height'].tolist())
        outerCurrent.extend(dfZ[quiet]['current'].to_list())
        outerLIX.extend(dfZ[quiet]['lIX'].to_list())
        
        # Plot the current for each exp
        plt.figure(1)
        plt.plot(1e12*dfZ[quiet]['height'], 1e12*dfZ[quiet]['current'], 'k', alpha=0.25, label='Measured')
        
        # and the conductance
        # convert from Siemens (1/V * I) to G_0 -> 1 G_0 = 7.75x10^-5 S
        cond = (1/dfZ[quiet]['bias'])*dfZ[quiet]['current']*7.75e-5
        plt.figure(4)
        plt.plot(1e12*dfZ[quiet]['height'], cond, 'k', alpha=0.25, label='Measured')
        # add the conductance for the averaged plots
        outerCond.extend(cond)
        
        # Plot kappa for each sweep 
        kappa = -dfZ[quiet]['lIX']/(2*dfZ[quiet]['current'])
        # Change sign of kappa depending on lockin phase sign
        # kappa = dfZ[quiet]['lIX']/(2*dfZ[quiet]['current'])
        # if data from a .dat file need to divide by dz (usually 10 pm)
        # dz = 10e-12
        # kappa = (-dfZ[quiet]['lIX']/dz)/(2*dfZ[quiet]['current'])
        plt.figure(2)
        plt.plot(1e12*dfZ[quiet]['height'], 1e-9*kappa, 'k', alpha=0.25, label='Measured')
        # add the kappa for the averaged plots
        outerKappa.extend(kappa)
        
        # and turn kappa into phi_eff (in eV)
        # kappa = sqrt(2 m_e phi)/ hbar
        # (kappa * hbar)^2 / 2 m_e
        phi_eff = (kappa * hbar)**2 / (2*m_e)
        plt.figure(3)                       # convert from J to eV
        plt.plot(1e12*dfZ[quiet]['height'], phi_eff/e, 'k', alpha=0.25, label='Measured')
        
    # Collect data for averaged plots
    # Check the data exists
    if not outerCurrent:
            print(f'\n#####\nNo data to plot. Exiting.\n#####\n')
            sys.exit(1)
    
    # Then make list of lists into a df
    outerdf = pd.DataFrame(np.column_stack([outerHeight, outerCurrent, outerCond, outerLIX, outerKappa]), columns=['height', 'current', 'conductance', 'lIX', 'kappa'])
    # print(outerdf.to_string())
    
    # Find average current
    plt.figure(99)
    sumCurrent, heightEdges, patches = plt.hist(outerdf['height'], weights=outerdf['current'], bins=70)
    binCount, heightEdges, patches = plt.hist(outerdf['height'], bins=70)
    plt.close()
    # can now normalise/find mean current in each bin
    normCurrent = sumCurrent / binCount 
    
    # and average conductance
    plt.figure(99)
    sumCond, heightEdges, patches =plt.hist(outerdf['height'], weights=outerdf['conductance'], bins=70)
    plt.close()
    # can now normalise/find mean conductance in each bin. Bin counts the same as above
    normCond = sumCond / binCount
    
    # Repeat for dIdz (lIX)
    plt.figure(99)
    sumLIX, heightEdges, patches = plt.hist(outerdf['height'], weights=outerdf['lIX'], bins=70)
    plt.close()
    # can now normalise/find mean dIdz in each bin. Bin counts the same as above
    normLIX = sumLIX / binCount 
    
    # Finally for kappa
    plt.figure(99)
    sumKappa, heightEdges, patches = plt.hist(outerdf['height'], weights=outerdf['kappa'], bins=70)
    plt.close()
    # can now normalise/find mean kappa in each bin, again using same bin counts
    normKappa = sumKappa / binCount  
    
    avedf = pd.DataFrame(np.column_stack([heightEdges[:-1], normCurrent, normCond, normLIX, normKappa]), columns=['binnedHeight', 'aveCurrent', 'aveCond', 'aveLIX', 'aveKappa'])
    
    ### Mask data to exp current/flat kappa region
    
    ### Fits in flat region
    
    # Make the averaged plots
    # starting with current
    plt.figure(1)
    plt.plot(1e12*avedf['binnedHeight'], 1e12*avedf['aveCurrent'], 'r', label='Average')
    # ...and the exp fit found in the flat region
    # plt.plot(1e12*avedf['binnedHeight'], 1e12*projCurrent, '--c', label='Fit in flat region')
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel('I (pA)', fontsize=20)
    plt.xlim(140, 890)
    # plt.tight_layout()
    plt.yscale('log')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('small')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper right', prop=fontP)
    plt.savefig(f'figures/ave_current_fit_{outstub}.png', bbox_inches='tight')
    plt.close()
    
    # ... and the averaged conductance
    # convert from Siemens (I/V) to G_0 -> 1 G_0 = 7.75x10^-5 S, V = 1 V
    plt.figure(4)
    plt.plot(1e12*avedf['binnedHeight'], avedf['aveCond'], 'r', label='Average')
    # ...and the exp fit found in the flat region
    # plt.plot(1e12*avedf['binnedHeight'], 1e12*projCurrent*7.75e-5, '--c', label='Fit in flat region')
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel(r'Conductance ($G_0 = \frac{2e^2}{h}$)', fontsize=20)
    plt.xlim(140, 890)
    plt.tight_layout()
    plt.yscale('log')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('small')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper right', prop=fontP)
    plt.savefig(f'figures/ave_cond_fit_{outstub}.png', bbox_inches='tight')
    plt.close()
    
    # now the averaged kappa
    plt.figure(2)
    plt.plot(1e12*avedf['binnedHeight'], 1e-9*avedf['aveKappa'], '--r', label='Average')
    # plt.plot(projRange, projKappa, '--y', label='Fit')
    # plt.plot(binnedHeight, projKappa, '--y', label='Fit')
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel(r'$\kappa$ (nm$^{-1}$)', fontsize=20)
    plt.xlim(140, 890)
    # plt.ylim(top=13.5)
    plt.ylim(0, 13)
    # plt.tight_layout()
    # plt.yscale('log')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('small')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='lower right', prop=fontP)
    plt.savefig(f'figures/ave_kappa_{outstub}.png', bbox_inches='tight')
    plt.close()
    
    # ... and the associated phi_eff
    phi_eff_ave = (avedf['aveKappa'] * hbar)**2 / (2*m_e) # in Joules
    plt.figure(3)                       # convert from J to eV
    plt.plot(1e12*avedf['binnedHeight'], phi_eff_ave/e, '--r', label='Average')
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel(r'$\phi_{eff}$ (eV)', fontsize=20)
    plt.xlim(140, 890)
    plt.ylim(0, 6.5)
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    # plt.tight_layout()
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('small')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='lower right', prop=fontP)
    plt.savefig(f'figures/ave_phiEff_{outstub}.png', bbox_inches='tight')
    plt.close()
    
    # separating the dIdz 
    plt.figure(5)
    plt.plot(1e12*avedf['binnedHeight'], avedf['aveLIX'], '--r', label='Average')
    # dx = avedf['binnedHeight'][1] - avedf['binnedHeight'][0]
    numDiffCurrent = np.gradient(avedf['aveCurrent'])
    # print(numDiffCurrent)
    plt.plot(1e12*avedf['binnedHeight'], 1e10*numDiffCurrent, 'c', label='Num diff')
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel(r'$\partial I / \partial z$', fontsize=20)
    plt.xlim(140, 890)
    # plt.tight_layout()
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('small')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='lower right', prop=fontP)
    plt.savefig(f'figures/ave_dIdz_{outstub}.png', bbox_inches='tight')
    plt.close()
    
    return

@click.command()
@click.option('--infilez', '-iz', help="Input file(s) with z lock-in modulated", type=str, prompt=True, multiple=True)
@click.option('--calibration', '-c', help="Calibration height to be added to the raw values, given in pm", type=float, multiple=False, default=0)
@click.option('--outstub', '-o', help="Out file name stub", type=str, prompt=True, multiple=False)

def main(infilez, calibration, outstub):
    
    izPlot(infilez, calibration, outstub)
    
if __name__=='__main__':
 	main()