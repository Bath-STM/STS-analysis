#########################################################
# Explainations
# ------------------------------------------------------- 
# Example command ``
# -------------------------------------------------------
# TODO - 
#########################################################
# from glob import glob
import nanonispy as nap
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import scipy
import sys
import numpy as np
import cycler
import click
from matplotlib.font_manager import FontProperties
from scipy import ndimage
from matplotlib.colors import LogNorm

import commonFunctions as cf

def kappaMap(infilez, outstub):
    
    if not infilez:
        return
    
    inFilesZ = cf.globFiles(infilez)
            
    outerHeight = [] # hold the relative tip height for each point in each experiment
    outerBias = [] # hold the bias for each point in each experiment
    outerKappa = [] # hold the kappa for each point in each experiment
    
    for inFileZ in inFilesZ:
        try:
            dfZ = cf.extract_file(inFileZ)
        except ValueError:
            print('\nPlease check the file path\n')
            
        # Check got a spectrum, not tip withdrawn
        # remove NaNs from the df and check have sufficient data points left
        if len(dfZ['current'][~np.isnan(dfZ['current'])]) < 30:
            print(f'\n### Not plotted {inFileZ}\n')
            continue
        
        # # Calibrate z data (0_raw = +1V, 100 pA) to Science paper (dark mol)
        # dfZ['height'] = dfZ['height'] + 700e-12
        # ... and then to the Gywdion measurements of heights 
        # # starting with a vacancy/noBromo dark
        # dfZ['height'] = dfZ['height'] + 720e-12
        # finally a bright adatom
        dfZ['height'] = dfZ['height'] + 750e-12
        
        # get the info for the kappa plot
        outerHeight.extend(dfZ['height'].tolist())
        outerBias.extend(dfZ['bias'].tolist())
        
        
        # These out by a factor? lIX -> delta I so need to / by delta z (x pm)....
        # ... this depends what is output in the multi_Inj_xxx.txt files and if 
        # this scaling by the amplitude is already done...?
        # Alredy taken care of in wvtxt files. Needs to be done for nanonis .dat files
        kappa = -dfZ['lIX']/(2*dfZ['current'])
        # delta z = 10pm for 03/08/23 data
        # kappa = (-dfZ['lIX']/10e-12)/(2*dfZ['current'])

        outerKappa.extend(kappa)
        
    # plot average kappa
    # first make list of lists into a df
    kappadf = pd.DataFrame(np.column_stack([outerHeight, outerBias, outerKappa]), columns=['height', 'bias', 'kappa'])
    # print(kappadf.to_string())
    
    # make two 2d hists. First sums the kappa values from each exp in the bin
    # second counts the exps in the bin
    plt.figure(1)
    # sumkappa, biasEdges, heightEdges, im = plt.hist2d(kappadf['bias'], 1e12*kappadf['height'], weights=1e9*kappadf['kappa'], bins=30)
    sumkappa, biasEdges, heightEdges, im = plt.hist2d(kappadf['bias'], 1e12*kappadf['height'], weights=1e-9*kappadf['kappa'], bins=30)
    binCount, biasEdges, heightEdges, im = plt.hist2d(kappadf['bias'], 1e12*kappadf['height'], bins=30)
    plt.close()
    # can now normalise/find mean kappa in each bin
    normkappa = sumkappa / binCount     
    
    # finally display the average kappa. 
    # Data needs to be transposed because of Python weirdness 
    plt.figure(2)
    plt.imshow(normkappa.T, origin='lower',  cmap='viridis', aspect='auto', extent=(biasEdges[0], biasEdges[-1], heightEdges[0], heightEdges[-1]))
    # plt.imshow(normkappa.T, origin='lower',  cmap='viridis', aspect='auto', vmin=3, vmax=12, extent=(biasEdges[0], biasEdges[-1], heightEdges[0], heightEdges[-1]))
    # plt.imshow(normkappa.T, origin='lower',  cmap='viridis', aspect='auto', extent=(biasEdges[0], biasEdges[-1], heightEdges[0], heightEdges[-1]), norm=LogNorm())
    # plt.xlim(0, 3)
    plt.ylabel('Relative tip height (pm)', fontsize=20)
    plt.xlabel('Bias (V)', fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    cbar = plt.colorbar(boundaries=(np.linspace(2, 12, 100)))
    cbar.set_label(r'$\kappa$ (nm$^{-1}$)')
    plt.savefig(f'figures/kappa2dMap_{outstub}.png')
    plt.close()
        
    return

def lDOSMap(infilev, outstub):
    
    if not infilev:
        return
    
    inFilesV = cf.globFiles(infilev)
            
    outerHeight = [] # hold the relative tip height for each point in each experiment
    outerBias = [] # hold the bias for each point in each experiment
    outerLDOS = [] # hold the LDOS for each point in each experiment
    outerDiffCond = [] # hold the differential conductance (dI/dV) for each point in each experiment
    outerCond = [] # hold the conductance (I/V) for each point in each experiment
    
    for inFileV in inFilesV:
        try:
            dfV = cf.extract_file(inFileV)
        except ValueError:
            print('\nPlease check the file path\n')
            
        if len(dfV['current'][~np.isnan(dfV['current'])]) < 30:
            print(f'\n### Not plotted {inFileV}\n')
            continue
           
        # plot the phase space probed by each exp
        plt.figure(3)
        plt.plot(dfV['bias'], 1e12*dfV['height'])
        
        # now get the info for the LDOS plot
        outerHeight.extend(dfV['height'].tolist())
        outerBias.extend(dfV['bias'].tolist())
        
        # find where current is below noise thresh 
        currentMask = abs(dfV['current']) < 5e-12
        # set current in this case to 0 as done in Feenstra, Phys. Rev. B, 1994 
        dfV['current'] = np.where(currentMask, 0, dfV['current'])
        
        # Still following Feenstra, convolve I/V with gaussian - width should be O(band gap)
        dfV['smoothCond'] = (ndimage.gaussian_filter1d(dfV['current']/dfV['bias'], 1.2))
        # print(dfV['smoothCond'].to_string())
        
        # Not following Feenstra, smooth dI/dV
        dfV['smoothLIX'] = (ndimage.gaussian_filter1d(dfV['lIX'], 1))
        
        lDOS = dfV['lIX']/dfV['smoothCond']
        # lDOS = dfV['smoothLIX']/dfV['smoothCond']
        diffCond = dfV['lIX']
        cond = dfV['smoothCond']
        
        outerLDOS.extend(lDOS.tolist())
        outerDiffCond.extend(diffCond.tolist())
        outerCond.extend(cond.tolist())
        
    # Save the phase space fig
    plt.figure(3)
    plt.xlabel('Bias (V)')
    plt.ylabel('Relative tip height (pm)')
    plt.savefig(f'figures/phaseSpace2dMap_{outstub}.png')
    
    # plot average LDOS
    # first make list of lists into a df
    lDOSdf = pd.DataFrame(np.column_stack([outerHeight, outerBias, outerLDOS]), columns=['height', 'bias', 'LDOS'])
    # find where LDOS is below 0
    lDOSMask = lDOSdf['LDOS'] < 0.0
    # set current in this case to 0 as -ve LDOS non-physical
    lDOSdf['LDOS'] = np.where(lDOSMask, 0, lDOSdf['LDOS'])
        
    # make two 2d hists. First sums the LDOS values from each exp in the bin
    # second counts the exps in the bin
    plt.figure(1)
    sumLDOS, biasEdges, heightEdges, im = plt.hist2d(lDOSdf['bias'], 1e12*lDOSdf['height'], weights=lDOSdf['LDOS'], bins=30)
    binCount, biasEdges, heightEdges, im = plt.hist2d(lDOSdf['bias'], 1e12*lDOSdf['height'], bins=30)
    plt.close()
    # can now normalise/find mean LDOS in each bin
    normLDOS = sumLDOS / binCount
    
    # finally display the average LDOS. 
    # Data needs to be transposed because of Python weirdness 
    plt.figure(2)
    plt.imshow(normLDOS.T, origin='lower',  cmap='viridis', aspect='auto', extent=(biasEdges[0], biasEdges[-1], heightEdges[0], heightEdges[-1]))
    plt.xlim(0, 3)
    plt.ylabel('Relative tip height (pm)', fontsize=20)
    plt.xlabel('Bias (V)', fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    cbar = plt.colorbar()
    cbar.set_label(r'(dI/dV)/$\overline{(I/V)}$')
    plt.savefig(f'figures/lDOS2dMap_{outstub}.png')
    plt.close()
    
    # plot average conductance
    # first make list of lists into a df
    conddf = pd.DataFrame(np.column_stack([outerHeight, outerBias, outerCond]), columns=['height', 'bias', 'cond'])
    # print(conddf.to_string())
    
    # make two 2d hists. First sums the cond values from each exp in the bin
    # Also convert from Siemens to G_0 -> 1 G_0 = 7.75x10^-5 S
    # second counts the exps in the bin
    plt.figure(1)
    sumcond, biasEdges, heightEdges, im = plt.hist2d(conddf['bias'], 1e12*conddf['height'], weights=7.75e5*conddf['cond'], bins=30)
    binCount, biasEdges, heightEdges, im = plt.hist2d(conddf['bias'], 1e12*conddf['height'], bins=30)
    plt.close()
    # can now normalise/find mean cond in each bin
    normcond = sumcond / binCount     
    
    # finally display the average cond. 
    # Data needs to be transposed because of Python weirdness 
    plt.figure(2)
    plt.imshow(normcond.T, origin='lower',  cmap='viridis', aspect='auto', extent=(biasEdges[0], biasEdges[-1], heightEdges[0], heightEdges[-1]))
    # plt.imshow(normcond.T, origin='lower',  cmap='viridis', aspect='auto', extent=(biasEdges[0], biasEdges[-1], heightEdges[0], heightEdges[-1]), norm=LogNorm())
    # plt.xlim(0, 3)
    plt.ylabel('Relative tip height (pm)')
    plt.xlabel('Bias (V)')
    cbar = plt.colorbar()
    cbar.set_label(r'Conductance ($G_0 = \frac{2e^2}{h}$)')
    plt.savefig(f'figures/cond2dMap_{outstub}.png')
    plt.close()

    # plot average diff conductance
    # first make list of lists into a df
    diffConddf = pd.DataFrame(np.column_stack([outerHeight, outerBias, outerDiffCond]), columns=['height', 'bias', 'diffCond'])
    # print(diffConddf.to_string())
    
    # make two 2d hists. First sums the diffCond values from each exp in the bin
    # second counts the exps in the bin
    plt.figure(1)
    sumdiffCond, biasEdges, heightEdges, im = plt.hist2d(diffConddf['bias'], 1e12*diffConddf['height'], weights=1e9*diffConddf['diffCond'], bins=30)
    binCount, biasEdges, heightEdges, im = plt.hist2d(diffConddf['bias'], 1e12*diffConddf['height'], bins=30)
    plt.close()
    # can now normalise/find mean diffCond in each bin
    normdiffCond = sumdiffCond / binCount     
    
    # finally display the average diffCond. 
    # Data needs to be transposed because of Python weirdness 
    plt.figure(2)
    plt.imshow(normdiffCond.T, origin='lower',  cmap='viridis', aspect='auto', extent=(biasEdges[0], biasEdges[-1], heightEdges[0], heightEdges[-1]))
    # plt.imshow(normdiffCond.T, origin='lower',  cmap='viridis', aspect='auto', extent=(biasEdges[0], biasEdges[-1], heightEdges[0], heightEdges[-1]), norm=LogNorm())
    # plt.xlim(0, 3)
    plt.ylabel('Relative tip height (pm)')
    plt.xlabel('Bias (V)')
    cbar = plt.colorbar()
    cbar.set_label(r'Differntial conductance (S)')
    plt.savefig(f'figures/diffCond2dMap_{outstub}.png')
    plt.close()
    return


@click.command()
@click.option('--infilez', '-iz', help="Input waveForm .txt file(s) with z lock-in modulated", type=str, prompt=False, multiple=True)
@click.option('--infilev', '-iv', help="Input waveForm .txt file(s) with V lock-in modulated", type=str, prompt=False, multiple=True)
@click.option('--outstub', '-o', help="Out file name stub", type=str, prompt=True, multiple=False)
# python3 constContour.py -i data/2023-06-20/IV_005.dat -i data/2023-06-20/IV_006.dat -i data/2023-06-20/IV_01\*.dat -o 200623
# cp -t target_directory foo_{0..54}.jpg
# python3 2dMap.py -iz data/2023-08-03/multi2_Inj_230803_\*_1_003.txt -iv data/2023-08-03/multi2_Inj_230803_\*_1_007.txt -o 230803_new

def main(infilez, infilev, outstub):
    kappaMap(infilez, outstub)
    lDOSMap(infilev, outstub)
    

if __name__=='__main__':
 	main()