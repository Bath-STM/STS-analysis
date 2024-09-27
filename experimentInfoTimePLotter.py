#########################################################
# Explainations
# ------------------------------------------------------- 
# Example command ``
# -------------------------------------------------------
# TODO - 
#########################################################
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import sys
import numpy as np
import math
from matplotlib.font_manager import FontProperties
import click
import matplotlib.ticker as ticker

import commonFunctions as cf

def infoPlot(infiles, outstub):
    
    if not infiles:
        return
    
    inFiles = cf.globFiles(infiles)
    
    # setup figure to hold the current, bias, lockIn & height plots
    fig1 = plt.figure(1, figsize=(7, 15))
    # fig1 = plt.figure(1, figsize=(7, 10))
    ax_I, ax_V, ax_lI, ax_z = fig1.subplots(4, 1, sharex=True)
        
    for inFile in inFiles:
        try:
            df = cf.extract_file(inFile)
        except ValueError:
            print('\nPlease check the file path\n')
            
        # Check got a spectrum, not tip withdrawn
        # remove NaNs from the df and check have sufficient data points left
        if len(df['current'][~np.isnan(df['current'])]) < 30:
            print(f'\n### Not plotted {inFile}\n')
            continue
        
        print(df)
        
        # Plot current
        ax_I.plot(df['time'], 1e12*df['current'], alpha=0.7)
        ax_I.set_ylabel('Current (pA)')
        
        # Plot bias
        ax_V.plot(df['time'], df['bias'], alpha=0.7)
        ax_V.set_ylabel('Bias (V)')
        
        # Plot lockin
        ax_lI.plot(df['time'], df['lIX'], alpha=0.7)
        ax_lI.set_ylabel(r'$\partial I / \partial z$')
        
        # Plot z height
        ax_z.plot(df['time'], 1e12*df['height'], alpha=0.7)
        ax_z.set_ylabel('z height (pm)')
        
    plt.xlabel('Time (s)')
    ax = plt.gca()
    yloc = ticker.MultipleLocator(base=50.0) # this locator puts ticks at regular intervals
    ax.yaxis.set_major_locator(yloc)
    n = 2  # Keeps every 2nd tick label - use 0 or 1 at end next line for where want start
    [l.set_visible(False) for (i,l) in enumerate(ax.yaxis.get_ticklabels()) if i % n != 0]
    plt.savefig(f'figures/expInfo_{outstub}.png')
    plt.close()
    return

@click.command()
@click.option('--infiles', '-i', help="Input file(s)", type=str, prompt=True, multiple=True)
@click.option('--outstub', '-o', help="Out file name stub", type=str, prompt=True, multiple=False)
# python3 pedestalPlotter.py -iz data/2024-03-27/multi2_Inj_240327_\[2-9]_\*.txt -iz data/2024-03-27/multi2_Inj_240327_1[0-9]_\*.txt -o 240327_2-19

def main(infiles, outstub):
    
    infoPlot(infiles, outstub)
    
if __name__=='__main__':
 	main()