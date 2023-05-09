##############################################
# Get same phase in same colour?
##############################################

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
from scipy import ndimage

import commonFunctions as cf


def get_phase(filename):
    wSuff = filename.split('/')[-1]
    noSuff = wSuff.split('.d')[0]
    # stub = noSuff.split('phase_')[-1]
    stub = noSuff.split('V_')[-1]
    return stub

def integrate_LI(lI, z):
    runningIntegrand = []
    pointPairIntegral = []
    lIpairs = []
    zpairs = []
    # iterate through lists in pairs so have integral over the pairwise step region 
    for previous, current in zip(lI, lI[1:]):
        lIpairs.append([previous, current])
    for previous, current in zip(z, z[1:]):
        zpairs.append([previous, current])
    # print(f'lIp: {lIpairs}\nzp: {zpairs}')
    total = 0
    for lIp, zp in zip(lIpairs, zpairs):
        pointPInt = (2*np.sqrt(2))*scipy.integrate.simpson(lIp, zp)
        pointPairIntegral.append(abs(pointPInt))
        total += pointPInt
        runningIntegrand.append(abs(total))
    # print(runningIntegrand)
    # print(pointPairIntegral)
    return runningIntegrand, pointPairIntegral

def lockIn_investigate():
    # specObj  = nap.read.Spec('2023-02-16/Z-Spec_phase_-110.3_001.dat')

    # header = specObj.header

    # print(see_channels(specObj))
    # print(header)
    # sys.exit()

    # datFiles1 = glob('2023-02-24/Z-Spec_phase_1V_7*.dat')
    # datFiles2 = glob('2023-02-24/Z-Spec_phase_1V_8*.dat')
    # datFiles3 = glob('2023-02-24/Z-Spec_phase_1V_9*.dat')
    # datFiles = datFiles1 +  datFiles2 + datFiles3
    datFiles = glob('2023-02-24/Z-Spec_phase_1V_80_002.dat')
    # datFiles = glob('2023-02-24/Z-Spec_phase_1V_54.75_002.dat')
    datFiles.append('2023-02-24/Z-Spec_phase_1V_54.75_002.dat')

    # sort by phase value as a number not string
    datFiles.sort(key=lambda x: float(get_phase(x).split('_')[0]))
    print(datFiles)

    nFiles = len(datFiles)
    for datF in datFiles:
        print(datF)
        phase = get_phase(datF)
        print(phase)
        # sys.exit()

        specObj  = nap.read.Spec(datF)

        zData = cf.extract_channel(specObj, 'Z rel (m)')
        lIXdata = cf.extract_channel(specObj, 'LIX 1 omega (A)')
        lIXintegratedTotal, lIXintegratedPointPairs = integrate_LI(lIXdata*(10**12), zData*(10**12))
        lIYdata = cf.extract_channel(specObj, 'LIY 1 omega (A)')
        Idata = cf.extract_channel(specObj, 'Current (A)')

        color = plt.cm.plasma(np.linspace(0.1, 0.9, nFiles))
        mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

        plt.figure(1)
        plt.plot(zData*(10**12), lIXdata*(10**12), label=phase)

        plt.figure(2)
        plt.plot(zData*(10**12), lIYdata*(10**12), label=phase)

        plt.figure(3)
        plt.plot(zData*(10**12), Idata*(10**12), label=phase)

        plt.figure(4)
        # plt.plot(zData[:-1]*(10**12), lIXintegratedPointPairs, label=phase)
        plt.plot(zData[:-1]*(10**12), lIXintegratedPointPairs/(Idata[:-1]*(10**12)), label=phase)
        
    plt.figure(1)
    plt.xlabel('z (pm)')
    plt.ylabel('LIX (pA)')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    plt.legend(loc='upper right', prop=fontP)
    plt.savefig('../figures/lIX_1V_cap.png')
    plt.close()

    plt.figure(2)
    plt.xlabel('z (pm)')
    plt.ylabel('LIY (pA)')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    plt.legend(loc='upper right', prop=fontP)
    plt.savefig('../figures/lIY_1V_cap.png')
    plt.close()

    plt.figure(3)
    plt.xlabel('z (pm)')
    plt.ylabel('I (pA)')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    plt.legend(loc='upper right', prop=fontP)
    plt.savefig('../figures/lI_I_1V_cap.png')
    plt.close()

    plt.figure(4)
    plt.xlabel('z (pm)')
    plt.ylabel('I_x (pA) ')
    # plt.ylabel('I_x / I')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    plt.legend(loc='upper right', prop=fontP)
    # plt.savefig('../figures/lI_Ix_1V_cap.png')
    plt.savefig('../figures/lI_Ix_I_1V_cap.png')
    plt.close()
    
    return
#################################### end lockIn_investigate

def get_osc_height(filename):
    wSuff = filename.split('/')[-1]
    noSuff = wSuff.split('.d')[0]
    # stub = noSuff.split('phase_')[-1]
    # stub = noSuff.split('I_')[-1]
    stub = noSuff.split('80.3_')[-1]
    print(stub)
    return stub

def capCurrent_investigate():

    # datFiles = glob('2023-03-13/Oscilloscope_I_10mV_001.dat')
    datFiles = glob('2023-03-13/Oscilloscope_I_*mV_*.dat')

    # sort by amplitude as a number not string
    datFiles.sort(key=lambda x: float(get_osc_height(x).split('m')[0]))
    print(datFiles)

    rmsVals = []
    heights = []

    nFiles = len(datFiles)
    for datF in datFiles:
        print(datF)
        height = get_osc_height(datF)

        specObj  = nap.read.Spec(datF)

        timeData = cf.extract_channel(specObj, 'Time (s)')
        currentData = cf.extract_channel(specObj, 'Current (A)')

        color = plt.cm.plasma(np.linspace(0.1, 0.9, nFiles))
        mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

        # plot raw currnt
        plt.figure(1)
        plt.plot(timeData*(10**3), currentData*(10**12), '-', label=height, lw=0.5)
        plt.xlabel('Time (ms)')
        plt.ylabel('Current (pA)')
        plt.ylim(-200, 200)
        fontP = FontProperties() # Making legend smaller
        fontP.set_size('x-small')
        plt.legend(loc='upper right', prop=fontP)
        plt.savefig(f'../figures/cap_I_{height}.png')
        plt.close()

        # find rms current as function height of oscillation
        heights.append(float(get_osc_height(datF).split('m')[0]))
        rmsVals.append((10**12)*np.sqrt(np.sum(currentData**2)/len(currentData)))

    plt.figure(2)
    plt.plot(heights, rmsVals, 'x')
    plt.xlabel('Lock-in oscillation height (mV)')
    plt.ylabel('rms(I) (pA)')
    plt.xlim(3, 22)
    plt.savefig('../figures/cap_I_rms_heights.png')
    plt.close()


    return

def kappa_plot():
    datFiles = glob('2023-03-22/Z-Spec_phase_80.3_-205*.dat')
    
    datFiles.sort(key=lambda x: float(get_osc_height(x).split('p')[0]))
    print(datFiles)
    
    kappaList = []
    
    nFiles = len(datFiles)
    for datF in datFiles:
        
        height = get_osc_height(datF)

        specObj  = nap.read.Spec(datF)

        zData = cf.extract_channel(specObj, 'Z rel (m)')
        lIXdata = cf.extract_channel(specObj, 'LIX 1 omega (A)')
        lIYdata = cf.extract_channel(specObj, 'LIY 1 omega (A)')
        Idata = cf.extract_channel(specObj, 'Current (A)')
        
        # Check got a spectrum, not tip withdrawn
        # remove NaNs from the df and check have sufficient data points left
        # Idata = Idata[~np.isnan(Idata)]
        # print(Idata)
        if len(Idata[~np.isnan(Idata)]) < 10:
            print(f'\n### Not plotted {datF}\n')
            continue

        color = plt.cm.plasma(np.linspace(0.1, 0.9, nFiles))
        mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)
        
        # kappa = -lIXdata/(2*Idata)
        kappa = lIXdata/(2*Idata)
        kappaList.append(kappa)
        
        plt.figure(1)
        plt.plot(zData*(10**12), kappa, label=height)
        
    plt.figure(1)
    plt.xlabel('z (pm)')
    plt.ylabel('kappa')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    plt.legend(loc='upper right', prop=fontP)
    plt.savefig('../figures/kappa.png')
    plt.close()

    print(kappaList)
    
    return


def lDOS_plot():
    
    waveTxtFiles = glob('data/2023-05-03/multi2_Inj_230503_4_*.txt')
    
    outerLDOS = [] # hold the LDOS for each point in each experiment
    
    for wvTxtF in waveTxtFiles:
        try:
            df = cf.extract_waveTxt(wvTxtF)
        except ValueError:
            print('Please check file path')
            # print(e)
            
        fPath = wvTxtF.split('.')[0]
        fName = fPath.split('/')[-1]
        
        # find where current is below noise thresh 
        currentMask = abs(df['current']) < 5e-12
        # set current in this case to 0 as done in Feenstra, Phys. Rev. B, 1994 
        df['current'] = np.where(currentMask, 0, df['current'])
        
        # setup figure to hold the IV, dI/dV, LDOS plots
        fig1 = plt.figure(1, figsize=(5.75, 8))
        ax_dIdV, ax_IV, ax_LDOS = fig1.subplots(3, 1, sharex=True)
        
        ax_IV.plot(df['bias'], 10e12*df['current']/df['bias'], label='Raw')
        
        # Still following Feenstra, convolve I/V with gaussian - width should be O(band gap)
        df['smoothCond'] = (ndimage.gaussian_filter1d(df['current']/df['bias'], 1.2))
        # print(df['smoothCond'].to_string())
        
        ax_IV.plot(df['bias'], 10e12*df['smoothCond'], label='Smoothed')
        ax_IV.set_ylabel('Conductance (pA/V)')
        ax_IV.set_ylim(bottom=0)
        fontP = FontProperties() # Making legend smaller
        fontP.set_size('x-small')
        ax_IV.legend(loc='lower right', prop=fontP)
        
        # Not following Feenstra, smooth dI/dV
        df['smoothLIX'] = (ndimage.gaussian_filter1d(df['lIX'], 1))
        
        ax_dIdV.plot(df['bias'], df['lIX'])
        ax_dIdV.set_ylabel('dI/dV (pA/V)')
        ax_dIdV.set_yscale('log')
        
        lDOS = df['lIX']/df['smoothCond']
        ax_LDOS.plot(df['bias'], lDOS)
        ax_LDOS.set_ylabel(r'(dI/dV)/$\overline{(I/V)}$')
        ax_LDOS.set_ylim(0, 10)
        
        outerLDOS.append(lDOS)
        
        
        plt.figure(1)
        plt.xlabel('Bias (V)')
        plt.savefig(f'figures/combined_lDOS_{fName}.png')
        plt.close()
        
        # print(df.to_string())
    
    # plot average LDOS
    # first make list of lists into a df
    # transpose df so each column is a single experiment
    # doesn't strictly have to be done (would average over the other axis) but made sense to me 
    lDOSdf = pd.DataFrame(outerLDOS).T
    lDOSdf['mean'] = lDOSdf.mean(axis=1)
    
    
    plt.figure(2)
    plt.plot(df['bias'], lDOSdf['mean'])
    plt.ylabel(r'(dI/dV)/$\overline{(I/V)}$')
    plt.xlabel('Bias (V)')
    # plt.ylim(0, 15)
    plt.savefig(f'figures/ave_lDOS.png')
    plt.close()
    
    return


def main():
    
    # lockIn_investigate()
    # capCurrent_investigate()
    # kappa_plot()
    lDOS_plot()
    

if __name__=='__main__':
 	main()