#########################################################
# Explainations
# ------------------------------------------------------- 
# Example command ``
# -------------------------------------------------------
# TODO - 
#########################################################

##############################################
# Get same phase in same colour?
# Use new file reading in commonFunctions
##############################################

from glob import glob
from turtle import right
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
    return np.array(runningIntegrand), np.array(pointPairIntegral)

def lockIn_investigate():
    # specObj  = nap.read.Spec('2023-02-16/Z-Spec_phase_-110.3_001.dat')

    # header = specObj.header

    # print(see_channels(specObj))
    # print(header)
    # sys.exit()

    # inFiles1 = glob('2023-02-24/Z-Spec_phase_1V_7*.dat')
    # inFiles2 = glob('2023-02-24/Z-Spec_phase_1V_8*.dat')
    # inFiles3 = glob('2023-02-24/Z-Spec_phase_1V_9*.dat')
    # inFiles = inFiles1 +  inFiles2 + inFiles3
    inFiles = ['../firstLooks/data/2023-02-24/Z-Spec_phase_1V_80_002.dat']
    # inFiles = glob('2023-02-24/Z-Spec_phase_1V_54.75_002.dat')
    # inFiles.append('2023-02-24/Z-Spec_phase_1V_54.75_002.dat')

    # sort by phase value as a number not string
    # inFiles.sort(key=lambda x: float(get_phase(x).split('_')[0]))
    
    # inFiles = ['data/2024-03-27/multi2_Inj_240327_7_1_004.txt']
    # print(inFiles)

    nFiles = len(inFiles)
    for inF in inFiles:
        print(inF)
        # phase = get_phase(inF)
        phase = 'data'
        # print(phase)
        # # sys.exit()

        # specObj  = nap.read.Spec(inF)

        # zData = cf.extract_channel(specObj, 'Z rel (m)')
        # lIXdata = cf.extract_channel(specObj, 'LIX 1 omega (A)')
        # lIXintegratedTotal, lIXintegratedPointPairs = integrate_LI(lIXdata*(10**12), zData*(10**12))
        # lIYdata = cf.extract_channel(specObj, 'LIY 1 omega (A)')
        # Idata = cf.extract_channel(specObj, 'Current (A)')
        
        try:
            df = cf.extract_file(inF)
        except ValueError:
            print('\n############\nPlease check the file path\n############\n')

        color = plt.cm.plasma(np.linspace(0.1, 0.9, nFiles))
        mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

        plt.figure(1)
        plt.plot(df['height']*(10**12), df['lIX']*(10**12), label=phase)

        # plt.figure(2)
        # plt.plot(df['height']*(10**12), df['lIY']*(10**12), label=phase)

        plt.figure(3)
        print(df['current']*(10**12))
        plt.plot(df['height']*(10**12), df['current']*(10**12), 'x', label=phase)

        lIXintegratedTotal, lIXintegratedPointPairs = integrate_LI(df['lIX']*(10**12), df['height']*(10**12))
        # print(1e-10*np.array(lIXintegratedPointPairs))
        plt.figure(4)
        # plt.plot(df['height'][:-1]*(10**12), lIXintegratedPointPairs, label=phase)
        plt.plot(df['height'][:-1]*(10**12), lIXintegratedPointPairs/(df['current'][:-1]*(10**12)), label=phase)
        
    plt.figure(1)
    plt.xlabel('z (pm)')
    plt.ylabel('LIX (pA)')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    plt.legend(loc='upper right', prop=fontP)
    plt.savefig('figures/lIX_1V_cap.png')
    plt.close()

    plt.figure(2)
    plt.xlabel('z (pm)')
    plt.ylabel('LIY (pA)')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    plt.legend(loc='upper right', prop=fontP)
    plt.savefig('figures/lIY_1V_cap.png')
    plt.close()

    plt.figure(3)
    plt.xlabel('z (pm)')
    plt.ylabel('I (pA)')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    # plt.legend(loc='upper right', prop=fontP)
    plt.savefig('figures/lI_I_1V_cap.png')
    plt.close()

    plt.figure(4)
    plt.xlabel('z (pm)')
    plt.ylabel('I_x (pA) ')
    # plt.ylabel('I_x / I')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    plt.legend(loc='upper right', prop=fontP)
    # plt.savefig('../figures/lI_Ix_1V_cap.png')
    plt.savefig('figures/lI_Ix_I_1V_cap.png')
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
    
    # inFiles = glob('data/2023-05-10/IV_modA_10mV_*.dat')
    # Adatom push
    inFiles = glob('data/2023-03-22/*.dat')
    # Transfer bromo
    # inFiles = glob('data/2022-12-14/multi2_Inj_221214_2_[2-4]_*.txt')
    # Transfer dark noBromo
    # inFiles = glob('data/2022-12-15/multi2_Inj_221215_1_[3-5]_*.txt')
    # TRansfer bright noBromo
    # inFiles = glob('data/2022-12-14/multi2_Inj_221214_1_[1-3]_*.txt')
    print(inFiles)
    nFiles = len(inFiles)
    
    kappaList = []
    
    if inFiles[0].split('.')[-1] == 'dat':
        print('dat')
        
        inFiles.sort(key=lambda x: float(get_osc_height(x).split('p')[0]))
        
        for datF in inFiles:
        
            height = get_osc_height(datF)

            specObj  = nap.read.Spec(datF)

            zData = cf.extract_channel(specObj, 'Z rel (m)')
            # Adatom calibration
            zData = zData + 750e-12
            lIXdata = cf.extract_channel(specObj, 'LIX 1 omega (A)')
            Idata = cf.extract_channel(specObj, 'Current (A)')
            
            # Check got a spectrum, not tip withdrawn
            # remove NaNs from the df and check have sufficient data points left
            if len(Idata[~np.isnan(Idata)]) < 30:
                print(f'\n### Not plotted {datF}\n')
                continue

            color = plt.cm.plasma(np.linspace(0.1, 0.9, nFiles))
            mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)
            
            print(f'\nraw lIX {lIXdata}')
            print(f'\nDivided lIX {lIXdata/10e-12}')
            lIXdata = lIXdata/10e-12
            # sys.exit()
            
            # These out by a factor? lIX -> delta I so need to / by delta z (x pm)....
            # ... this depends what is output in the multi_Inj_xxx.txt files and if 
            # this scaling by the amplitude is already done...?
            # delta z = 10pm for most data
            # kappa = -lIXdata/(2*Idata)
            # kappa = lIXdata/(2*Idata)
            kappa = (lIXdata)/(2*Idata)
            # kappa = -lIXdata/(2*Idata)
            # kappa = lIXdata/(2*Idata)
            kappaList.append(kappa)
            
            plt.figure(1)
            plt.plot(1e12*zData, 1e-9*kappa, label=height)
            
    elif inFiles[0].split('.')[-1] == 'txt':
        print('txt')
        
        for wvTxtF in inFiles:
            try:
                df = cf.extract_waveTxt(wvTxtF)
            except ValueError:
                print('Please check file path')
                # print(e)
                
            zData = df['height']
            # Optional offset to calibrated height when above dark mol (Bromo).
            zData = zData + 750e-12
            # Optional offset to calibrated height when above dark contaminant (noBromo).
            # zData = zData - 70e-12
            lIXdata = df['lIX']
            Idata = df['current']
            
            # Check got a spectrum, not tip withdrawn
            # remove NaNs from the df and check have sufficient data points left
            if len(Idata[~np.isnan(Idata)]) < 10:
                print(f'\n### Not plotted {datF}\n')
                continue

            # color = plt.cm.plasma(np.linspace(0.1, 0.9, nFiles))
            # mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)
            
            # These out by a factor? lIX -> delta I so need to / by delta z (x pm)....
            # ... this depends what is output in the multi_Inj_xxx.txt files and if 
            # this scaling by the amplitude is already done...?
            # delta z = 10pm for most data
            # kappa = -lIXdata/(2*Idata)
            kappa = lIXdata/(2*Idata)
            kappa = (lIXdata/10e-12)/(2*Idata)
            kappaList.append(kappa)
            
            plt.figure(1) # plot kappa
            plt.plot(1e12*zData, 1e-9*kappa, 'k')
            plt.figure(2) # plot raw current
            plt.plot(1e12*zData, 1e12*Idata)
            plt.figure(3) # plot effective barrier/ work function
            # kappa = sqrt(2 m_e phi)/ hbar
            # (kappa * hbar)^2 / 2 m_e
            hbar = 6.626e-34 / (2*np.pi)
            m_e = 9.109e-31
            e = 1.602e-19
            phi_eff = (kappa * hbar)**2 / (2*m_e)
            plt.plot(1e12*zData, phi_eff/e)
        
    else:
        print('\nUnknown file type - exiting\n')
        sys.exit(1)
    
    plt.figure(1)
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel(r'$\kappa$ (nm$^{-1}$)', fontsize=20)
    # plt.xlim(-180, 60) # scaled for darkBromo
    plt.xlim(520, 800)
    plt.ylim(bottom=0)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    # plt.yscale('log')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    # plt.legend(loc='upper right', prop=fontP)
    plt.savefig('figures/kappa_230322_push.png')
    # plt.savefig('figures/kappa_cleanFC_Bromo_221214.png')
    # plt.savefig('figures/kappa_darkFC_noBromo_221215.png')
    plt.close()
    
    plt.figure(2)
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel('I (pA)', fontsize=20)
    # plt.xlim(-180, 60) # scaled for darkBromo
    plt.xlim(-190, 110)
    plt.ylim(bottom=3)
    plt.tight_layout()
    plt.yscale('log')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    # plt.legend(loc='upper right', prop=fontP)
    # plt.savefig('figures/kappa_230322_push.png')
    # plt.savefig('figures/kappa_cleanFC_Bromo_221214.png')
    # plt.savefig('figures/current_darkFC_Bromo_221214.png')
    plt.close()
    
    plt.figure(3)
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel(r'$\phi_{eff}$ (eV)', fontsize=20)
    # plt.xlim(-180, 60) # scaled for darkBromo
    plt.xlim(-190, 110)
    # plt.ylim(bottom=3)
    plt.tight_layout()
    # plt.yscale('log')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('x-small')
    # plt.legend(loc='upper right', prop=fontP)
    # plt.savefig('figures/kappa_230322_push.png')
    # plt.savefig('figures/kappa_cleanFC_Bromo_221214.png')
    # plt.savefig('figures/phiEff_darkFC_Bromo_221214.png')
    plt.close()

    # print(kappaList)
    
    return


def lDOS_plot():
    
    waveTxtFiles = glob('data/2023-07-11/multi2_Inj_230711_1_*.txt')
    waveTxtFiles += glob('data/2023-07-11/multi2_Inj_230711_5_*.txt')
    print(waveTxtFiles)
    outerLDOS = [] # hold the LDOS for each point in each experiment
    
    for wvTxtF in waveTxtFiles:
        try:
            df = cf.extract_waveTxt(wvTxtF)
        except ValueError:
            print('Please check file path')
            # print(e)
            
        fPath = wvTxtF.split('.')[0]
        fName = fPath.split('/')[-1]
        aveFileName = fName.split("_")[2] + '_' + fName.split("_")[3]
        
        print(aveFileName)
        
        # find where current is below noise thresh 
        currentMask = abs(df['current']) < 5e-12
        # set current in this case to 0 as done in Feenstra, Phys. Rev. B, 1994 
        df['current'] = np.where(currentMask, 0, df['current'])
        
        # setup figure to hold the IV, dI/dV, LDOS plots
        fig1 = plt.figure(1, figsize=(5.75, 8))
        ax_dIdV, ax_IV, ax_LDOS = fig1.subplots(3, 1, sharex=True)
        
        ax_IV.plot(df['bias'], 1e12*df['current']/df['bias'], label='Raw')
        
        # Still following Feenstra, convolve I/V with gaussian - width should be O(band gap)
        df['smoothCond'] = (ndimage.gaussian_filter1d(df['current']/df['bias'], 1.2))
        # print(df['smoothCond'].to_string())
        
        ax_IV.plot(df['bias'], 1e12*df['smoothCond'], label='Smoothed')
        ax_IV.set_ylabel('Conductance (pA/V)')
        ax_IV.set_ylim(bottom=0)
        fontP = FontProperties() # Making legend smaller
        fontP.set_size('x-small')
        ax_IV.legend(loc='lower right', prop=fontP)
        
        # Not following Feenstra, optionally smooth dI/dV
        df['smoothLIX'] = (ndimage.gaussian_filter1d(df['lIX'], 1))
        
        ax_dIdV.plot(df['bias'], df['lIX'])
        # ax_dIdV.plot(df['bias'], df['smoothLIX'], label='Smoothed')
        ax_dIdV.set_ylabel('dI/dV (pA/V)')
        # ax_dIdV.set_yscale('log')
        # ax_dIdV.legend(loc='lower right', prop=fontP)
        
        lDOS = df['lIX']/df['smoothCond']
        # lDOS = df['smoothLIX']/df['smoothCond']
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
    
    print(lDOSdf)
    
    # plt.figure(2)
    plt.figure(2, figsize=(5,3), dpi=600)
    plt.plot(df['bias'], lDOSdf['mean'])
    plt.ylabel(r'(dI/dV)/$\overline{(I/V)}$')
    plt.xlabel('Bias (V)')
    # plt.ylim(0, 15)
    plt.ylim(top=8)
    plt.tight_layout()
    plt.savefig(f'figures/ave_lDOS_{aveFileName}.png')
    plt.savefig(f'figures/ave_lDOS_{aveFileName}.pdf')
    plt.close()
    
    return


def main():
    
    lockIn_investigate()
    # capCurrent_investigate()
    # kappa_plot()
    # lDOS_plot()
    

if __name__=='__main__':
 	main()