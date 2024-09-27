#########################################################
# Holds functions used across many analysis scripts.
# These will be mostly related to extracting data from files into DataFrames.
# The DataFrames are then returned to the main analysis scripts.
# Should gracefully deal with missing/incorrect files etc.
# ------------------------------------------------------- 
# Usually used as `import commonFunctions as cf`
# -------------------------------------------------------
# TODO - 
#########################################################

import sys
import numpy as np
import pandas as pd
from glob import glob
import nanonispy as nap

def see_channels(napObj):
    """List the available data channels in a nanonis (.dat) file.   
	
	Arguments
	---------
	napObj : nanonispy.read.Spec
		nanonispy spectrum object from the dat file to be read
	
	Returns
	------
	napObj.signals.keys() : dict_keys([str])
        The available data channels in the dat files
	"""
 
    return napObj.signals.keys()

def extract_channel(napObj, chan):
    """Extract the data from a given channel in a nanonis (.dat) file.   
	
	Arguments
	---------
	napObj : nanonispy.read.Spec
		nanonispy spectrum object from the dat file to be read
    chan : str
        The data channel to be read
	
    Returns
	------
	napObj.signals.get(chan) : numpy.ndarray[float]
        The extracted data from the channel
	"""
 
    return napObj.signals.get(chan)

class UnhandledFileError(Exception):
    """
    To be raised when unknown file extension or wrong file format is passed.
    """
    pass

def extract_file(file):
    """Extract data from LabView (waveTxt) and nanonis (.dat) files.    
	
	Arguments
	---------
	file : str
		Filename string to be read
	
	Returns
	------
	df : pd.DataFrame
        Data frame with extracted data
	"""
 
    if file.split('.')[-1] == 'dat':
        try:
            df = extract_datFile(file)
        except ValueError as err:
            raise ValueError(f'\n#####\nFile `{file}` is potentially damaged\nValueError: {err}\n#####\n') # e.g. missing columns or values
    elif file.split('.')[-1] == 'txt':
        print('its wvTxt')
        try:
            df = extract_waveTxt(file)
        except ValueError as err:
            raise ValueError(f'\n#####\nFile `{file}` is potentially damaged\nValueError: {err}\n#####\n') # e.g. missing columns or values
    else: # handle non standard file type 
        raise UnhandledFileError(f"\n#####\n`{file}` is not of a supported file type. Please check your input.\n#####\n")
     
    return df

def extract_datFile(datFile):
    """Extract data from a nanonis (.dat) file. 
	
	Arguments
	---------
	file : str
		Filename string to be read
	
	Returns
	------
	df : pd.DataFrame
        Data frame with extracted data
	"""
 
    print(f'Extracting {datFile}')
    
    # Get the specObject
    specObj  = nap.read.Spec(datFile)
    chans = see_channels(specObj)
    df = pd.DataFrame()
    for chan in chans:
        df[f'{chan}'] = extract_channel(specObj, chan)
    
    # Rename columns to be consistent with wvTxt output
    df.columns = ['height', 'current', 'bias', 'lIX', 'lIY']
    # switch sign of lI data to be consistent with wvTxt
    df['lIX'] = -df['lIX']
    df['lIY'] = -df['lIY']
    
    return df

def extract_waveTxt(wvFile):
    """Extract data from a LabView (waveTxt) file.    
	
	Arguments
	---------
	file : str
		Filename string to be read
	
	Returns
	------
	df : pd.DataFrame
        Data frame with extracted data
	"""
    
    print(f'Extracting {wvFile}')
    datetime, height, bias, current, lIX = np.loadtxt(wvFile, usecols=(0,1,2,3,4), skiprows=5, unpack=True, delimiter='\t', dtype=str)
    
    #  Time info comes in a weird date and clock time format...
    # ...have to read in all data as a string to use string format methods. 
    # Get time in seconds from start of exp and convert to float
    date_times = np.char.split(datetime)
    ftr = [3600,60,1] # needed to convert hms to s
    timeList = []
    for dt in date_times:
        time_hms = dt[-1]
        time_s = sum([a*b for a,b in zip(ftr, map(float, time_hms.split(':')))])
        timeList.append(float(time_s))
    timeArr = np.array(timeList)
    # find 'zero' time and express others relative to that
    time = timeArr - min(timeArr)
    
    #  Convert rest of data to floats 
    height = np.asarray(height, dtype=float)
    bias = np.asanyarray(bias, dtype=float)
    current = np.asanyarray(current, dtype=float)
    lIX = np.asanyarray(lIX, dtype=float)
    
    # and then put together into a df
    df = pd.DataFrame(data={'time': time, 'height': height, 'bias': bias, 'current': current, 'lIX': lIX})
    
    # print(time)
    # print(height)
    # print(bias)
    # print(current)
    # print(lIX)
    
    return(df)

def globFiles(infiles):
    """Extract the file names to be read from cmd line input.
        Designed to handle combinations of multiple fully named files and inputs including wildcards et al. 
        Exit if input files do not exist.   
        Globbing done in helper function (glob_fileStr) to allow error handling 
	
	Arguments
	---------
	infiles : (str)
		Tuple of cmd line file string(s) to be extracted
	
	Returns
	------
	inFiles : [str]
        List with each element a file to be read from 
	"""

    inFiles = []
    print(f'infile(s): {infiles}')
    
    for inFileStr in infiles:
        try: # attempt extract for each of the file strings and add to output list if found
            inFiles.extend(glob_fileStr(inFileStr))
        except FileNotFoundError as err:
            print(err)
        
    return inFiles

def glob_fileStr(fileStr):
    """Helper function to globFiles
        Extract the file names to be read from a given file string
        Designed to handle wildcards et al. 
        Raise error if input files do not exist.     
	
	Arguments
	---------
	fileStr : str
		A single cmd line file string(s) to be extracted
	
	Returns
	------
	inFs : [str]
        List with each element a file to be read from 
	"""
 
    inFs = glob(fileStr)
    if not inFs:
        raise FileNotFoundError(f"\n#####\nNo file matching the path `{fileStr}` was found. Please check your input.\n#####\n")
    
    return inFs
