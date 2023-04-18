##################################################
# TODO sort error handling
# TODO Add function defintion documentaion blocks
##################################################

import sys
import numpy as np
import pandas as pd

def see_channels(napObj):
    return napObj.signals.keys()

def extract_channel(napObj, chan):
   return napObj.signals.get(chan)

def apply_mask_to_all(mask, *arrays):
    
    assert all([arr.shape == mask.shape for arr in arrays]), "All Arrays need to have the same shape as the mask"
    
    return tuple([arr[mask] for arr in arrays])

def check_waveTxt(wvFile):
    if wvFile.split('.')[-1] == 'txt':
        # sys.exit(1) # or better raise an exception...
        # ...file not found/wrong file format class with file name printed to screen
        # raise UnhandledFileError(f"{wvFile} is not a .txt file as required. Convert to txt with labview")
        pass
    else:
        raise ValueError(f"{wvFile} is not a .txt file as required. Convert to txt with labview")
    
    # return

def extract_waveTxt(wvFile):
    
    # try:
    #     datetime, height, bias, current, lIX = np.loadtxt(wvFile, usecols=(0,1,2,3,4), skiprows=5, unpack=True, delimiter='\t', dtype=str)
    # except Exception as e:
    #     print(f"{wvFile} is not a .txt file as required. Convert to txt with labview")
    #     raise e
         
    # if wvFile.split('.')[-1] == 'txt':
    #     # sys.exit(1) # or better raise an exception...
    #     # ...file not found/wrong file format class with file name printed to screen
    #     # raise UnhandledFileError(f"{wvFile} is not a .txt file as required. Convert to txt with labview")
    #     pass
    # else:
    #     raise ValueError(f"{wvFile} is not a .txt file as required. Convert to txt with labview")
    
    check_waveTxt(wvFile)
    
    print(f'Extracting {wvFile}')
    datetime, height, bias, current, lIX = np.loadtxt(wvFile, usecols=(0,1,2,3,4), skiprows=5, unpack=True, delimiter='\t', dtype=str)
    
    #  Time info comes in a weird date and clock time format...
    # ..have to read in all data as a string to use string format methods. 
    # Get time in seconds and convert to float
    times = np.char.split(datetime)
    timeList = []
    for t in times:
        timeList.append(float(t[-1].split(':')[-1]))
    time = np.array(timeList)
    
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
    
    # print(df)
        
    return(df)


class UnhandledFileError(Exception):

    """
    To be raised when unknown file extension or wrong file format is passed.
    """
    pass