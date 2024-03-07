import numpy as np

'''
All timescales are given in Gyr
'''

def time(tstart, t0, dtout):
     
    """

    Input arguments (illustrative values in [  ]):
      tstart = star formation start time [0.5]
      t0 = maximum time [14]
      dtout = output timestep [0.02]

    Returns:
      tmod + tstart, a 1-d time array that goes from t_start + dtout to t0 with a
                     timestep of dtout
      tmod, a 1-d time array 
      
      

    Time unit is Gyr
    
    """
    # Define tmax 
    tmax = t0 - tstart
    
    # Define tmod
    tmod = np.arange(dtout, tmax+1.e-5,dtout)

    # return time array
    return tmod + tstart, tmod
    