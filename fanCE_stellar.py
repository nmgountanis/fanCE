import numpy as np


def weight(t,t0):
    
    """
    Compute weighting function
    
    Input arguments (illustrative values in [  ]):
      t = time array for desired outputs [0.0,0.02,0.04, ..., 12.0]
      t0 = age of universe at observed redshift [14]
           t0=0 returns weight=1, corresponding to mass-weighting

    Returns:
      weighting function, a 1-d array evaluated at the times in t. 

    Time unit is Gyr
    
    """
    
    if t0 == 0:
        return 1                               # weight=1 corresponds to mass-weighting
    
    
    else:
        return 1.5 * ((t0 - t + 0.1)**-0.78)   # light-weighting based off Gountanis+24 eq. 9
    
    
    
def stellar(t,t0,sfr,weight,OH,MgH,FeH,OFe,MgFe):
    
    """
    Compute light-weighted, log-averaged age and abundances.
    
    Input arguments (illustrative values in [  ]):
      t      = time array for desired outputs [0.0,0.02,0.04, ..., 12.0]
      t0     = age of universe at observed redshift [14]
      sfr    = star formation rate [4.03277751e-05,...7.50721612e-04]
      weight = light-weighting [0.19607826, 0.1963038, ..., 9.03839379]
      OH     = [O/H] evaluated at the times in t
      MgH    = [Mg/H] evaluated at the times in t
      FeH    = [Fe/H] evaluated at the times in t
      OFe    = [O/Fe] evaluated at the times in t
      MFe    = [Mg/Fe] evaluated at the times in t

    Returns:
      mean_age,OHStar,MgHStar,FeHStar,OFeStar, MgFeStar
      each a 1-d array evaluated at the times in t

    Time unit is Gyr
    
    """
                
    OHStar = np.cumsum(OH*sfr*weight)/np.cumsum(weight*sfr)   # light-weighted, log-averaged [O/H]
    
    MgHStar = np.cumsum(MgH*sfr*weight)/np.cumsum(weight*sfr) # light-weighted, log-averaged [Mg/H]
    
    FeHStar = np.cumsum(FeH*sfr*weight)/np.cumsum(weight*sfr) # light-weighted, log-averaged [Fe/H]
            
    OFeStar = OHStar - FeHStar                                # light-weighted, log-averaged [O/Fe]
    
    MgFeStar = MgHStar - FeHStar                              # light-weighted, log-averaged [Mg/Fe]
                
    
    # light-weighted, log-averaged age
    mean_age = 10**(np.cumsum(np.log10(t0-t+1e-5)*sfr*weight)/np.cumsum(sfr*weight)) 
    
    return mean_age,OHStar,MgHStar,FeHStar,OFeStar,MgFeStar
    
    
     