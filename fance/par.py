import numpy as np

def DefaultfanCEParSet():
    return {
        't_start': 0.5,
        't0': 14.0,
        'dt': 0.02,
        'tauSFE': 1.0,
        'tauSFH1': 2.0,
        'tauSFH2': 8.0,
        'aFeCC': 0.45, # plateau
        'aFeEq': 0.0, # late-time equilibrium
        'mu': 1.1, 
        'yOCC': None, # computed in fance using aFeEq and mu if set to None
        'yMgCC': None, # computed in fance using aFeEq and mu if set to None
        'yFeCC': None, # computed in fance using aFeEq and mu if set to None
        'yFeIa': None, # computed in fance using aFeEq and mu if set to None
        'fRetCC': 1.0,
        'r': 0.4,
        'eta': 0.3,
        'tauIa': 0,
        'tDminIa': 0.15,
        'SolarO': 7.33e-3,
        'SolarMg': 6.71e-4,
        'SolarFe': 1.37e-3,
        'IaDTD_fn': "powerlaw",
    }
