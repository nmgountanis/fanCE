import numpy as np

def DefaultWAFParSet(dt=0.02, t0=14, t_start=0.5):
    tmod = np.arange(dt, t0 - t_start + dt, dt)

    return {
        'tmod': tmod, 
        'tauSFE': 1.0,
        'tauSFH': 8.0,
        'yOCC': 0.0071,
        'yMgCC': 0.0007,
        'yFeCC': 0.0005,
        'yFeIa': 0.0007,
        'r': 0.4,
        'eta': 0.3,
        'tauIa': 1.5,
        'tDminIa': 0.15,
        'SolarO': 0.0073,
        'SolarMg': 0.0007,
        'SolarFe': 0.0014,
        'SFH_fn': "exponential",
        'IaDTD_fn': "powerlaw",
    }

    
def DefaultfanCEParSet():
    return {
        "t_start": 0.5,
        "t0": 14.0,
        "dt": 0.02,
        "tauSFE": 1.0,
        "tauSFH1": 2.0,
        "tauSFH2": 8.0,
        "yOCC": 0.0071,
        "yMgCC": 0.0007,
        "yFeCC": 0.0005,
        "yFeIa": 0.0007,
        "aFeCC": 0.45,
        "mu": 1.1,
        "fRetCC": 1.0,
        "r": 0.4,
        "eta": 0.3,
        "tauIa": 1.5,
        "tDminIa": 0.05,
        "SolarO": 0.0073,
        "SolarMg": 0.0007,
        "SolarFe": 0.0014,
        "IaDTD_fn": "exponential",
    }
