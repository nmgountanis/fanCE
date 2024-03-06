'''
fanCE default parameters
These parameters describe the reference model in Gountanis et al. 2024
All timescales are given in Gyr
'''

tstart=0.5              # star formation start time
t0=14.0                 # maximum time
dtout=0.02              # output timestep

tauStar=1               # star formation efficiency timescale
eta=0.3                 # outflow mass loading factor
r=0.4                   # recyling fraction
fret=1.0                # CCSN retention factor (1-fret is direct loss)

'''
The SFR is proportional to (1-e^{-t/tau1})*e^{-t/tau2}.
For typical choices, this rises linearly on a timescale tau1 and then falls
exponentially on a timescale tau2.
tau1=0 is an acceptable limit.
tau2=0 is interpreted as tau2=infinity.
tau1=tau2=0 is a constant SFR. 
'''
tau1=2                  # rise timescale
tau2=8                  # fall timescale

'''
SNIa DTD parameters
tauIa=0 is interpreted as a double exponential instead of a single. 
If tauIa=0, MUST choose tdmin=0.05 or 0.15.
'''
tauIa=0                 # e-folding time for SNIa exponential DTD
tdmin=0.15              # minimum time delay for SNIa

'''
Solar abundances by mass
Based on Magg et al. (2022) Table 5 + 0.04 dex to correct for diffusion
'''
SolarO=7.33e-3
SolarMg=6.71e-4
SolarFe=1.37e-3

'''
IMF-averaged CCSN yields
The yield calibration is based on Weinberg++ 2023, eq. 11
'''
afecc=0.45              # plateau value for [alpha/Fe]

mocc=0.973*SolarO*(0.00137/SolarFe)*(10**(afecc-0.45))  # CCSN oxygen
mfecc=mocc*(SolarFe/SolarO)*(10**(-afecc))              # CCSN iron
mmgcc=mocc*SolarMg/SolarO                               # CCSN magnesium

'''
Population averaged SNIa Fe yield
Will evolve to afeeq at late times when integrated to t=infity for a constant SFR
'''
afeeq=0.05
mfeIa=mfecc*(10.**(afecc-afeeq) - 1.)
