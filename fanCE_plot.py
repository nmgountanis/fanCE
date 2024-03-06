import matplotlib as mpl
import matplotlib.pyplot as plt
import fance_par as p
from fance import *

mpl.rcParams["font.family"] = "serif"
mpl.rcParams["text.usetex"] = True
mpl.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
mpl.rcParams["errorbar.capsize"] = 5
mpl.rcParams["axes.linewidth"] = 2
mpl.rcParams["xtick.major.size"] = 16
mpl.rcParams["xtick.major.width"] = 2
mpl.rcParams["xtick.minor.size"] = 8
mpl.rcParams["xtick.minor.width"] = 1
mpl.rcParams["ytick.major.size"] = 16
mpl.rcParams["ytick.major.width"] = 2
mpl.rcParams["ytick.minor.size"] = 8
mpl.rcParams["ytick.minor.width"] = 1
mpl.rcParams["axes.labelsize"] = 28
mpl.rcParams["xtick.labelsize"] = 25
mpl.rcParams["ytick.labelsize"] = 25
mpl.rcParams["legend.fontsize"] = 25
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["ytick.right"] = True
mpl.rcParams["xtick.top"] = True
mpl.rcParams["xtick.minor.visible"] = True
mpl.rcParams["ytick.minor.visible"] = True
mpl.rcParams['legend.title_fontsize'] = 18
mpl.rcParams["errorbar.capsize"] = 3

def plotSFR(t,sfr,label=True, axes=None):
    
    """
    Plots the observed velocity function (no completeness corrections).
    """
    if axes==None: axes = plt.gca()
    obssigmas = np.array([ dwarf['sigma'] for dwarf_name,dwarf in dwarfs.items() if dwarf['sigma'] != None ])
    axes.plot(t,obsvfxn,color='0.5',label='observed' if label==True else '') # 'k' if plotband else 'C5' # obs

    
    