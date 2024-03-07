import matplotlib as mpl
import matplotlib.pyplot as plt
import fanCE_par as p
from fanCE import *

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

def plotSFH(t,sfr,dtout,labels,axes=None):
    
    """
    Plots the observed velocity function (no completeness corrections).
    """
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    plt.style.use('seaborn-v0_8-colorblind')

    linew = 3
    
    if axes==None: axes = plt.gca()
        
    ax.plot(t,sfr/dtout,color="C0",label=labels, linewidth=linew)
    
    ax.set_ylabel('SFR')
    ax.set_xlabel('time [Gyr]')
    ax.set_xlim([-0.5, 14.5])
    ax.set_ylim([0, .17])
    ax.yaxis.set_ticks([0,.05,.1, .15])
    ax.legend(title= "$\\tau_1$, $\\tau_2$, t$_{\mathrm{start}}$", fontsize="18", loc = "upper right")
    plt.tight_layout()