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
    Plots normalized star formation rate versus time
    
    Input arguments:
      t = 1-d time array
      sfr = star formation rate, 1-d array evaluated at the times in t
      labels = string label
    
    Returns:
     A figure ploting the star formation history normalized to give the 
     galaxy's stellar mass formed per Gyr.
      
     """
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
    plt.style.use('seaborn-v0_8-colorblind')

    linew = 3
    
    if axes==None: axes = plt.gca()
        
    ax.plot(t,sfr/dtout,color="C2",label=labels, linewidth=linew)
    
    ax.set_xlabel('time [Gyr]')
    ax.set_ylabel('SFR')
    ax.set_xlim([-0.5, 14.5])
    ax.set_ylim([0, .17])
    ax.yaxis.set_ticks([0,.05,.1, .15])
    ax.legend(fontsize="18", loc = "upper right")
    
    plt.tight_layout()
    
def plotGas(t,MgH,FeH,MgFe,labels,axes=None):
    
    """
    Plots [Mg/H], [Fe/H], and [Mg/Fe] versus time.
    [O/H] and [O/Fe] are identical to [Mg/H] and [Mg/Fe].
    
    Input arguments:
      t = 1-d time array
      MgH = [Mg/H], 1-d array evaluated at the times in t
      FeH = [Fe/H], 1-d array evaluated at the times in t
      MFe = [Mg/Fe], 1-d array evaluated at the times in t
      labels = string label
    
    Returns:
     A figure ploting [Mg/H], [Fe/H], and [Mg/Fe] versus time
     in seperate subfigures.
     
    """
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))
    plt.style.use('seaborn-v0_8-colorblind')

    linew = 3
    
    if axes==None: axes = plt.gca()
    
    ax[0].plot(t,MgH,color="C2", linewidth=linew)
    ax[0].set_xlabel('time [Gyr]')
    ax[0].set_ylabel('[Mg/H]')
    ax[0].set_xlim([-0.5, 14.5])
    
    ax[1].plot(t,FeH,color="C2", linewidth=linew)
    ax[1].set_xlabel('time [Gyr]')
    ax[1].set_ylabel('[Fe/H]')
    ax[1].set_xlim([-0.5, 14.5])

    ax[2].plot(t,MgFe,color="C2",label=labels,linewidth=linew)
    ax[2].set_ylabel('[Mg/Fe]')
    ax[2].set_xlabel('time [Gyr]')
    ax[2].set_xlim([-0.5, 14.5])
    ax[2].legend(fontsize="18", loc = "upper right")
    
    plt.tight_layout()
    