import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

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
mpl.rcParams["legend.title_fontsize"] = 18
mpl.rcParams["errorbar.capsize"] = 3


def plotSFH(
        t: np.ndarray,
        sfr: np.ndarray,
        labels: str,
        linewidth: float = 3,
):
    """
    Plots normalized star formation rate versus time

    :param np.ndarray t: Model output (not input!) time array (in Gyr).
    :param np.ndarray sfr: Star formation rate at times in t
    :param str labels: Model label
    :param float linewidth: Linewidth of plot
    :return: None
    """
    dt = t[1] - t[0]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
    plt.style.use("seaborn-v0_8-colorblind")
    ax.plot(t, sfr / dt, color="C2", label=labels, linewidth=linewidth)
    ax.set_xlabel("time [Gyr]")
    ax.set_ylabel("SFR")
    ax.set_xlim([-0.5, 14.5])
    ax.set_ylim([0, 0.17])
    ax.yaxis.set_ticks([0, 0.05, 0.1, 0.15])
    ax.legend(fontsize="18", loc="upper right")
    plt.tight_layout()
    plt.show()

def plotGas(
        t: np.ndarray,
        MgH: np.ndarray,
        FeH: np.ndarray,
        MgFe: np.ndarray,
        labels: str,
        linewidth: float = 3,
):
    """
    Plots [Mg/H], [Fe/H], and [Mg/Fe] versus time in seperate subfigures.
    [O/H] and [O/Fe] are identical to [Mg/H] and [Mg/Fe].

    :param np.ndarray t: Model output (not input!) time array (in Gyr).
    :param np.ndarray sfr: Star formation rate at times in t
    :param np.ndarray MgH: [Mg/H], 1-d array evaluated at the times in t
    :param np.ndarray FeH: [Fe/H], 1-d array evaluated at the times in t
    :param np.ndarray MgFe: [Mg/Fe], 1-d array evaluated at the times in t
    :param str labels: Model label
    :param float linewidth: Linewidth of plot
    :return: None
    """
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))
    plt.style.use("seaborn-v0_8-colorblind")

    ax[0].plot(t, MgH, color="C2", linewidth=linewidth)
    ax[0].set_xlabel("time [Gyr]")
    ax[0].set_ylabel("[Mg/H]")
    ax[0].set_xlim([-0.5, 14.5])

    ax[1].plot(t, FeH, color="C2", linewidth=linewidth)
    ax[1].set_xlabel("time [Gyr]")
    ax[1].set_ylabel("[Fe/H]")
    ax[1].set_xlim([-0.5, 14.5])

    ax[2].plot(t, MgFe, color="C2", label=labels, linewidth=linewidth)
    ax[2].set_ylabel("[Mg/Fe]")
    ax[2].set_xlabel("time [Gyr]")
    ax[2].set_xlim([-0.5, 14.5])
    ax[2].legend(fontsize="18", loc="upper right")

    plt.tight_layout()
    plt.show()


def plotStellar(
        t: np.ndarray,
        sfr: np.ndarray,
        mean_age: np.ndarray,
        MgH: np.ndarray,
        MgFe: np.ndarray,
        MgHStar: np.ndarray,
        MgFeStar: np.ndarray,
        labels: str,
        linewidth: float = 3,
):
    """
    Plots normalized SFH, <age>, [Mg/H], <[Mg/H]>, [Mg/Fe]
    and <[Mg/Fe]> versus time and with <[Mg/H]> and <[Mg/Fe]>
    versus <age>.

    :param np.ndarray t: Model output (not input!) time array (in Gyr).
    :param np.ndarray sfr: Star formation rate at times in t
    :param np.ndarray mean_age: <age>, 1-d array evaluated at the times in t
    :param np.ndarray MgH: [Mg/H], 1-d array evaluated at the times in t
    :param np.ndarray MgFe: [Mg/Fe], 1-d array evaluated at the times in t
    :param np.ndarray MgHStar: <[Mg/H]>, 1-d array evaluated at the times in t
    :param np.ndarray MgFeStar: <[Mg/Fe]>, 1-d array evaluated at the times in t
    :param str labels: Model label
    :param float linewidth: Linewidth of plot
    :return: None
    """
    dt = t[1] - t[0]
    fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(11, 18))
    plt.style.use("seaborn-v0_8-colorblind")

    ax[0, 0].plot(t, sfr / dt, c="C2", linewidth=linewidth)
    ax[0, 0].set_xlabel("time (Gyr)", labelpad=-2)
    ax[0, 0].set_ylabel("SFR")
    ax[0, 0].set_xlim([-0.5, 14.5])

    ax[0, 1].plot(t, mean_age, c="C2", label=labels, linewidth=linewidth)
    ax[0, 1].set_xlabel("time (Gyr)", labelpad=-2)
    ax[0, 1].set_ylabel("$\langle$age$\\rangle$ (Gyr)")
    ax[0, 1].set_xlim([-0.5, 14.5])
    ax[0, 1].yaxis.set_ticks([5, 10, 15])
    ax[0, 1].legend(fontsize="16", loc="upper right")

    ax[1, 0].plot(t, MgH, c="C2", linewidth=linewidth)
    ax[1, 0].set_xlabel("time (Gyr)", labelpad=-2)
    ax[1, 0].set_ylabel("[Mg/H]")
    ax[1, 0].set_xlim([-0.5, 14.5])
    ax[1, 0].set_ylim([-2.5, 0.75])

    ax[1, 1].plot(t, MgHStar, c="C2", linewidth=linewidth)
    ax[1, 1].set_xlabel("time (Gyr)", labelpad=-2)
    ax[1, 1].set_ylabel("$\langle$[Mg/H]$\\rangle$")
    ax[1, 1].set_xlim([-0.5, 14.5])
    ax[1, 1].set_ylim([-2.5, 0.75])

    ax[2, 0].plot(t, MgFe, c="C2", linewidth=linewidth)
    ax[2, 0].set_xlabel("time (Gyr)", labelpad=-2)
    ax[2, 0].set_ylabel("[Mg/Fe]", labelpad=-15)
    ax[2, 0].set_xlim([-0.5, 14.5])
    ax[2, 0].set_ylim([-0.2, 0.575])

    ax[2, 1].plot(t, MgFeStar, c="C2", linewidth=linewidth)
    ax[2, 1].set_xlabel("time (Gyr)", labelpad=-2)
    ax[2, 1].set_ylabel("$\langle$[Mg/Fe]$\\rangle$", labelpad=-15)
    ax[2, 1].set_xlim([-0.5, 14.5])
    ax[2, 1].set_ylim([-0.2, 0.575])

    ax[3, 0].plot(mean_age, MgHStar, c="C2", linewidth=linewidth)
    ax[3, 0].set_xlabel("$\langle$age$\\rangle$ (Gyr)", labelpad=-1)
    ax[3, 0].set_ylabel("$\langle$[Mg/H]$\\rangle$")
    ax[3, 0].set_xlim([1, 14.5])
    ax[3, 0].set_ylim([-2.5, 0.75])

    ax[3, 1].plot(mean_age, MgFeStar, c="C2", linewidth=linewidth)
    ax[3, 1].set_xlabel("$\langle$age$\\rangle$ (Gyr)", labelpad=-1)
    ax[3, 1].set_ylabel("$\langle$[Mg/Fe]$\\rangle$", labelpad=-15)
    ax[3, 1].set_xlim([1, 14.5])
    ax[3, 1].set_ylim([-0.2, 0.5])

    plt.tight_layout()
    fig.subplots_adjust(hspace=0.28)
    fig.subplots_adjust(wspace=0.3)
    plt.show()
