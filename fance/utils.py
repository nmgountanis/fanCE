import numpy as np


def stellar_props(
        t: np.ndarray,
        t0: float,
        sfr: np.ndarray,
        OH: np.ndarray,
        MgH: np.ndarray,
        FeH: np.ndarray,
        light_weighted: bool = True,
):
    """
    Compute light-/mass-weighted, log-averaged age and abundances.

    :param np.ndarray t: Model output (not input!) time array (in Gyr).
    :param float t0: Age of Universe at observed redshift [14]
    :param np.ndarray sfr: Star formation rate at times in t
    :param np.ndarray OH: Stellar [O/H] at times in t
    :param np.ndarray MgH: Stellar [Mg/H] at times in t
    :param np.ndarray FeH: Stellar [Fe/H] at times in t
    :param light_weighted: Use light weighting from Gountanis+24 Eq. 1.
        If false, outputs are mass-weighted.
    :return: Tuple[float, float, float, float, float, float]:
        mean_age, OHStar, MgHStar, FeHStar, OFeStar, MgFeStar --- each a 1-d array corresponding to the time array
    """

    if light_weighted:  # light-weighting based off Gountanis+24 Eq. 1
        weights = 1.5 * ((t0 - t + 0.1)**-0.9)
    else:  # mass weighted
        weights = np.ones_like(t)
    OHStar = np.cumsum(OH * sfr * weights) / np.cumsum(weights * sfr)
    MgHStar = np.cumsum(MgH * sfr * weights) / np.cumsum(weights * sfr)
    FeHStar = np.cumsum(FeH * sfr * weights) / np.cumsum(weights * sfr)
    OFeStar = OHStar - FeHStar
    MgFeStar = MgHStar - FeHStar
    # light-/mass-weighted, log-averaged age
    mean_age = 10**(np.cumsum(np.log10(t0 - t + 1e-5) * sfr * weights) / np.cumsum(sfr * weights))
    return mean_age, OHStar, MgHStar, FeHStar, OFeStar, MgFeStar


def add_burst():
    raise NotImplementedError
