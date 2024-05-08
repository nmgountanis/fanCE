import numpy as np


def waf2017(
    tmod: np.ndarray,
    tauSFE: float = 1.0,
    tauSFH: float = 8.0,
    yOCC: float = 0.0071,
    yMgCC: float = 0.0007,
    yFeCC: float = 0.0005,
    yFeIa: float = 0.0007,
    r: float = 0.4,
    eta: float = 0.3,
    tauIa: float = 1.5,
    tDminIa: float = 0.15,
    SolarO: float = 0.0073,
    SolarMg: float = 0.0007,
    SolarFe: float = 0.0014,
    SFH_fn: str = "exponential",
    IaDTD_fn: str = "powerlaw",
):
    """
    Compute [O/H], [Mg/H], [Fe/H], [O/Fe] and [Mg/Fe] tracks using WAF analytic model for a constant,
    exponentially declining, or linear-exponential SFR.
    Reference: Weinberg, Andrews, & Freudenburg 2017, ApJ 837, 183 (Particularly Appendix C)

    :param np.ndarray tmod: Input time array (in Gyr) for model.
        Should be output time array - star formation start time.
    :param float tauSFE: SFE efficiency timescale [1.0]
    :param float tauSFH: e-folding timescale of SFH [8.0]
        (Ignored if SFH_fn == 'constant')
    :param float yOCC: IMF-averaged CCSN O yield [0.0071]
    :param float yMgCC: IMF-averaged CCSN Mg yield [0.0007]
    :param float yFeCC: IMF-averaged CCSN iron yield [0.0005]
    :param float yFeIa: IMF-averaged SNIa iron yield over 10 Gyr [0.0007]
    :param float r: recycling fraction [0.4, based on Kroupa IMF]
    :param float eta: outflow mass loading factor [0.3]
    :param float tauIa: = e-folding time for SNIa DTD [1.5]
        (Ignored if IaDTD_fn == 'powerlaw')
    :param float tDminIa: minimum time delay for SNIa [0.15]
        (If IaDTD_fn == 'powerlaw', then tDminIa must be either 0.05 or 0.15)
    :param float SolarO: Solar O mass fraction [0.0073]
    :param float SolarMg: Solar Mg mass fraction [0.0007]
    :param float SolarFe: Solar iron mass fraction [0.0014]
    :param str SFH_fn: functional form of the SFH. ['exponential']
        Must be one of 'constant', 'exponential', or 'linexp'
    :param IaDTD_fn: functional form of the SN Ia DTD ['powerlaw']
        Must be either 'exponential' or 'powerlaw'
    :return Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        SFR, OH, MgH, FeH, OFe, MgFe --- each a 1-d array corresponding to the time array, tmod
        (SFR is normalized such that for evenly spaced tmod, it is the fraction of stars formed in each time step.)
    """
    # Prophylactically preventing divide-by-zero errors
    xfac = 1.0013478
    yfac = 1.0016523
    # Parse Ia DTD
    if IaDTD_fn == "exponential":
        tauIa1 = tauIa
        tauIa2 = tauIa
        Ia_norm1 = 1.0
        Ia_norm2 = 0.0
    elif IaDTD_fn == "powerlaw":
        if tDminIa == 0.05:
            tauIa1 = 0.25
            tauIa2 = 3.5
            Ia_norm1 = 0.493
            Ia_norm2 = 0.507
        elif tDminIa == 0.15:
            tauIa1 = 0.5
            tauIa2 = 5.0
            Ia_norm1 = 0.478
            Ia_norm2 = 0.522
        else:
            RuntimeError(
                "tDminIa must be either 0.05 or 0.15 if IaDTD_fn == 'powerlaw' "
            )
    else:
        raise RuntimeError("IaDTD_fn must be either 'exponential' or 'powerlaw' ")
    # Compute "Harmonic Difference" Timescales
    tauSFE *= xfac
    tauSFH *= yfac
    tauDep = tauSFE / (1 + eta - r)
    tauDepIa1 = 1.0 / (1.0 / tauDep - 1.0 / tauIa1)
    tauDepIa2 = 1.0 / (1.0 / tauDep - 1.0 / tauIa2)
    if SFH_fn == "constant":
        tauDepSFH = tauDep
        tauIaSFH1 = tauIa1
        tauIaSFH2 = tauIa2
        tauSFH = 1e8
    elif SFH_fn in ["exponential", "linexp"]:
        tauDepSFH = 1.0 / (1.0 / tauDep - 1.0 / tauSFH)
        tauIaSFH1 = 1.0 / (1.0 / tauIa1 - 1.0 / tauSFH)
        tauIaSFH2 = 1.0 / (1.0 / tauIa2 - 1.0 / tauSFH)
    else:
        RuntimeError("SFH_fn must be one of 'constant', 'exponential', or 'linexp' ")
    # Compute equilibrium abundances, WAF equations 28-30
    ZOEq = yOCC * tauDepSFH / tauSFE
    ZMgEq = yMgCC * tauDepSFH / tauSFE
    ZFeCCEq = yFeCC * tauDepSFH / tauSFE
    ZFeIaEq1 = (
        Ia_norm1
        * yFeIa
        * ((tauDepSFH / tauSFE) * (tauIaSFH1 / tauIa1) * np.exp(tDminIa / tauSFH))
    )
    ZFeIaEq2 = (
        Ia_norm2
        * yFeIa
        * ((tauDepSFH / tauSFE) * (tauIaSFH2 / tauIa2) * np.exp(tDminIa / tauSFH))
    )
    # Compute non-equilibrium abundances
    if SFH_fn in ["constant", "exponential"]:  # WAF equations 50, 52, and 53
        delta_t = tmod - tDminIa
        ZO = ZOEq * (1.0 - np.exp(-tmod / tauDepSFH))
        ZMg = ZMgEq * (1.0 - np.exp(-tmod / tauDepSFH))
        ZFeCC = ZO * ZFeCCEq / ZOEq
        ZFeIa1 = ZFeIaEq1 * (
            1.0
            - np.exp(-delta_t / tauDepSFH)
            - (tauDepIa1 / tauDepSFH)
            * (np.exp(-delta_t / tauIaSFH1) - np.exp(-delta_t / tauDepSFH))
        )
        ZFeIa2 = ZFeIaEq2 * (
            1.0
            - np.exp(-delta_t / tauDepSFH)
            - (tauDepIa2 / tauDepSFH)
            * (np.exp(-delta_t / tauIaSFH2) - np.exp(-delta_t / tauDepSFH))
        )
    elif SFH_fn == "linexp":  # WAF equations 56-58
        delta_t = tmod - tDminIa
        ZO = ZOEq * (1.0 - (tauDepSFH / tmod) * (1.0 - np.exp(-tmod / tauDepSFH)))
        ZMg = ZMgEq * (1.0 - (tauDepSFH / tmod) * (1.0 - np.exp(-tmod / tauDepSFH)))
        ZFeCC = ZO * ZFeCCEq / ZOEq
        ZFeIa1 = (
            ZFeIaEq1
            * (tauIaSFH1 / tmod)
            * (
                delta_t / tauIaSFH1
                + (tauDepIa1 / tauDepSFH) * np.exp(-delta_t / tauIaSFH1)
                + (1.0 + (tauDepSFH / tauIaSFH1) - (tauDepIa1 / tauDepSFH))
                * np.exp(-delta_t / tauDepSFH)
                - (1.0 + (tauDepSFH / tauIaSFH1))
            )
        )
        ZFeIa2 = (
            ZFeIaEq2
            * (tauIaSFH2 / tmod)
            * (
                delta_t / tauIaSFH2
                + (tauDepIa2 / tauDepSFH) * np.exp(-delta_t / tauIaSFH2)
                + (1.0 + (tauDepSFH / tauIaSFH2) - (tauDepIa2 / tauDepSFH))
                * np.exp(-delta_t / tauDepSFH)
                - (1.0 + (tauDepSFH / tauIaSFH2))
            )
        )
    else:
        RuntimeError("SFH_fn must be one of 'constant', 'exponential', or 'linexp' ")
    ZFeIa1[tmod < tDminIa] = 0
    ZFeIa2[tmod < tDminIa] = 0
    ZFe = ZFeCC + ZFeIa1 + ZFeIa2
    # Compute SFR
    if SFH_fn == "constant":
        SFR = np.ones_like(tmod)
    elif SFH_fn == "exponential":
        SFR = np.exp(-tmod / tauSFH)
    elif SFH_fn == "linexp":
        SFR = tmod * np.exp(-tmod / tauSFH)
    else:
        RuntimeError("SFH_fn must be one of 'constant', 'exponential', or 'linexp' ")
    SFR /= np.sum(SFR)
    # Convert to [Alpha/H], [Fe/H]
    OH = np.log10(ZO / SolarO)
    MgH = np.log10(ZMg / SolarMg)
    FeH = np.log10(ZFe / SolarFe)
    OFe = OH - FeH
    MgFe = MgH - FeH
    return SFR, OH, MgH, FeH, OFe, MgFe


def fance(
    aFeCC: float = 0.45,
    tauSFE: float = 1.0,
    tauSFH1: float = 2.0,
    tauSFH2: float = 8.0,
    yOCC: float = 0.0071,
    yMgCC: float = 0.0007,
    yFeCC: float = 0.0005,
    yFeIa: float = 0.0007,
    fRetCC: float = 1.0,
    r: float = 0.4,
    eta: float = 0.3,
    tauIa: float = 1.5,
    tDminIa: float = 0.05,
    SolarO: float = 0.0073,
    SolarMg: float = 0.0007,
    SolarFe: float = 0.0014,
    IaDTD_fn: str = "exponential",
    t_start: float = 0.5,
    t0: float = 14.0,
    dt: float = 0.02,
):
    
    """
    Compute [O/H], [Mg/H], [Fe/H], [O/Fe] and [Mg/Fe] tracks using WAF analytic model for a rise-fall SFR.
    Reference: Weinberg, Andrews, & Freudenburg 2017, ApJ 837, 183 (Particularly Appendix C)

    :param np.ndarray tmod: Input time array (in Gyr) for model.
        Should be output time array - star formation start time.
    :param float tauSFE: SFE efficiency timescale [1.0]
    :param float tauSFH1: Along with tauSFH2, determines the SFR of the form
        SFR \propto (1-exp(-tmod/tau1))*exp(-tau2) [2.0]
        (If tauSFH1 <= 0, then the model reverts to an exponentially declining SFH with tauSFH = tauSFH2)
    :param float tauSFH2: Along with tauSFH1, determines the SFR of the form
        SFR \propto (1-exp(-tmod/tau1))*exp(-tau2) [8.0]
        (If tauSFH1 <= 0 and tauSFH2 <= 0, then the model reverts to a constant SFH)
    :param float yOCC: IMF-averaged CCSN O yield [0.0071]
    :param float yMgCC: IMF-averaged CCSN Mg yield [0.0007]
    :param float yFeCC: IMF-averaged CCSN iron yield [0.0005]
    :param float yFeIa: IMF-averaged SNIa iron yield over 10 Gyr [0.0007]
    :param float fRetCC: Fraction of metals returned to the ISM by CCSNe [1.0]
    :param float r: recycling fraction [0.4, based on Kroupa IMF]
    :param float eta: outflow mass loading factor [0.3]
    :param float tauIa: = e-folding time for SNIa DTD [1.5]
        (Ignored if IaDTD_fn == 'powerlaw')
    :param float tDminIa: minimum time delay for SNIa [0.15]
        (If IaDTD_fn == 'powerlaw', then tDminIa must be either 0.05 or 0.15)
    :param float SolarO: Solar O mass fraction [0.0073]
    :param float SolarMg: Solar Mg mass fraction [0.0007]
    :param float SolarFe: Solar iron mass fraction [0.0014]
    :param IaDTD_fn: functional form of the SN Ia DTD ['powerlaw']
        Must be either 'exponential' or 'powerlaw'
    :return Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        SFR, OH, MgH, FeH, OFe, MgFe --- each a 1-d array corresponding to the time array, tmod
        (SFR is normalized such that for evenly spaced tmod, it is the fraction of stars formed in each time step.)
    """

    # define intermediate time array for normalization reasons. 
    # Cannot start at zero. Also must go to t0-t_start
    tmod = np.arange(dt, t0 - t_start + dt, dt)
    
    # Define t used for plotting
    t = tmod + t_start

    # Set the Yields (CC and Ia) by CC plateu defined by user (aFeCC)
    # Followig eq. 17 (?) in Weinberg++23
    SolarO = 0.00733
    SolarMg = 0.000671
    SolarFe = 0.00137
    aFeEq = 0.0 # late time evolution equilibrium [a/fe]
    yOCC = 0.973 * SolarO * (0.00137 / SolarFe) * (10**(aFeCC-0.45))
    yFeCC = yOCC * (SolarFe / SolarO) * (10**(-aFeCC))
    yMgCC = yOCC * SolarMg / SolarO

    delta = aFeCC - aFeEq
    mu = 1.1
    yFeIa = yFeCC * (10.**(delta) - 1.) / mu
    
    
    # Modulate Yields by Return Fraction
    yOCC *= fRetCC
    yMgCC *= fRetCC
    yFeCC *= fRetCC
    # Parse SFH Parameters & Run WAF models
    if (tauSFH1 <= 0) and (tauSFH2 <= 0):
        SFR, OH, MgH, FeH, OFe, MgFe = waf2017(
            tmod=tmod,
            tauSFE=tauSFE,
            yOCC=yOCC,
            yMgCC=yMgCC,
            yFeCC=yFeCC,
            yFeIa=yFeIa,
            r=r,
            eta=eta,
            tauIa=tauIa,
            tDminIa=tDminIa,
            SolarO=SolarO,
            SolarMg=SolarMg,
            SolarFe=SolarFe,
            SFH_fn="constant",
            IaDTD_fn=IaDTD_fn,
        )
    elif (tauSFH1 <= 0) and (tauSFH2 > 0):
        SFR, OH, MgH, FeH, OFe, MgFe = waf2017(
            tmod=tmod,
            tauSFE=tauSFE,
            tauSFH=tauSFH2,
            yOCC=yOCC,
            yMgCC=yMgCC,
            yFeCC=yFeCC,
            yFeIa=yFeIa,
            r=r,
            eta=eta,
            tauIa=tauIa,
            tDminIa=tDminIa,
            SolarO=SolarO,
            SolarMg=SolarMg,
            SolarFe=SolarFe,
            SFH_fn="exponential",
            IaDTD_fn=IaDTD_fn,
        )
    else:
        # SFR = constant*[exp(-tmod/tau2)-exp(-tmod/tauh)], with tauh = (1/tau1+1/tau2)
        # Evaluate each exponential term with waf2017() and combine, follow WAF Eq. 117
        if tauSFH2 <= 0:
            tauh = tauSFH1
            SFR2, OH2, MgH2, FeH2, OFe2, MgFe2 = waf2017(
                tmod=tmod,
                tauSFE=tauSFE,
                yOCC=yOCC,
                yMgCC=yMgCC,
                yFeCC=yFeCC,
                yFeIa=yFeIa,
                r=r,
                eta=eta,
                tauIa=tauIa,
                tDminIa=tDminIa,
                SolarO=SolarO,
                SolarMg=SolarMg,
                SolarFe=SolarFe,
                SFH_fn="constant",
                IaDTD_fn=IaDTD_fn,
            )
            SFR2 = np.ones_like(SFR2)
        else:
            tauh = 1 / (1 / tauSFH1 + 1 / tauSFH2)
            SFR2, OH2, MgH2, FeH2, OFe2, MgFe2 = waf2017(
                tmod=tmod,
                tauSFE=tauSFE,
                tauSFH=tauSFH2,
                yOCC=yOCC,
                yMgCC=yMgCC,
                yFeCC=yFeCC,
                yFeIa=yFeIa,
                r=r,
                eta=eta,
                tauIa=tauIa,
                tDminIa=tDminIa,
                SolarO=SolarO,
                SolarMg=SolarMg,
                SolarFe=SolarFe,
                SFH_fn="exponential",
                IaDTD_fn=IaDTD_fn,
            )
            SFR2 = np.exp(-tmod / tauSFH2)
        SFRh, OHh, MgHh, FeHh, OFeh, MgFeh = waf2017(
            tmod=tmod,
            tauSFE=tauSFE,
            tauSFH=tauh,
            yOCC=yOCC,
            yMgCC=yMgCC,
            yFeCC=yFeCC,
            yFeIa=yFeIa,
            r=r,
            eta=eta,
            tauIa=tauIa,
            tDminIa=tDminIa,
            SolarO=SolarO,
            SolarMg=SolarMg,
            SolarFe=SolarFe,
            SFH_fn="exponential",
            IaDTD_fn=IaDTD_fn,
        )
        SFRh = -1 * np.exp(-tmod / tauh)
        SFR = SFR2 + SFRh
        # Convert to Linear Metallicity and Combine with WAF Eq. 117
        ZO2 = SolarO * 10 ** OH2
        ZMg2 = SolarMg * 10 ** MgH2
        ZFe2 = SolarFe * 10 ** FeH2
        ZOh = SolarO * 10 ** OHh
        ZMgh = SolarMg * 10 ** MgHh
        ZFeh = SolarFe * 10 ** FeHh
        ZO = (SFR2 * ZO2 + SFRh * ZOh) / SFR
        ZMg = (SFR2 * ZMg2 + SFRh * ZMgh) / SFR
        ZFe = (SFR2 * ZFe2 + SFRh * ZFeh) / SFR
        # Convert Back to Log Metallicity
        OH = np.log10(ZO / SolarO)
        MgH = np.log10(ZMg / SolarMg)
        FeH = np.log10(ZFe / SolarFe)
        OFe = OH - FeH
        MgFe = MgH - FeH
    # Normalize SFR
    SFR /= np.sum(SFR)
    return t, SFR, OH, MgH, FeH, OFe, MgFe

