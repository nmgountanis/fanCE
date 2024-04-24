"""
Functions for analytic one-zone chemical evolution.

The first, waf, implements the analytic solution of Weinberg, Andrews, & Freudenburg 
 (2017) for a constant or exponential declinining star formation history.

The second, waf_linexp, implements the WAF solution for a linear-exponential
(a.k.a. delayed tau) star formation history.

The third, fance (flexible analytic chemical evolution), implements a two-parameter
star formation history proportional to (1-exp(-t/tau1))*exp(-t/tau2). For tau2>tau1>0, 
this form rises linearly at early times, reaches a maximum at 
t = tau1*ln(1+tau2/tau1) ~ tau1, and falls exponentialy with a timescale tau2 at late times.

The solution for this case can be implemented using equation (117) of WAF and 
calls to the waf function.

fance also allows for CCSN-enhanced winds by incorporating a retention factor,
fret, for CCSN ejecta.  If fret<1, winds are assumed to be comprised of a 
fraction (1-fret) of the CCSN ejecta that entrain sufficient ambient ISM
to raise the overall mass-loading factor to eta.  If fret=1 this just becomes
the usual case in which the galactic wind has the ambient ISM metallicity.
It can be shown that choosing fret<1 is identical to reducing the CCSN yield
by a factor fret, with the same eta.

"""

def waf(tstart,t0,dtout,mocc,mmgcc,mfecc,mfeIa,r,SolarO,SolarMg,SolarFe,tauIa,tdmin,eta,tauStar,tauSfh):
    """
    Compute [O/H], [Mg/H], [Fe/H], [O/Fe], and [Mg/Fe] tracks using WAF analytic model for
    constant SFR or exponentially declining SFR

    Input arguments (illustrative values in [  ]):
      tstart = star formation start time [0.5]
      t0 = age of universe at observed redshift [14]
      dtout = output timestep [0.02]
      mocc = IMF-averaged CCSN oxygen yield [0.0071]
      mmgcc = IMF-averaged CCSN magnessium yield [0.0007] 
      mfecc = IMF-averaged CCSN iron yield [0.0005]
      mfeIa = IMF-averaged SNIa iron yield over 10 Gyr [0.0007]
      r = recycling fraction [0.4, based on Kroupa IMF]
      SolarO = solar oxygen mass fraction [0.0073]
      SolarMg = solar magnessium mass fraction [0.0007]
      SolarFe = solar iron mass fraction [0.0014]
      tauIa = e-folding time for SNIa DTD [1.5]
              tauIa <= 0, will use double-exp approximation to power-law
      tdmin = minimum time delay for SNIa [0.05]
              for tauIa <= 0 MUST choose tdmin = 0.05 or 0.15
      eta = outflow mass loading factor [0.3]
      tauStar = SFE efficiency timescale [1.0]
      tauSfh = e-folding timescale of SFH [8.0] 
               if tauSfh <= 0, constant SFR will be used

    Returns:
      t, SFR, OH, MgH, FeH, OFe, and MgFe
      each a 1-d array evaluated at the times in t
      SFR is normalized to sum to 1.0, so for evenly spaced t it is the
      fraction of stars formed in each bin

    Time unit is Gyr
    Reference: Weinberg, Andrews, & Freudenburg 2017, ApJ 837, 183
               Particularly Appendix C
    
    """

    import numpy as np
    
    tmax = t0 - tstart                        # maximum time
    
    tmod = np.arange(dtout, tmax+1.e-5,dtout) # fance() input time 
    
    t = tmod + tstart                         # actual time

    # these near-unit multiplications will be used to ensure that timescales
    # don't actually become equal to each other causing divide-by-zero errors
    xfac=1.0013478
    yfac=1.0016523

    # For tauIa=0, use double-exponential that approximates t^{-1.1} power-law
    # see WAF section 5.1 and Appx C for details
    # In this case, we require tdmin=0.05 or tdmin=0.15
    if (tauIa<0.01):
        if (tdmin==0.15):
            tauIa1=0.5
            tauIa2=5.0
            norm1=0.478
            norm2=0.522
        elif (tdmin==0.05):
            tdmin=0.05
            tauIa1=0.25
            tauIa2=3.5
            norm1=0.493
            norm2=0.507
        else:
            print ("can't set DTD parameters for this tdmin\n")
            sys.exit()
    # otherwise set normalization of second exponential to zero
    else:
        tauIa1=tauIa
        tauIa2=tauIa
        norm1=1.0
        norm2=0.0

    # compute harmonic difference timescales
    tauStar*=xfac
    tauSfh*=yfac
    tauDep=tauStar/(1+eta-r)
    tauDepIa1=1./(1./tauDep-1./tauIa1)
    tauDepIa2=1./(1./tauDep-1./tauIa2)
    if (tauSfh>0.0):
        tauDepSfh=1./(1./tauDep-1./tauSfh)
        tauIaSfh1=1./(1./tauIa1-1./tauSfh)
        tauIaSfh2=1./(1./tauIa2-1./tauSfh)
    else:
        tauDepSfh=tauDep
        tauIaSfh1=tauIa1
        tauIaSfh2=tauIa2
        tauSfh=1.e8

    # compute equilibrium abundances, WAF equations 28-30
    ZoEq=mocc*tauDepSfh/tauStar
    ZmgEq=mmgcc*tauDepSfh/tauStar
    ZfeccEq=mfecc*tauDepSfh/tauStar
    ZfeIaEq1=norm1*mfeIa*((tauDepSfh/tauStar)*(tauIaSfh1/tauIa1)*
                           np.exp(tdmin/tauSfh))
    ZfeIaEq2=norm2*mfeIa*((tauDepSfh/tauStar)*(tauIaSfh2/tauIa2)*
                           np.exp(tdmin/tauSfh))

    # deltat rather than t enters for SNIa calculations
    deltat=tmod-tdmin
    # equations 50, 52, and 53 of WAF
    Zo=ZoEq*(1.-np.exp(-tmod/tauDepSfh))
    Zmg=Zo*ZmgEq/ZoEq
    Zfecc=Zo*ZfeccEq/ZoEq
    ZfeIa1=ZfeIaEq1*(1.-np.exp(-deltat/tauDepSfh)-(tauDepIa1/tauDepSfh)* 
                     (np.exp(-deltat/tauIaSfh1)-np.exp(-deltat/tauDepSfh)))
    ZfeIa2=ZfeIaEq2*(1.-np.exp(-deltat/tauDepSfh)-(tauDepIa2/tauDepSfh)* 
                     (np.exp(-deltat/tauIaSfh2)-np.exp(-deltat/tauDepSfh)))

    # zero SNIa contribution before tdmin    
    ZfeIa1[tmod<tdmin]=0.
    ZfeIa2[tmod<tdmin]=0.
    Zfe=Zfecc+ZfeIa1+ZfeIa2 # sum iron contributions

    sfr=np.ones(len(tmod))
    if (tauSfh>0.0):
        sfr=np.exp(-tmod/tauSfh)
    sfrtot=sum(sfr)
    sfr=sfr/sfrtot          # fraction of star formation per bin

    # convert to [O/H], [Mg/H], [Fe/H], [O/Fe], [Mg/Fe]
    OH=np.log10(Zo/SolarO)
    MgH=np.log10(Zmg/SolarMg)
    FeH=np.log10(Zfe/SolarFe)
    OFe=OH-FeH
    MgFe=MgH-FeH
    return t,sfr,OH,MgH,FeH,OFe,MgFe

def waf_linexp(tstart,t0,dtout,mocc,mmgcc,mfecc,mfeIa,r,SolarO,SolarMg,SolarFe,tauIa,tdmin,eta,tauStar,tauSfh):
    """
    Compute [O/H], [Mg/H], [Fe/H], [O/Fe], and [Mg/Fe] tracks using WAF analytic model for 
     linear-exponential (a.k.a. "delayed tau") star formation history

    Input arguments (illustrative values in [  ]):
      tstart = star formation start time [0.5]
      t0 = age of universe at observed redshift [14]
      dtout = output timestep [0.02]
      mocc = IMF-averaged CCSN oxygen yield [0.0071]
      mmgcc = IMF-averaged CCSN magnessium yield [0.0007] 
      mfecc = IMF-averaged CCSN iron yield [0.0005]
      mfeIa = IMF-averaged SNIa iron yield over 10 Gyr [0.0007]
      r = recycling fraction [0.4, based on Kroupa IMF]
      SolarO = solar oxygen mass fraction [0.0073]
      SolarMg = solar magnessium mass fraction [0.0007]
      SolarFe = solar iron mass fraction [0.0014]
      tauIa = e-folding time for SNIa DTD [1.5]
              tauIa <= 0, will use double-exp approximation to power-law
      tdmin = minimum time delay for SNIa [0.05]
              for tauIa <= 0 MUST choose tdmin = 0.05 or 0.15
      eta = outflow mass loading factor [0.3]
      tauStar = SFE efficiency timescale [1.0]
      tauSfh = e-folding timescale of SFH [8.0] 

    Returns:
      t, SFR, OH, MgH, FeH, OFe, and MgFe
      each a 1-d array evaluated at the times in t
      SFR is normalized to sum to 1.0, so for evenly spaced t it is the
      fraction of stars formed in each bin

    Time unit is Gyr
    Reference: Weinberg, Andrews, & Freudenburg 2017, ApJ 837, 183
               Particularly Appendix C
    
    """

    import numpy as np
    
    tmax = t0 - tstart                        # maximum time
    
    tmod = np.arange(dtout, tmax+1.e-5,dtout) # fance() input time 
    
    t = tmod + tstart                         # actual time

    # these near-unit multiplications will be used to ensure that timescales
    # don't actually become equal to each other causing divide-by-zero errors
    xfac=1.0013478
    yfac=1.0016523

    # For tauIa=0, use double-exponential that approximates t^{-1.1} power-law
    # see WAF section 5.1 and Appx C for details
    # In this case, we require tdmin=0.05 or tdmin=0.15
    if (tauIa<0.01):
        if (tdmin==0.15):
            tauIa1=0.5
            tauIa2=5.0
            norm1=0.478
            norm2=0.522
        elif (tdmin==0.05):
            tdmin=0.05
            tauIa1=0.25
            tauIa2=3.5
            norm1=0.493
            norm2=0.507
        else:
            print ("can't set DTD parameters for this tdmin\n")
            sys.exit()
    # otherwise set normalization of second exponential to zero
    else:
        tauIa1=tauIa
        tauIa2=tauIa
        norm1=1.0
        norm2=0.0

    # compute harmonic difference timescales
    tauStar*=xfac
    tauSfh*=yfac
    tauDep=tauStar/(1+eta-r)
    tauDepIa1=1./(1./tauDep-1./tauIa1)
    tauDepIa2=1./(1./tauDep-1./tauIa2)
    if (tauSfh>0.0):
        tauDepSfh=1./(1./tauDep-1./tauSfh)
        tauIaSfh1=1./(1./tauIa1-1./tauSfh)
        tauIaSfh2=1./(1./tauIa2-1./tauSfh)
    else:
        tauDepSfh=tauDep
        tauIaSfh1=tauIa1
        tauIaSfh2=tauIa2
        tauSfh=1.e8

    # compute equilibrium abundances, WAF equations 28-30
    ZoEq=mocc*tauDepSfh/tauStar
    ZmgEq=mmgcc*tauDepSfh/tauStar
    ZfeccEq=mfecc*tauDepSfh/tauStar
    ZfeIaEq1=norm1*mfeIa*((tauDepSfh/tauStar)*(tauIaSfh1/tauIa1)*
                           np.exp(tdmin/tauSfh))
    ZfeIaEq2=norm2*mfeIa*((tauDepSfh/tauStar)*(tauIaSfh2/tauIa2)*
                           np.exp(tdmin/tauSfh))

    # deltat rather than t enters for SNIa calculations
    deltat=tmod-tdmin
    # equations 56-58 of WAF
    Zo=ZoEq*(1.-(tauDepSfh/tmod)*(1.-np.exp(-tmod/tauDepSfh)))
    Zmg=Zo*ZmgEq/ZoEq
    Zfecc=Zo*ZfeccEq/ZoEq
    ZfeIa1=ZfeIaEq1*(tauIaSfh1/t)*(
           deltat/tauIaSfh1 + 
           (tauDepIa1/tauDepSfh)*np.exp(-deltat/tauIaSfh1) + 
           (1.+(tauDepSfh/tauIaSfh1)-(tauDepIa1/tauDepSfh))*
           np.exp(-deltat/tauDepSfh) -
           (1.+(tauDepSfh/tauIaSfh1)) )
    ZfeIa2=ZfeIaEq2*(tauIaSfh2/tmod)*(
           deltat/tauIaSfh2 + 
           (tauDepIa2/tauDepSfh)*np.exp(-deltat/tauIaSfh2) + 
           (1.+(tauDepSfh/tauIaSfh2)-(tauDepIa2/tauDepSfh))*
           np.exp(-deltat/tauDepSfh) -
           (1.+(tauDepSfh/tauIaSfh2)) )

    # zero SNIa contribution before tdmin    
    ZfeIa1[tmod<tdmin]=0.
    ZfeIa2[tmod<tdmin]=0.
    Zfe=Zfecc+ZfeIa1+ZfeIa2  # sum iron contributions

    sfr=np.ones(len(tmod))
    if (tauSfh>0.0):
        sfr=tmod*np.exp(-tmod/tauSfh)
    sfrtot=sum(sfr)
    sfr=sfr/sfrtot           # fraction of star formation per bin

    # convert to [O/H], [Fe/H], [O/Fe]
    OH=np.log10(Zo/SolarO)
    MgH=np.log10(Zmg/SolarMg)
    FeH=np.log10(Zfe/SolarFe)
    OFe=OH-FeH
    MgFe=MgH-FeH
    return t, sfr,OH,MgH,FeH,OFe,MgFe

def fance(tstart,t0,dtout,mocc,mmgcc,mfecc,fret,mfeIa,r,SolarO,SolarMg,SolarFe,tauIa,tdmin,eta,
          tauStar,tau1,tau2):
    """
    Compute [O/H], [Mg/H], [Fe/H], [O/Fe] and [Mg/Fe] tracks using analytic model for rise-fall SFH

    Input arguments (illustrative values in [  ]):
      tstart = star formation start time [0.5]
      t0 = age of universe at observed redshift [14]
      dtout = output timestep [0.02]
      mocc = IMF-averaged CCSN oxygen yield [0.0071]
      mmgcc = IMF-averaged CCSN magnessium yield [0.0007] 
      mfecc = IMF-averaged CCSN iron yield [0.0005]
      mfeIa = IMF-averaged SNIa iron yield over 10 Gyr [0.0007]
      r = recycling fraction [0.4, based on Kroupa IMF]
      SolarO = solar oxygen mass fraction [0.0073]
      SolarMg = solar magnessium mass fraction [0.0007]
      SolarFe = solar iron mass fraction [0.0014]
      tauIa = e-folding time for SNIa DTD [1.5]
              tauIa <= 0, will use double-exp approximation to power-law
      tdmin = minimum time delay for SNIa [0.05]
              for tauIa <= 0 MUST choose tdmin = 0.05 or 0.15
      eta = outflow mass loading factor [0.3]
      tauStar = SFE efficiency timescale [1.0]
      tau1, tau2 = parameters of rise-fall SFH [2.0, 8.0]
                   SFR \propto (1-exp(-t/tau1))*exp(-tau2)
                   tau1=0 -> pure exponential decline
                   tau1=tau2=0 -> constant SFR
    Returns:
      t, SFR, OH, MgH, FeH, OFe, and MgFe
      each a 1-d array evaluated at the times in t
      SFR is normalized to sum to 1.0, so for evenly spaced t it is the
      fraction of stars formed in each bin

    Time unit is Gyr
    
    """

    import numpy as np
   
    tmax = t0 - tstart                        # maximum time
    
    tmod = np.arange(dtout, tmax+1.e-5,dtout) # fance() input time 
    
    t = tmod + tstart                         # actual time

    # fret < 1 is equivalent to reducing CCSN yield for the same eta
    mocc*=fret
    mmgcc*=fret
    mfecc*=fret

    if (tau1<=0):  # rise timescale = 0, so just exp with tauSfh=tau2
        t,sfr,OH,MgH,FeH,OFe,MgFe=waf(tstart,t0,dtout,mocc,mmgcc,mfecc,mfeIa,r,SolarO,SolarMg,SolarFe,tauIa,tdmin,eta,tauStar,tau2)
    
    else:
        
        """
        For tau1>0, the SFH can be written as 
        SFR = constant*[exp(-t/tau2)-exp(-t/tauh)], with tauh = (1/tau1+1/tau2).
        The exponential cases follow WAF, and they can be combined using WAF
        equation (117), even though the second exponential has a negative 
        pre-factor.
        """

        if (tau2<=0):
            tauh=tau1
        else:
            tauh=1/(1/tau1 + 1/tau2)

        t,sfr1,OH1,MgH1,FeH1,OFe1,MgFe1=waf(tstart,t0,dtout,mocc,mmgcc,mfecc,mfeIa,r,SolarO,SolarMg,SolarFe,tauIa,tdmin,
                           eta,tauStar,tau2)
        t,sfr2,OH2,MgH2,FeH2,OFe2,MgFe2=waf(tstart,t0,dtout,mocc,mmgcc,mfecc,mfeIa,r,SolarO,SolarMg,SolarFe,tauIa,tdmin,
                           eta,tauStar,tauh)

        # returned sfr arrays are separately normalized, so we need to recompute
        if (tau2<=0):
            sfr1=np.ones(len(tmod))  # constant SFR
        else:
            sfr1=np.exp(-tmod/tau2)
        sfr2= -1.*np.exp(-tmod/tauh) # note negative pre-factor
        sfr=sfr1+sfr2

    # Now we use WAF eq. 117, but we need to convert to linear metallicity
        Zo1= SolarO*(10**OH1)
        Zmg1=SolarMg*(10**MgH1)
        Zfe1=SolarFe*(10**FeH1)
        Zo2= SolarO*(10**OH2)
        Zmg2=SolarMg*(10**MgH2)
        Zfe2=SolarFe*(10**FeH2)
        Zo=  (sfr1*Zo1+sfr2*Zo2)/sfr
        Zmg= (sfr1*Zmg1+sfr2*Zmg2)/sfr
        Zfe= (sfr1*Zfe1+sfr2*Zfe2)/sfr

    # convert back to log quantities and return
        OH=np.log10(Zo/SolarO)
        MgH=np.log10(Zmg/SolarMg)
        FeH=np.log10(Zfe/SolarFe)
        OFe=OH-FeH
        MgFe=MgH-FeH

    # normalize sfr
        sfrtot=sum(sfr)
        sfr=sfr/sfrtot            # fraction of star formation per bin

    return t,sfr,OH,MgH,FeH,OFe,MgFe

def mdf(xmin,xmax,dx,xh,sfr):
    import numpy as np

    x=np.arange(xmin,xmax,dx)
    y=np.zeros(len(x))
    for i in np.arange(len(xh)):
    # which metallicity bin does this timestep fall into?
        ibin=np.int((xh[i]-xmin)/dx+0.5)
    # values outside range are added to first or last bin
        ibin=np.clip(ibin,0,len(x)-1)
        # add to mdf a quantity proportional to stars formed in that time bin
        y[ibin]+=sfr[i]	
    # normalize to unit integral
    y=y/(np.sum(y)*dx)
    # return bin centers instead of bin lower boundaries
    x+=dx/2.
    return x,y
