import numpy as np

SolarO = 0.00733
SolarMg = 0.000671
SolarFe = 0.00137

aFeCC = 0.45
aFeEq = 0.05
yOCC = 0.973 * SolarO * (0.00137 / SolarFe) * (10**(aFeCC-0.45))
yFeCC = yOCC * (SolarFe / SolarO) * (10**(-aFeCC))
yMgCC = yOCC * SolarMg / SolarO
yFeIa = yFeCC * (10.**(aFeCC - aFeEq) - 1.)


class DefaultWAFParSet:
    def __init__(
        self,
    ):
        self.tauSFE = 1.0
        self.tauSFH = 8.0
        self.yOCC = yOCC
        self.yMgCC = yMgCC
        self.yFeCC = yFeCC
        self.yFeIa = yFeIa
        self.r = 0.4
        self.eta = 0.3
        self.tauIa = 1.5
        self.tDminIa = 0.05
        self.SolarO = SolarO
        self.SolarMg = SolarMg
        self.SolarFe = SolarFe
        self.SFH_fn = 'exponential'
        self.IaDTD_fn = 'powerlaw'
        self.t_start = 0.5
        self.t0 = 14.0
        self.dt = 0.02


class DefaultfanCEParSet:
    def __init__(
        self,
    ):
        self.tauSFE = 1.0
        self.tauSFH1 = 2.0
        self.tauSFH2 = 8.0
        self.fRetCC = 1.0
        self.yOCC = yOCC
        self.yMgCC = yMgCC
        self.yFeCC = yFeCC
        self.yFeIa = yFeIa
        self.r = 0.4
        self.eta = 0.3
        self.tauIa = 1.5
        self.tDminIa = 0.15
        self.SolarO = SolarO
        self.SolarMg = SolarMg
        self.SolarFe = SolarFe
        self.IaDTD_fn = 'powerlaw'
        self.t_start = 0.5
        self.t0 = 14.0
        self.dt = 0.02
