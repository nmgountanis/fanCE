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
        self.tmod = np.arange(self.dt, self.t0 - self.t_start + self.dt, self.dt)
        self.t = self.tmod + self.t_start
        self.model_kwargs = dict(
            tauSFE=self.tauSFE,
            tauSFH=self.tauSFH,
            yOCC=self.yOCC,
            yMgCC=self.yMgCC,
            yFeCC=self.yFeCC,
            yFeIa=self.yFeIa,
            r=self.r,
            eta=self.eta,
            tauIa=self.tauIa,
            tDminIa=self.tDminIa,
            SolarO=self.SolarO,
            SolarMg=self.SolarMg,
            SolarFe=self.SolarFe,
            SFH_fn=self.SFH_fn,
            IaDTD_fn=self.IaDTD_fn,
            tmod=self.tmod,
        )

    def update(self, p):
        for p_name, p_val in p.items():
            if p_name[:3] == 'log':
                setattr(self, p_name[3:], 10 ** p_val)
            else:
                setattr(self, p_name, p_val)
        self.tmod = np.arange(self.dt, self.t0 - self.t_start + self.dt, self.dt)
        self.t = self.tmod + self.t_start
        self.model_kwargs = dict(
            tauSFE=self.tauSFE,
            tauSFH=self.tauSFH,
            yOCC=self.yOCC,
            yMgCC=self.yMgCC,
            yFeCC=self.yFeCC,
            yFeIa=self.yFeIa,
            r=self.r,
            eta=self.eta,
            tauIa=self.tauIa,
            tDminIa=self.tDminIa,
            SolarO=self.SolarO,
            SolarMg=self.SolarMg,
            SolarFe=self.SolarFe,
            SFH_fn=self.SFH_fn,
            IaDTD_fn=self.IaDTD_fn,
            tmod=self.tmod,
        )


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
        self.tmod = np.arange(self.dt, self.t0 - self.t_start + self.dt, self.dt)
        self.t = self.tmod + self.t_start
        self.model_kwargs = dict(
            tauSFE=self.tauSFE,
            tauSFH1=self.tauSFH1,
            tauSFH2=self.tauSFH2,
            fRetCC=self.fRetCC,
            yOCC=self.yOCC,
            yMgCC=self.yMgCC,
            yFeCC=self.yFeCC,
            yFeIa=self.yFeIa,
            r=self.r,
            eta=self.eta,
            tauIa=self.tauIa,
            tDminIa=self.tDminIa,
            SolarO=self.SolarO,
            SolarMg=self.SolarMg,
            SolarFe=self.SolarFe,
            IaDTD_fn=self.IaDTD_fn,
            tmod=self.tmod,
        )

    def update(self, p):
        for p_name, p_val in p.items():
            if p_name[:3] == 'log':
                setattr(self, p_name[3:], 10 ** p_val)
            else:
                setattr(self, p_name, p_val)
        self.tmod = np.arange(self.dt, self.t0 - self.t_start + self.dt, self.dt)
        self.t = self.tmod + self.t_start
        self.model_kwargs = dict(
            tauSFE=self.tauSFE,
            tauSFH1=self.tauSFH1,
            tauSFH2=self.tauSFH2,
            fRetCC=self.fRetCC,
            yOCC=self.yOCC,
            yMgCC=self.yMgCC,
            yFeCC=self.yFeCC,
            yFeIa=self.yFeIa,
            r=self.r,
            eta=self.eta,
            tauIa=self.tauIa,
            tDminIa=self.tDminIa,
            SolarO=self.SolarO,
            SolarMg=self.SolarMg,
            SolarFe=self.SolarFe,
            IaDTD_fn=self.IaDTD_fn,
            tmod=self.tmod,
        )
