from libc.math cimport exp
from sector cimport *

import numpy as np
from numpy import isnan


def get_dust_params_dtype():
    cdef dust_params_t dp[1]
    return np.asarray(<dust_params_t[:1]>dp).dtype


class dust_params:
    @property
    def dtype(self):
        return self._dtype


    @property
    def propNames(self):
        return self._propNames

    
    @property
    def paramNames(self):
        return self._paramNames


    @property
    def fixedParams(self):
        return self._fixedParams
    
    
    def __init__(self, propNames, paramNames, **kwargs):
        self._dtype = get_dust_params_dtype()
        self._propNames = tuple(propNames)
        self._paramNames = tuple(paramNames)
        self._fixedParams = np.full(len(paramNames), np.nan, dtype = np.double)
        self.set_fixed_params(**kwargs)


    def set_fixed_params(self, **kwargs):
        paramNames = self.paramNames
        for name in kwargs.keys():
            if name not in paramNames:
                raise KeyError("%s is not in %s!"%(name, paramNames.__str__()))
        p = np.array(self._fixedParams, copy = True)
        for iN, name in enumerate(paramNames):
            if name in kwargs:
                p[iN] = kwargs[name]
        p.flags.writeable = False
        self._fixedParams = p


    def fixed_names(self):
        return np.array(self.paramNames)[~isnan(self._fixedParams)]


    def full_params(self, params):
        p = np.array(self._fixedParams, copy = True)
        p[isnan(p)] = params
        return p


class SFR_model_params(dust_params):
    def __init__(self, **kwargs):
        super().__init__(('Sfr',), ('tauISM', 'tauBC', 's1', 'n', 'tBC', 'a'), **kwargs)


    def __call__(self, params, props, z):
        tauISM, tauBC, s1, n, tBC, a = self.full_params(params)
        #   -Convert M_solar/yr to 100 M_solar/yr
        SFR = props['Sfr']/100.
        nGal = len(SFR)
        #
        dustParams = np.zeros(nGal, dtype = self._dtype)
        # Total optical depth
        factor = SFR**s1*exp(-a*z)
        # ISM optical depth
        dustParams['tauUV_ISM'] = tauISM*factor
        # ISM reddening slope
        dustParams['nISM'] = np.full(nGal, n, dtype = np.double)
        # Birth cloud optical depth
        dustParams['tauUV_BC'] = tauBC*factor
        # Birth cloud reddening slope
        dustParams['nBC'] = np.full(nGal, n, dtype = np.double)
        # Birth cloud lifetime
        dustParams['tBC'] = np.full(nGal, tBC, dtype = np.double)
        return dustParams


def compute_mags_mhysa(
    double[:] inBCFlux, double[:] outBCFlux,
    double[:] centreWaves, double[:] logWaves, int nBeta, int nFlux,
    dust_params_t[:] dustParams):
    # Add dust
    cdef:
        int iG
        int nGal = dustParams.shape[0]
        double *pIF = &inBCFlux[0]
        double *pOF = &outBCFlux[0]
        double *waves = &centreWaves[0]
        dust_params_t *pDustParams = &dustParams[0]

    for iG in range(nGal):
        dust_absorption_approx(pIF, pOF, waves, nFlux, pDustParams)
        pIF += nFlux
        pOF += nFlux
        pDustParams += 1

    # Fit UV slopes
    beta = np.zeros(nGal, dtype = np.double)
    flux = np.asarray(inBCFlux) + np.asarray(outBCFlux)
    cdef:   
        double[:] mvBeta = beta
        double[:] mvFlux = flux
        double *pBeta = &mvBeta[0]
        double *pFlux = &mvFlux[0]
    waves = &logWaves[0]
    fit_UV_slope(pBeta, pFlux, nGal, nFlux, waves, nBeta, 1)

    # Get M1600
    #   -flux should be in a unit of Jansky
    flux.resize(nGal, nFlux)
    M1600 = -2.5*np.log10(flux[:, nBeta]) + 8.9

    return M1600, beta
