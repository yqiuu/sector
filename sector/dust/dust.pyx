from libc.math cimport exp
from ..utils cimport *

import numpy as np
from numpy import isnan, isscalar, vectorize


__all__ = [
    'get_dust_params_dtype',
    'dust_params',
    'SFR_model_params',
    'gas_model_params',
    'DTG_model_params',
    'compute_mags_mhysa',
    'reddening_curve_Calzetti00',
    'reddening_Mason15'
]


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


class gas_model_params(dust_params):
    def __init__(self, **kwargs):
        super().__init__(
            ('ColdGas', 'DiskScaleLength'),
            ('tauISM', 'tauBC', 's1', 's2', 'n', 'tBC', 'a'),
            **kwargs
        )


    def __call__(self, params, props, z):
        tauISM, tauBC, s1, s2, n, tBC, a = self.full_params(params)
        #   -Unit 10^10 h^-1 M_solar
        gasMass = props['ColdGas']
        #   -Convert h^-1 Mpc to h^-1 kpc
        radius = props['DiskScaleLength']*1e3
        nGal = len(radius)
        #
        dustParams = np.zeros(nGal, dtype = self._dtype)
        # Total optical depth
        factor = gasMass**s1*radius**s2*exp(-a*z)
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


class DTG_model_params(dust_params):
    def __init__(self, **kwargs):
        super().__init__(
            ('MetalsColdGas', 'ColdGas', 'DiskScaleLength'),
            ('tauISM', 'tauBC', 's1', 's2', 'n', 'tBC', 'a'),
            **kwargs
        )


    def __call__(self, params, props, z):
        tauISM, tauBC, s1, s2, n, tBC, a = self.full_params(params)
        #   -Unit 10^10 h^-1 M_solar
        metalsMass = props['MetalsColdGas']
        gasMass = props['ColdGas']
        metallicity = self.calc_metallicity(metalsMass, gasMass)
        #   -Convert h^-1 Mpc to h^-1 kpc
        radius = props['DiskScaleLength']*1e3
        nGal = len(radius)
        #
        dustParams = np.zeros(nGal, dtype = self._dtype)
        # Total optical depth
        factor = metallicity**s1*gasMass*radius**s2*exp(-a*z)
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


    def calc_metallicity(self, metalsMass, gasMass):
        cond = gasMass > 0.
        metallicity = np.zeros(len(gasMass))
        #   -Convert to Z_solar
        metallicity[cond] = metalsMass[cond]/gasMass[cond]/0.02
        return metallicity


def compute_mags_mhysa(
    double[:] inBCFlux, double[:] outBCFlux,
    double[:] centreWaves, double[:] logWaves, int nBeta, int nFlux,
    dust_params_t[:] dustParams
):
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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Dust model of Mason et al . 2015                                              #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
from warnings import warn

from scipy.interpolate import interp1d
from scipy.optimize import brentq


DEF DUST_C = -2.33
DEF DUST_M0 = -19.5
DEF DUST_SIGMA = .34
DEF DUST_BRIGHTER = -35.
DEF DUST_FAINTER = 0.
DEF DUST_BOUND = -5.


cdef double beta_MUV(double obsMag, double slope, double inter):
    if obsMag >= DUST_M0:
        return (inter - DUST_C)*exp(slope*(obsMag - DUST_M0)/(inter - DUST_C)) \
               + DUST_C
    else:
        return slope*(obsMag - DUST_M0) + inter


cdef dust_equation(double obsMag, double slope, double inter,
                   double insMag, double noise):
    return obsMag - insMag \
           - (4.43 + 1.99*(beta_MUV(obsMag, slope, inter) + noise))


def dust_extinction(M1600, double z, double scatter):
    #=====================================================================
    # Calculate the dust extinction at rest frame 1600 angstrom
    #
    # M1600: rest frame 1600 angstrom magnitudes. It can be an array.
    # z: redshift
    #
    # Returns: dust extinction at rest frame 1600 angstrom
    #          M1600_obs = M1600 + A1600,
    #          where M1600_obs is the dust attenuated magnitudes
    # Reference Mason et al. 2015, equation 4
    #           Bouwens 2014 et al. 2014, Table 3
    #=====================================================================

    cdef:
        int iM
        int nM
    if isscalar(M1600):
        nM = 1
        M1600 = np.array([M1600], dtpye = 'f8')
    else:
        nM = len(M1600)
        M1600 = np.asarray(M1600, dtype = 'f8')
    cdef:
        double insMag
        double[:] mvM1600 = M1600
        double[:] mvA1600 = np.zeros(nM)
        double[:] mvScatter
        double slope = interp1d([2.5, 3.8, 5., 5.9, 7., 8.],
                                [-.2, -.11, -.14, -.2, -.2, -.15],
                                fill_value = 'extrapolate')(z)
        double inter = interp1d([2.5, 3.8, 5., 5.9, 7., 8.],
                                [-1.7, -1.85, -1.91, -2., -2.05, -2.13],
                                fill_value = 'extrapolate')(z)

    if scatter != 0.:
        mvScatter = np.random.normal(0., scatter, nM)
        for iM in xrange(nM):
            insMag = mvM1600[iM]
            if insMag < DUST_BOUND:
                mvA1600[iM] = brentq(dust_equation, DUST_BRIGHTER, DUST_FAINTER,
                                     args = (slope, inter, mvM1600[iM], mvScatter[iM])) \
                              - mvM1600[iM]
            else:
                mvA1600[iM] = 0.
    else:
        for iM in xrange(nM):
            insMag = mvM1600[iM]
            if insMag < DUST_BOUND:
                mvA1600[iM] = brentq(dust_equation, DUST_BRIGHTER, DUST_FAINTER,
                                     args = (slope, inter, mvM1600[iM], scatter)) \
                              - mvM1600[iM]
            else:
                mvA1600[iM] = 0.
    A1600 = np.asarray(mvA1600)
    A1600[A1600 < 0.] = 0.
    return A1600


@vectorize
def reddening_curve_Calzetti00(lam):
    #=====================================================================
    # Function of the reddening curve of Calzetti et al. 2000
    #
    # lam: wavelengths in a unit of angstrom
    # Reference Calzetti et al. 2000, Liu et al. 2016
    #=====================================================================
    lam *= 1e-4 # Convert angstrom to mircometer
    if lam < .12 or lam > 2.2:
        warn("Warning: wavelength is beyond the range of the reddening curve")
    if lam < .12:
        return -92.44949*lam + 23.21331
    elif lam < .63:
        return 2.659*(-2.156 + 1.509/lam - 0.198/lam**2 + 0.011/lam**3) + 4.05
    elif lam < 2.2:
        return  2.659*(-1.857 + 1.040/lam) + 4.05
    else:
        return max(0., -.57136*lam + 1.62620)


def reddening_Mason15(waves, M1600, z, scatter = 0.):
    """
    Compute the dust extinction at given wavelengths.

    Parameters
    ----------
    waves: array_like
        Wavelength in a unit of :math:`\\unicode{x212B}`.
    M1600: array_like
        Magnitudes at rest-frame 1600 :math:`\\unicode{x212B}`.
    z: float
        redshift.
    scatter: float
        Add a Gaussian scatter to the Meurer relation. If 0, no
        scatter is applied.

    Returns
    -------
    A: array_like
        Dust extinction at given wavelengths, which is additive to AB
        magnitudes. It has a dimension of ``(len(M1600), len(waves))``.
    """
    A1600 = dust_extinction(M1600, z, scatter)
    reddening = reddening_curve_Calzetti00(waves)/reddening_curve_Calzetti00(1600.)
    #
    waves = np.ravel(waves)
    if len(waves) == 1:
        return reddening*A1600
    else:
        return reddening*A1600.reshape(-1, 1)
