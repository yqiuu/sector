# cython: c_string_encoding = ascii

import os, sys
from warnings import warn
from copy import deepcopy

from libc.stdlib cimport malloc, free
from libc.stdio cimport *
from libc.math cimport exp

import numpy as np, h5py
from numpy import isscalar
from pandas import DataFrame

from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from dragons import meraxes


__all__ = [
    'save_star_formation_history',
    'Lyman_absorption',
    'get_wavelength',
    'HST_filters',
    'beta_filters',
    'composite_spectra',
    'calibration',
    'dust_extinction',
    'reddening_curve',
    'reddening'
]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Functions to compute the IGM absorption                                       #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def Lyman_absorption(obsWaves, z):
    """
    Compute the IGM transmission curve from Inoue et al. 2014

    Parameters
    ----------
    obsWaves: array_like
        Wavelengths in observer frames.
    z: float
        redshift.

    Returns
    -------
    trans: ndarray
        1-D array containing the IGM transmission.
    """
    obsWaves = np.array(obsWaves).flatten()
    trans = np.ones(len(obsWaves))
    cdef:
        double[::1] mvTrans = trans
        double[::1] mvObsWaves = obsWaves
        int nWaves = len(trans)
    add_Lyman_absorption(&mvTrans[0], &mvObsWaves[0], nWaves, <double>z)
    return trans


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Functions about ISM absorptionb                                               #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
cdef dust_params_t *init_dust_parameters(dust):
    cdef:
        int iG
        int nGal = len(dust)
        double[:, ::1] mvDustParams = np.array(dust)
        dust_params_t *dustParams = <dust_params_t*>malloc(nGal*sizeof(dust_params_t))
        dust_params_t *pDustParams = dustParams

    for iG in xrange(nGal):
        pDustParams.tauUV_ISM = mvDustParams[iG, 0]
        pDustParams.nISM = mvDustParams[iG, 1]
        pDustParams.tauUV_BC = mvDustParams[iG, 2]
        pDustParams.nBC = mvDustParams[iG, 3]
        pDustParams.tBC = mvDustParams[iG, 4]
        pDustParams += 1

    return dustParams


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Functions to read SED templates                                               #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def get_wavelength(path):
    #=====================================================================
    # Return wavelengths of SED templates in a unit of angstrom
    #=====================================================================
    return np.array(h5py.File(os.path.join(path, "sed_library.hdf5"), "r").get("waves"))


cdef void free_raw_spectra(sed_params_t *spectra):
    free(spectra.age)
    free(spectra.waves)
    free(spectra.raw)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Functions to process filters                                                  #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
from filters import filterDict

def HST_filters(filterNames):
    """
    Quick access HST filters.

    Parameters
    ----------
    filterNames: list
        Available filters: B435, V606, i775, I814, z850, Y098, Y105,
        J125, H160, 3.6.

    Returns
    -------
    obsBands: list
        For each row, the first element is the filter name, and the
        second element is the transmission curve. The output can be
        passed to ``composite_spectra``.
    """
    obsBands = []
    for name in filterNames:
        obsBands.append([name, np.load(filterDict[name])])
    return obsBands


def beta_filters():
    #=====================================================================
    # return the filters defined by Calzetti et al. 1994, which is used to
    # calculate the UV continuum slope
    #=====================================================================
    windows = np.array([[1268., 1284.],
                        [1309., 1316.],
                        [1342., 1371.],
                        [1407., 1515.],
                        [1562., 1583.],
                        [1677., 1740.],
                        [1760., 1833.],
                        [1866., 1890.],
                        [1930., 1950.],
                        [2400., 2580.]])
    return windows


cdef void generate_filters(sed_params_t *spectra, outType,
                           betaBands, restBands, obsBands, z, obsFrame):
    cdef:
        double *c_betaBands = NULL
        double *c_restBands = NULL
        double *obsTrans = NULL
        double *obsWaves = NULL
        int *nObsWaves = NULL
    # Set redshift
    spectra.z = z
    #
    if outType == "ph":
        nRest = len(restBands)
        if nRest > 0:
            centre, width = np.array(restBands).T
            restBands = np.vstack([centre - width/2., centre + width/2.]).T.flatten()
            c_restBands = init_1d_double(restBands)
        nObs = len(obsBands)
        if nObs > 0:
            trans = np.array([])
            waves = np.array([])
            nObsWaves = <int*>malloc(nObs*sizeof(int))
            for iF in range(nObs):
                trans = np.append(trans, obsBands[iF][1][1])
                waves = np.append(waves, obsBands[iF][1][0])
                nObsWaves[iF] = len(obsBands[iF][1][0])
            obsTrans = init_1d_double(trans)
            obsWaves = init_1d_double(waves)
        init_filters(spectra, NULL, 0, c_restBands, nRest, obsTrans, obsWaves, nObsWaves, nObs, z)
        free(c_restBands)
        free(obsTrans)
        free(obsWaves)
        free(nObsWaves)
    elif outType == "sp":
        spectra.nFlux = spectra.nWaves
        if obsFrame:
            spectra.nObs = 1
        else:
            spectra.nObs = 0
        spectra.nFilterWaves = NULL
        spectra.filterWaves = NULL
        spectra.filters = NULL
        spectra.centreWaves = NULL
        spectra.logWaves = NULL
    elif outType == "UV slope":
        if betaBands == []:
            betaBands = beta_filters()
        c_betaBands = init_1d_double(betaBands.flatten())
        c_restBands = init_1d_double(np.array([1550., 1650.]))
        init_filters(spectra, c_betaBands, len(betaBands), c_restBands, 1, NULL, NULL, NULL, 0, 0.)
        free(c_betaBands)
        free(c_restBands)


cdef void free_filters(sed_params_t *spectra):
    free(spectra.nFilterWaves)
    free(spectra.filterWaves)
    free(spectra.filters)
    free(spectra.centreWaves)
    free(spectra.logWaves)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Primary functions                                                             #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def get_output_name(prefix, postfix, snap, path):
    #=====================================================================
    # Function to generate the name of the output
    #=====================================================================
    fname = prefix + "_%03d"%snap + postfix
    # Avoid repeated name
    idx = 2
    fileList = os.listdir(path)
    while fname in fileList:
        fname = prefix + "_%03d_%d"%(snap, idx) + postfix
        idx += 1
    return os.path.join(path, fname)


def save_star_formation_history(fname, snapList, idxList, h, prefix = 'sfh', outPath = './'):
    """
    Store star formation history to the disk.

    Parameters
    ----------
    fname: str
        Full path to input hdf5 master file.
    snapList: list
        List of snapshots to be computed.
    gals: list
        List of arraies of galaxy indices.
    h: float
        Dimensionless Hubble constant. This is substituded into all
        involved functions in ``meraxes`` python package.
    prefix: str
        The name of the output file is 'prefix_XXX.bin', where XXX is
        number of the snapshot.
    outPath: str
        Path to the output.
    """
    cdef:
        int iS, nSnap
        int snap, snapMax, snapMin
        stellar_population sfh

    if isscalar(snapList):
        snapMax = snapList
        nSnap = 1
        snapList = [snapList]
        idxList = [idxList]
    else:
        snapMax = max(snapList)
        nSnap = len(snapList)
    cdef galaxy_tree_meraxes galData = galaxy_tree_meraxes(fname, snapMax, h)
    # Read and save galaxy merge trees
    for iS in xrange(nSnap):
        outName = get_output_name(prefix, '.bin', snapList[iS], outPath)
        stellar_population(galData, snapList[iS], idxList[iS]).save(outName)


def composite_spectra(fname, snapList, gals, h, Om0, sedPath,
                      dust = None, IGM = 'I2014', approx = False,
                      outType = 'ph',
                      betaBands = [], restBands = [[1600, 100],], obsBands = [], obsFrame = False,
                      timeGrid = 0,
                      prefix = 'mags', outPath = './',
                      nThread = 1):
    """
    Main function to calculate galaxy magnitudes and spectra.

    Parameters
    ----------
    fname: str
        Full path to input hdf5 master file.
    snapList: list
        List of snapshots to be computed.
    gals: list
        Each element of the list can be an array of galaxy indices or
        a path to stored star formation history.
    h: float
        Dimensionless Hubble constant. This is substituded into all
        involved functions in meraxes python package. It is also used
        to calculate the luminosity distance.
    Om0: float
        Current day matter content of the Universe. It is used to
        calculate the luminosity distance.
    sedPath: str
        Full path to SED templates.
    dust: ndarray
        Parameters for the dust model. It should have a shape of
        ``(len(snapList), len(gals), 5)``. The five parameters are
        tauUV_ISM, nISM, tauUV_BC, nBC, tBC.
    IGM: str
        Method to calculate the transmission due to the Lyman
        absorption. It can only be 'I2014'. It is only applicable
        to observer frame quantities.
    outType: str
        If 'ph', output AB magnitudes in filters given by restBands
        and obsBands.

        If 'sp', output full spectra in unit of
        :math:`erg/s/\\unicode{x212B}/cm^2`. if obsFrame is true, flux
        densities is normlised by the luminosity distance;otherwise,
        it is normlised by :math:`10 pc`. Wavelengths are in a unit of
        :math:`\\unicode{x212B}`.

        If 'UV slope', output slopes, normalisations, and correlation
        cofficients by a power law fit at UV range using 10 windows
        given by Calzetti et al. 1994. It also outputs flux densities
        in these windows in a unit of :math:`erg/s/\\unicode{x212B}/cm^2`
        normlised by :math:`10 pc`. Wavelengths are in a unit of
        :math:`\\unicode{x212B}`.
    restBands: list
        List of doublets to specify rest frame filters. The first
        element of the doublet is the centre wavelength, and
        the second one is band width.
    obsBands: list
        List of doublets to specify observer frame filters. The first
        element of the doublet is the filter name, and the second one
        is a 2-D array. The first row of the array is the wavelength
        in a unit of :math:`\\unicode{x212B}`, and the second row gives
        the transmission curve.
    obsFrame: bool
        See ``outType``.
    prefix: str
        The name of the output file is 'prefix_XXX.hdf5', where XXX is
        number of the snapshot.
    outPath: str
        Path to the output.
    nThread: int
        Number of threads used by the OpenMp.

    Returns
    -------
    mags: pandas.DataFrame
        If ``snapList`` is a scalar, it returns the output according to
        ``outType``.

        This function always generates at least one output in the
        directory defined by ``outPath``. The output, whose name is
        defined by ``prefix``, are a ``pandas.DataFrame`` object. Its
        ``index`` is the same with that given in the input. In additon,
        this function never overwrites an output which has the same name;
        instead it generates an output with a different name.
    """
    cosmo = FlatLambdaCDM(H0 = 100.*h, Om0 = Om0)

    cdef:
        int iS, iF, iG
        int nSnap
        int sanpMin = 1
        int snapMax
        galaxy_tree_meraxes galData = None

    if isscalar(snapList):
        snapMax = snapList
        nSnap = 1
        snapList = [snapList]
        gals = [gals]
    else:
        snapMax = max(snapList)
        nSnap = len(snapList)

    # If SFHs are not from files, load outputs from meraxes.
    if type(gals[0]) is not str:
        galData = galaxy_tree_meraxes(fname, snapMax, h)

    waves = get_wavelength(sedPath)
    cdef:
        sed_params_t spectra
        gal_params_t *galParams
        double z
        int nGal

        int nWaves = len(waves)
        int nRest = 0
        int nObs = 0
        int nFlux = 0
        int c_outType = 0

        int nR = 3

        dust_params_t *dustParams = NULL

        double *c_output
        double[:] mvOutput

    for iS in xrange(nSnap):
        # Read star formation rates and metallcities form galaxy merger trees
        sfh = stellar_population(galData, snapList[iS], gals[iS])
        if timeGrid != 0:
            sfh.reconstruct(timeGrid)
        galParams = sfh.pointer()
        z = galParams.z
        nGal = galParams.nGal
        # Convert the format of dust parameters
        if dust is not None:
            dustParams = init_dust_parameters(dust[iS])
        # Compute the transmission of the IGM
        if IGM == 'I2014':
            spectra.igm = 1
        else:
            spectra.igm = 0
        # Read raw SED templates
        init_templates_raw(&spectra, os.path.join(sedPath, "sed_library.hdf5"))
        # Generate Filters
        if outType == 'ph':
            generate_filters(&spectra, outType, [], restBands, obsBands, z, False)
            nRest = len(restBands)
            nObs = len(obsBands)
            nFlux = nRest + nObs
            c_outType = 0
        elif outType == 'sp':
            generate_filters(&spectra, outType, [], [], [], z, obsFrame)
            nFlux = nWaves
            c_outType = 1
        elif outType == 'UV slope':
            generate_filters(&spectra, outType, betaBands, [], [], z, False)
            nFlux = spectra.nFlux
            centreWaves = np.array(<double[:nFlux]>spectra.centreWaves)
            c_outType = 2
        else:
            raise KeyError("outType can only be 'ph', 'sp' and 'UV Slope'")
        shrink_templates_raw(&spectra, galParams.ageStep[galParams.nAgeStep - 1])
        # Compute spectra
        c_output = composite_spectra_cext(&spectra, galParams, dustParams,
                                         c_outType, <short>approx, nThread)
        # Save the output to a numpy array
        if outType == 'UV slope':
            mvOutput = <double[:nGal*(nFlux + nR)]>c_output
            output = np.hstack([np.asarray(mvOutput[nGal*nFlux:],
                                           dtype = 'f4').reshape(nGal, -1),
                                np.asarray(mvOutput[:nGal*nFlux],
                                           dtype = 'f4').reshape(nGal, -1)])
        else:
            mvOutput = <double[:nGal*nFlux]>c_output
            output = np.asarray(mvOutput, dtype = 'f4').reshape(nGal, -1)
        # Convert apparent magnitudes to absolute magnitudes
        if outType == 'ph' and nObs > 0:
            output[:, nRest:] += cosmo.distmod(z).value
        # Convert to observed frame fluxes
        if outType == 'sp' and obsFrame:
            factor = 10./cosmo.luminosity_distance(z).to(u.parsec).value
            output *= factor*factor
        # Set output column names
        if outType == 'ph':
            columns = []
            for iF in xrange(nRest):
                columns.append("M%d-%d"%(restBands[iF][0], restBands[iF][1]))
            for iF in xrange(nObs):
                columns.append(obsBands[iF][0])
        elif outType == 'sp':
            columns = (1. + z)*waves if obsFrame else waves
        elif outType == 'UV slope':
            columns = np.append(["beta", "norm", "R"], centreWaves)
            columns[-1] = "M1600-100"
        # Save the output to the disk
        if type(gals[0]) is str:
            indices = np.asarray(<int[:nGal]>galParams.indices, dtype = 'i4')
        else:
            indices = gals[iS]
        DataFrame(output, index = indices, columns = columns).\
        to_hdf(get_output_name(prefix, ".hdf5", snapList[iS], outPath), "w")

        if len(snapList) == 1:
            mags = DataFrame(deepcopy(output), index = indices, columns = columns)

        #free_gal_params(&galParams)
        free(dustParams)
        free_filters(&spectra)
        free(c_output)

    free_raw_spectra(&spectra)

    if len(snapList) == 1:
        return mags


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Calibration                                                                   #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
cdef class calibration:
    cdef:
        int nSnap
        gal_params_t *galParams
        sed_params_t *spectra
        short approx
        short nThread
        tuple args

    def __cinit__(self, sfhList, sedPath, betaBands = [], approx = False, nThread = 1):
        cdef:
            int iS
            int nSnap = len(sfhList)
            gal_params_t *pGalParams
            sed_params_t *pSpectra
            int nBeta

        self.nSnap = nSnap
        self.galParams = <gal_params_t*>malloc(nSnap*sizeof(gal_params_t))
        self.spectra = <sed_params_t*>malloc(nSnap*sizeof(sed_params_t))
        self.approx = <short>approx
        self.nThread = <short>nThread
        self.args = (sfhList, sedPath, betaBands, approx, nThread)

        pGalParams = self.galParams
        pSpectra = self.spectra
        for iS in xrange(nSnap):
            # Read star formation rates and metallcities form galaxy merger trees
            read_gal_params(pGalParams, sfhList[iS])
            #
            pSpectra.igm = 0
            # Read raw SED templates
            init_templates_raw(pSpectra, os.path.join(sedPath, "sed_library.hdf5"))
            # Generate filters
            generate_filters(pSpectra, "UV slope", betaBands, [], [], pGalParams.z, False)
            #
            shrink_templates_raw(pSpectra, pGalParams.ageStep[pGalParams.nAgeStep - 1])
            #
            pGalParams += 1
            pSpectra += 1


    def __dealloc__(self):
        cdef:
            int iS
            int nSnap = self.nSnap
            gal_params_t *pGalParams = self.galParams
            sed_params_t *pSpectra = self.spectra

        for iS in xrange(nSnap):
            free_gal_params(pGalParams)
            free_filters(pSpectra)
            free_raw_spectra(pSpectra)
            #
            pGalParams += 1
            pSpectra += 1

        free(self.galParams)
        free(self.spectra)


    def __reduce__(self):
        return (rebuild_calibration, self.args)


    def run(self, dust):
        cdef:
            int iS
            int nSnap = self.nSnap
            gal_params_t *pGalParams = self.galParams
            sed_params_t *pSpectra = self.spectra
            int iG, nGal
            dust_params_t *dustParams = NULL
            double *output
            double *pOutput
            int nFlux = pSpectra.nFlux
            int iM1600 = nFlux - 1
            int nR = 3
            int iBeta = 0
            double[:] mvM1600
            double[:] mvBeta
        M1600 = np.empty(nSnap, object)
        beta = np.empty(nSnap, object)
        for iS in xrange(nSnap):
            # Compute spectra
            dustParams = init_dust_parameters(dust[iS])
            output = composite_spectra_cext(pSpectra, pGalParams, dustParams, 2,
                                            self.approx, self.nThread)
            pOutput = output
            nGal = pGalParams.nGal
            mvM1600 = np.zeros(nGal, dtype = 'f8')
            mvBeta = np.zeros(nGal, dtype = 'f8')
            for iG in xrange(nGal):
                mvM1600[iG] = pOutput[iM1600]
                pOutput += nFlux
            for iG in xrange(nGal):
                mvBeta[iG] = pOutput[iBeta]
                pOutput += nR
            M1600[iS] = np.asarray(mvM1600)
            beta[iS] = np.asarray(mvBeta)
            free(dustParams)
            free(output)
            #
            pGalParams += 1
            pSpectra += 1
        return M1600, beta


def rebuild_calibration(sfhList, sedPath, betaBands, approx, nThread):
    return calibration(sfhList, sedPath, betaBands, approx, nThread)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Dust model of Mason et al . 2015                                              #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
from numpy import vectorize
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
def reddening_curve(lam):
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


def reddening(waves, M1600, z, scatter = 0.):
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
    if isscalar(waves):
        return reddening_curve(waves)/reddening_curve(1600.)*A1600
    else:
        waves = np.asarray(waves)
        return reddening_curve(waves)/reddening_curve(1600.)*A1600.reshape(-1, 1)


