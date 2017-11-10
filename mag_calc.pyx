import os
from warnings import warn
from time import time
from struct import pack, unpack

from cython import boundscheck, wraparound
from libc.stdlib cimport malloc, free
from libc.math cimport exp, log

import numpy as np
from numpy import isnan, isscalar, vectorize
from pandas import DataFrame

from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from dragons import meraxes

# Global variables
#========================================================================================
global sTime, inputDict
inputDict = "/lustre/projects/p113_astro/yqiu/magcalc/input/"
filterList = {"B435":os.path.join(inputDict, "HST_ACS_F435W.npy"), 
              "V606":os.path.join(inputDict, "HST_ACS_F606W.npy"), 
              "i775":os.path.join(inputDict, "HST_ACS_F775W.npy"), 
              "I814":os.path.join(inputDict, "HST_ACS_F814W.npy"), 
              "z850":os.path.join(inputDict, "HST_ACS_F850LP.npy"), 
              "Y098":os.path.join(inputDict, "HST_IR_F098M.npy"), 
              "Y105":os.path.join(inputDict, "HST_IR_F105W.npy"), 
              "J125":os.path.join(inputDict, "HST_IR_F125W.npy"), 
              "H160":os.path.join(inputDict, "HST_IR_F160W.npy"), 
              "3.6":os.path.join(inputDict,  "HST_IRAC_3.6.npy")}

cdef:
    int **g_firstProgenitor
    int **g_nextProgenitor
    float **g_metals
    float **g_sfr
#========================================================================================


cdef int *init_1d_int(int[:] memview):
    cdef:
        int nSize = memview.shape[0]
        int *p = <int*>malloc(nSize*sizeof(int))
        int[:] cMemview = <int[:nSize]>p
    cMemview[...] = memview
    return p


cdef float *init_1d_float(float[:] memview):
    cdef:
        int nSize = memview.shape[0]
        float *p = <float*>malloc(nSize*sizeof(float))
        float[:] cMemview = <float[:nSize]>p
    cMemview[...] = memview
    return p


cdef double *init_1d_double(double[:] memview):
    cdef:
        int nSize = memview.shape[0]
        double *p = <double*>malloc(nSize*sizeof(double))
        double[:] cMemview = <double[:nSize]>p
    cMemview[...] = memview
    return p


def timing_start(text):
    global sTime
    sTime = time()
    print "#**********************************************************"
    print text
 

def timing_end():
    global sTime
    elapsedTime = int(time() - sTime)
    print "# Done!"
    print "# Elapsed time: %i min %i sec"%(elapsedTime/60, elapsedTime%60)
    print "#**********************************************************\n"


def get_wavelength():
    """
    Return wavelengths of SED templates in a unit of angstrom
    """
    global inputDict
    fp = open(os.path.join(inputDict, "sed_waves.bin"), "rb")
    nWaves = unpack('i', fp.read(sizeof(int)))[0]
    waves = np.array(unpack('%dd'%nWaves, fp.read(nWaves*sizeof(double))))
    fp.close()
    return waves


def HST_filters(filterNames):
    """
    Quick access of transmission curves of HST filters

    filterNames can be a list of filter names. (B435, V606, i775,
    I814, z850, Y098, Y105, J125, H160, 3.6)
    
    Output is a 2-d list. Rows are different filters. First column
    is the filter names. Second column is the transmission curve.

    The output is to be passed to galaxy_mags(...)
    """
    global filterList
    obsBands = []
    for name in filterNames:
        obsBands.append([name, np.load(filterList[name])])
    return obsBands


def read_filters(restFrame, obsBands, z):
    """
    This function is to generate transmission curves that has the 
    same wavelengths with SED templates. It is called by galaxy_mags(...). 
    The input format refer to galaxy_mags(...). 

    Before integration over the filters, the fluxes must be a function of wavelength.
    After integration over the filters, the fluxex becomes a function of frequency.
    """
    waves = get_wavelength()
    nRest = len(restFrame)
    nObs = len(obsBands)
    filters = np.zeros([nRest + nObs, len(waves)])
    obsWaves = (1 + z)*waves
    for i in xrange(nRest):
        centre, bandWidth = restFrame[i]
        lower = centre - bandWidth/2.
        upper = centre + bandWidth/2.
        filters[i] = np.interp(waves, [lower, upper], [1., 1.], left = 0., right = 0.)
        filters[i] /= np.trapz(filters[i]/waves, waves)
        filters[i] *= 3.34e4*waves
    for i in xrange(nObs):
        fWaves, trans = obsBands[i][1]
        filters[nRest + i] = np.interp(obsWaves, fWaves, trans, left = 0., right = 0.)
        filters[nRest + i] /= np.trapz(filters[nRest + i]/waves, waves)
        filters[nRest + i] *= 3.34e4*obsWaves
    return filters.flatten()


def beta_filters():
    """
    return the filters defined by Calzetti et al. 1994, which is used to calculate
    the UV continuum slope
    """
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
    waves = get_wavelength()
    nFilter = len(windows)
    filters = np.zeros([nFilter + 1, len(waves)])
    for iF in xrange(nFilter):
        filters[iF] = np.interp(waves, windows[iF], [1., 1.], left = 0., right = 0.)
        filters[iF] /= np.trapz(filters[iF], waves)
    filters[-1] = read_filters([[1600., 100.]], [], 0.)
    centreWaves = np.append(windows.mean(axis = 1), 1600.)
    return centreWaves, filters.flatten()
        

def get_output_name(prefix, postfix, snap, path):
    """
    Function to generate the name of the output
    """
    fname = prefix + "_%03d"%snap + postfix
    # Avoid repeated name
    idx = 2
    fileList = os.listdir(path)
    while fname in fileList:
        fname = prefix + "_%03d_%d"%(snap, idx) + postfix
        idx += 1
    return os.path.join(path, fname)


def read_meraxes(fname, int snapMax, h):
    """
    This function reads meraxes output. It is called by galaxy_mags(...).
    Meraxes output is stored by g_firstProgenitor, g_nextProgenitor, g_metals
    and g_sfr. They are external variables of mag_calc_cext.c

    fname: path of the meraxes output
    snapMax: start snapshot
    h: liitle h

    Return: the smallest snapshot number that contains a galaxy
    """
    timing_start("# Read meraxes output")
    cdef:
        int snapNum = snapMax+ 1
        int snapMin = snapMax
        int snap, N
        int[:] intMemview1, intMemview2
        float[:] floatMemview1, floatMemview2
    global g_firstProgenitor 
    global g_nextProgenitor
    global g_metals
    global g_sfr
    g_firstProgenitor = <int**>malloc(snapNum*sizeof(int*))
    g_nextProgenitor = <int**>malloc(snapMax*sizeof(int*))
    g_metals = <float**>malloc(snapNum*sizeof(float*))
    g_sfr = <float**>malloc(snapNum*sizeof(float*))

    meraxes.set_little_h(h = h)
    for snap in xrange(snapMax, -1, -1):
        try:
            # Copy metallicity and star formation rate to the pointers
            gals = meraxes.io.read_gals(fname, snap, 
                                        ["StellarMass", "MetalsStellarMass", "Sfr"])
            metals = gals["MetalsStellarMass"]/gals["StellarMass"]
            metals[isnan(metals)] = 0.001
            g_metals[snap] = init_1d_float(metals)
            g_sfr[snap] = init_1d_float(gals["Sfr"])
            snapMin = snap
            gals = None
        except IndexError:
            print "# There is no galaxies in snapshot %d"%snap
            break;
    print "# snapMin = %d"%snapMin
    for snap in xrange(snapMin, snapNum):
        # Copy first progenitor indices to the pointer
        g_firstProgenitor[snap] = \
        init_1d_int(meraxes.io.read_firstprogenitor_indices(fname, snap))
        # Copy next progenitor indices to the pointer
        if snap < snapMax:
            g_nextProgenitor[snap] = \
            init_1d_int(meraxes.io.read_nextprogenitor_indices(fname, snap))

    timing_end()    
    return snapMin


cdef void free_meraxes(int snapMin, int snapMax):
    """
    Function to free g_firstProgenitor, g_nextProgenitor, g_metals and g_sfr
    """
    cdef int i
    # There is no indices in g_nextProgenitor[snapMax]
    for i in xrange(snapMin, snapMax):
        free(g_nextProgenitor[i])

    snapMax += 1
    for i in xrange(snapMin, snapMax):
        free(g_firstProgenitor[i])
        free(g_metals[i])
        free(g_sfr[i])

    free(g_firstProgenitor)
    free(g_nextProgenitor)
    free(g_metals)
    free(g_sfr)

cdef extern from "mag_calc_cext.h":
    struct props:
        short index
        short metals
        float sfr

    struct prop_set:
        props *nodes
        int nNode

    prop_set *read_properties_by_progenitors(int **firstProgenitor, int **nextProgenitor,
                                             float **galMetals, float **galSFR,
                                             int tSnap, int *indices, int nGal)


def trace_star_formation_history(fname, snap, galIndices, h):
    # Read galaxy properties from Meraxes outputs
    cdef int snapMin = read_meraxes(fname, snap, h)
    # Trace galaxy merge trees
    cdef:
        int iG
        int nGal = len(galIndices)
        int *indices = init_1d_int(np.asarray(galIndices, dtype = 'i4'))
        prop_set *galProps = \
        read_properties_by_progenitors(g_firstProgenitor, g_nextProgenitor, g_metals, g_sfr,
                                       snap, indices, nGal)
    free(indices)
    free_meraxes(snapMin, snap)
    # Convert output to numpy array
    cdef:
        int iN
        int nNode
        props *nodes
        double[:, ::1] mvNodes
    output = np.empty(nGal, dtype = object)
    for iG in xrange(nGal):
        nNode = galProps[iG].nNode
        nodes = galProps[iG].nodes
        mvNodes = np.zeros([nNode, 3])
        for iN in xrange(nNode):
            mvNodes[iN][0] = nodes[iN].index
            mvNodes[iN][1] = nodes[iN].metals
            mvNodes[iN][2] = nodes[iN].sfr
        output[iG] = np.asarray(mvNodes)
    return output


def save_star_formation_history(fname, snapList, idxList, h, 
                                prefix = 'sfh', path = './'):
    # Read galaxy properties form Meraxes outputs
    cdef:
        int iS, nSnap
        int snap, snapMax, snapMin
    if isscalar(snapList):
        snapMax = snapList
        nSnap = 1
        snapList = [snapList]
        idxList = [idxList]
    else:
        snapMax = max(snapList)
        nSnap = len(snapList)
    snapMin = read_meraxes(fname, snapMax, h)
    # Read and save galaxy merge trees
    cdef:
        int iG, nGal
        int *indices
        prop_set *galProps

        int iN, nNode
        props *pNodes
    for iS in xrange(nSnap):
        snap = snapList[iS]
        fp = open(get_output_name(prefix, ".bin", snap, path), "wb")
        galIndices = idxList[iS]
        nGal = len(galIndices)
        fp.write(pack('i', nGal))
        fp.write(pack('%di'%nGal, *galIndices))
        indices = init_1d_int(np.asarray(galIndices, dtype = 'i4'))
        galProps = read_properties_by_progenitors(g_firstProgenitor, g_nextProgenitor, 
                                                  g_metals, g_sfr, snap, indices, nGal)
        free(indices)
        for iG in xrange(nGal):
            nNode = galProps[iG].nNode
            fp.write(pack('i', nNode))
            pNodes = galProps[iG].nodes
            for iN in xrange(nNode):
                fp.write(pack('hhf', pNodes.index, pNodes.metals, pNodes.sfr))
                pNodes += 1
        fp.close()
    free_meraxes(snapMin, snapMax)


cdef prop_set *read_properties_by_file(name):
    fp = open(name, "rb")
    cdef:
        int iG
        int nGal = unpack('i', fp.read(sizeof(int)))[0]
        prop_set *galProps = <prop_set*>malloc(nGal*sizeof(prop_set))
        prop_set *pGalProps = galProps

        int iN, nNode
        props *pNodes
    fp.read(nGal*sizeof(int)) # Skip galaxy indices
    for iG in xrange(nGal):
        pGalProps = galProps + iG
        nNode = unpack('i', fp.read(sizeof(int)))[0]
        pGalProps.nNode = nNode
        pNodes = <props*>malloc(nNode*sizeof(props))
        pGalProps.nodes = pNodes
        for iN in xrange(nNode):
            pNodes.index = unpack('h', fp.read(sizeof(short)))[0]
            pNodes.metals = unpack('h', fp.read(sizeof(short)))[0]
            pNodes.sfr = unpack('f', fp.read(sizeof(float)))[0]
            pNodes += 1
    fp.close()
    return galProps


def read_galaxy_indices(name):
    fp = open(name, "rb")
    nGal = unpack('i', fp.read(sizeof(int)))[0]
    indices = np.array(unpack('%di'%nGal, fp.read(nGal*sizeof(int))))
    fp.close()
    return indices

def Lyman_absorption_Fan(double[:] obsWaves, double z):
    """
    Depreciate function. It is original to calculate the optical depth of
    Fan et al. 2006 
    """
    cdef:
        int i
        int nWaves = obsWaves.shape[0]
        double[:] absorption = np.zeros(nWaves)
        double tau
        double ratio
    for i in xrange(nWaves):
        ratio = obsWaves[i]/1216.
        if ratio < 1. + z:
            if ratio < 6.5:
                tau = .85*(ratio/5.)**4.3
            else:
                tau = .15*(ratio/5.)**10.9
        else:
            tau = 0.
        absorption[i] = exp(-tau)
    
    return np.asarray(absorption)


DEF NLYMAN = 39 # Inoue calculated the absorption of 40th Lyman series

def Lyman_absorption_Inoue(double[:] obsWaves, double z):
    """
    Function to calculate the optical depth of Inoue et al. 2014
    It is called by galaxy_mags(...).

    obsWaves: wavelength in unit of angstrom
    z: redshift

    Return: transmission (dimensionless)
    """
    # Reference Inoue et al. 2014
    cdef:
        double LymanSeries[NLYMAN]
        double LAF1[NLYMAN]
        double LAF2[NLYMAN]
        double LAF3[NLYMAN]
        double DLA1[NLYMAN]
        double DLA2[NLYMAN]

    LymanSeries[:] = [1215.67, 1025.72, 972.537, 949.743, 937.803,
                      930.748, 926.226, 923.150, 920.963, 919.352,
                      918.129, 917.181, 916.429, 915.824, 915.329,
                      914.919, 914.576, 914.286, 914.039, 913.826,
                      913.641, 913.480, 913.339, 913.215, 913.104,
                      913.006, 912.918, 912.839, 912.768, 912.703,
                      912.645, 912.592, 912.543, 912.499, 912.458,
                      912.420, 912.385, 912.353, 912.324]
    LAF1[:] = [1.690e-02, 4.692e-03, 2.239e-03, 1.319e-03, 8.707e-04,
               6.178e-04, 4.609e-04, 3.569e-04, 2.843e-04, 2.318e-04,
               1.923e-04, 1.622e-04, 1.385e-04, 1.196e-04, 1.043e-04,
               9.174e-05, 8.128e-05, 7.251e-05, 6.505e-05, 5.868e-05,
               5.319e-05, 4.843e-05, 4.427e-05, 4.063e-05, 3.738e-05,
               3.454e-05, 3.199e-05, 2.971e-05, 2.766e-05, 2.582e-05,
               2.415e-05, 2.263e-05, 2.126e-05, 2.000e-05, 1.885e-05,
               1.779e-05, 1.682e-05, 1.593e-05, 1.510e-05]
    LAF2[:] = [2.354e-03, 6.536e-04, 3.119e-04, 1.837e-04, 1.213e-04,
               8.606e-05, 6.421e-05, 4.971e-05, 3.960e-05, 3.229e-05,
               2.679e-05, 2.259e-05, 1.929e-05, 1.666e-05, 1.453e-05,
               1.278e-05, 1.132e-05, 1.010e-05, 9.062e-06, 8.174e-06,
               7.409e-06, 6.746e-06, 6.167e-06, 5.660e-06, 5.207e-06,
               4.811e-06, 4.456e-06, 4.139e-06, 3.853e-06, 3.596e-06,
               3.364e-06, 3.153e-06, 2.961e-06, 2.785e-06, 2.625e-06,
               2.479e-06, 2.343e-06, 2.219e-06, 2.103e-06]
    LAF3[:] = [1.026e-04, 2.849e-05, 1.360e-05, 8.010e-06, 5.287e-06,
               3.752e-06, 2.799e-06, 2.167e-06, 1.726e-06, 1.407e-06,
               1.168e-06, 9.847e-07, 8.410e-07, 7.263e-07, 6.334e-07,
               5.571e-07, 4.936e-07, 4.403e-07, 3.950e-07, 3.563e-07,
               3.230e-07, 2.941e-07, 2.689e-07, 2.467e-07, 2.270e-07,
               2.097e-07, 1.943e-07, 1.804e-07, 1.680e-07, 1.568e-07,
               1.466e-07, 1.375e-07, 1.291e-07, 1.214e-07, 1.145e-07,
               1.080e-07, 1.022e-07, 9.673e-08, 9.169e-08]
    DLA1[:] = [1.617e-04, 1.545e-04, 1.498e-04, 1.460e-04, 1.429e-04,
               1.402e-04, 1.377e-04, 1.355e-04, 1.335e-04, 1.316e-04,
               1.298e-04, 1.281e-04, 1.265e-04, 1.250e-04, 1.236e-04,
               1.222e-04, 1.209e-04, 1.197e-04, 1.185e-04, 1.173e-04,
               1.162e-04, 1.151e-04, 1.140e-04, 1.130e-04, 1.120e-04,
               1.110e-04, 1.101e-04, 1.091e-04, 1.082e-04, 1.073e-04,
               1.065e-04, 1.056e-04, 1.048e-04, 1.040e-04, 1.032e-04,
               1.024e-04, 1.017e-04, 1.009e-04, 1.002e-04]
    DLA2[:] = [5.390e-05, 5.151e-05, 4.992e-05, 4.868e-05, 4.763e-05, 
               4.672e-05, 4.590e-05, 4.516e-05, 4.448e-05, 4.385e-05, 
               4.326e-05, 4.271e-05, 4.218e-05, 4.168e-05, 4.120e-05,
               4.075e-05, 4.031e-05, 3.989e-05, 3.949e-05, 3.910e-05, 
               3.872e-05, 3.836e-05, 3.800e-05, 3.766e-05, 3.732e-05,
               3.700e-05, 3.668e-05, 3.637e-05, 3.607e-05, 3.578e-05,
               3.549e-05, 3.521e-05, 3.493e-05, 3.466e-05, 3.440e-05,
               3.414e-05, 3.389e-05, 3.364e-05, 3.339e-05]

    cdef:
        int i, j
        int nWaves = obsWaves.shape[0]
        double[:] absorption = np.zeros(nWaves)
        double tau
        double lamObs, ratio

    for i in xrange(nWaves):
        tau = 0.
        lamObs = obsWaves[i]
        # Lyman series
        for j in xrange(NLYMAN):
            ratio = lamObs/LymanSeries[j]
            if ratio < 1. + z:
                # LAF terms
                if ratio < 2.2:
                    tau += LAF1[j]*ratio**1.2
                elif ratio < 5.7:
                    tau += LAF2[j]*ratio**3.7
                else:
                    tau += LAF3[j]*ratio**5.5
                # DLA terms
                if ratio < 3.:
                    tau += DLA1[j]*ratio**2.
                else:
                    tau += DLA2[j]*ratio**3.
        # Lyman continuum
        ratio = lamObs/912.
        # LAF terms
        if z < 1.2:
            if ratio < 1. + z:
                tau += .325*(ratio**1.2 - (1. + z)**-.9*ratio**2.1)
        elif z < 4.7:
            if ratio < 2.2:
                tau += 2.55e-2*(1. + z)**1.6*ratio**2.1 + .325*ratio**1.2 - .25*ratio**2.1
            elif ratio < 1. + z:
                tau += 2.55e-2*((1. + z)**1.6*ratio**2.1 - ratio**3.7)
        else:
            if ratio < 2.2:
                tau += 5.22e-4*(1. + z)**3.4*ratio**2.1 + .325*ratio**1.2 - 3.14e-2*ratio**2.1
            elif ratio < 5.7:
                tau += 5.22e-4*(1. + z)**3.4*ratio**2.1 + .218*ratio**2.1 - 2.55e-2*ratio**3.7
            elif ratio < 1. + z:
                tau += 5.22e-4*((1. + z)**3.4*ratio**2.1 - ratio**5.5)
        # DLA terms
        if z < 2.:
            if ratio < 1. + z:
                tau += .211*(1. + z)**2. - 7.66e-2*(1. + z)**2.3*ratio**-.3 - .135*ratio**2.
        else:
            if ratio < 3.:
                tau += .634 + 4.7e-2*(1. + z)**3. - 1.78e-2*(1. + z)**3.3*ratio**-.3 \
                       -.135*ratio**2. - .291*ratio**-.3
            elif ratio < 1. + z:
                tau += 4.7e-2*(1. + z)**3. - 1.78e-2*(1. + z)**3.3*ratio**-.3 \
                       -2.92e-2*ratio**3.
        absorption[i] = exp(-tau)

    return np.asarray(absorption)


def get_age_list(fname, snap, nAgeList, h):
    """
    Function to generate an array of stellar ages. It is called by galaxy_mags(...).
    """
    travelTime = meraxes.io.read_snaplist(fname, h)[2]*1e6 # Convert Myr to yr
    ageList = np.zeros(nAgeList)
    for i in xrange(nAgeList):
        ageList[i] = travelTime[snap - i - 1] - travelTime[snap]
    return ageList


cdef extern from "mag_calc_cext.h":
    void free_raw_spectra()

    void free_int_spectra()

    struct dust_params:
        double tauUV_ISM
        double nISM
        double tauUV_BC
        double nBC
        double tBC

    float *composite_spectra_cext(prop_set *galProps, int nGal,
                                  double z, double *ageList, int nAgeList,
                                  double *filters, int nRest, int nObs, int mAB,
                                  double *absorption, dust_params *dustArgs)

    float *UV_slope_cext(prop_set *galProps, int nGal,
                         double z, double *ageList, int nAgeList,
                         double *logWaves, double *filters, int nFilter,
                         dust_params *dustArgs)


cdef dust_params *dust_parameters(dustParams):
    cdef:
        int iG
        int nGal = len(dustParams)
        double[:, ::1] mvDustParams = np.array(dustParams)
        dust_params *dustArgs = <dust_params*>malloc(nGal*sizeof(dust_params))
        dust_params *pDustArgs 

    for iG in xrange(nGal):
        pDustArgs = dustArgs + iG
        pDustArgs.tauUV_ISM = mvDustParams[iG, 0]
        pDustArgs.nISM = mvDustParams[iG, 1]
        pDustArgs.tauUV_BC = mvDustParams[iG, 2]
        pDustArgs.nBC = mvDustParams[iG, 3]
        pDustArgs.tBC = mvDustParams[iG, 4]

    return dustArgs
    


def composite_spectra(fname, snapList, idxList, h, Om0, 
                      sfhPath = None,
                      restFrame = [], obsBands = [],
                      obsFrame = False,
                      IGM = 'I2014',
                      dustParams = None,
                      prefix = "mags", path = "./"):
    """
    Main function to calculate galaxy magnitudes
    
    fname: path of meraxes output

    snapList & idxList example:

    snapList = [100, 78]
    idxList = [[0, 1, 2], [100, 101]]
    
    The above means that the function will compute the magnitudes of galaxy 0, 1, 2
    at snapshot 100, and galaxy 100, 101 at snapshot 78.

    h: little h
    Om0: matter content of the universe (necessary to calculate luminosity distance)

    restFrame example:

    restFrame = [[1500., 100], [1600., 50.], [1700., 150.]]
    
    The above means that the function will compute three kinds rest frame magnitudes
    They are centred at 1500 angstrom with filter width 100 anstrom, 
    centred at 1600 angstrom with filter width 50. angstrom, 
    and centred at 1700 angstrom with filter width 150. angstrom

    obsBands: observed frmae magnitudes to be calculated. It can be the output of 
    HST_filters(...).

    Return: the function will store the output as a pandas hdf file. All results are 
    in the AB magnitude.
    """
    cosmo = FlatLambdaCDM(H0 = 100.*h, Om0 = Om0)
   
    cdef:
        int i, iG
        int snap, nSnap
        int sanpMin, snapMax

    if isscalar(snapList):
        snapMax = snapList
        nSnap = 1
        snapList = [snapList]
        idxList = [idxList]
        if sfhPath is not None:
            sfhPath = [sfhPath]
    else:
        snapMax = max(snapList)
        nSnap = len(snapList)

    if sfhPath is None:
        snapMin = read_meraxes(fname, snapMax, h)
    else:
        snapMin = 1

    waves = get_wavelength()
    cdef:
        int nWaves = len(waves)
        int nGal
        int *indices

        prop_set *galProps

        int nAgeList
        double *ageList
        
        double z

        int nRest = len(restFrame)
        int nObs = len(obsBands)
        int nFilter = nRest + nObs
        double *filters = NULL
        int mAB

        double *absorption = NULL

        dust_params *dustArgs = NULL

        float *cOutput 
        float[:] mvOutput
        float[:] mvMags

    for i in xrange(nSnap):
        snap = snapList[i]

        if sfhPath is None:
            galIndices = idxList[i]
            nGal = len(galIndices)
            indices = init_1d_int(np.asarray(galIndices, dtype = 'i4'))
            galProps = read_properties_by_progenitors(g_firstProgenitor, g_nextProgenitor, 
                                                      g_metals, g_sfr,
                                                      snap, indices, nGal)
            free(indices)
        else:
            galIndices = read_galaxy_indices(sfhPath[i])
            nGal = len(galIndices)
            galProps = read_properties_by_file(sfhPath[i])


        nAgeList = snap - snapMin + 1
        ageList= init_1d_double(get_age_list(fname, snap, nAgeList, h))
        z = meraxes.io.grab_redshift(fname, snap)

        if dustParams is not None:
            dustArgs = dust_parameters(dustParams[i])
        if nFilter != 0:
            filters = init_1d_double(read_filters(restFrame, obsBands, z))
            mAB = 1
        else:
            nFilter = nWaves
            if obsFrame:
                nObs = 1
            mAB = 0

        if IGM == 'I2014':
            absorption = init_1d_double(Lyman_absorption_Inoue((1. + z)*waves, z))


        cOutput = composite_spectra_cext(galProps, nGal,
                                         z, ageList, nAgeList,
                                         filters, nRest, nObs, mAB,
                                         absorption, dustArgs)
        mvOutput = <float[:nGal*nFilter]>cOutput
        output = np.asarray(mvOutput, dtype = 'f4').reshape(nGal, -1)

        if filters != NULL and nObs > 0:
            # Convert apparent magnitudes to absolute magnitudes
            output[:, nRest:] += cosmo.distmod(z)
        
        if filters == NULL and obsFrame:
            # Convert to observed frame fluxes
            factor = 10./cosmo.luminosity_distance(z).to(u.parsec).value
            output *= factor*factor
  
        if filters != NULL:
            names = []
            for i in xrange(nRest):
                names.append("M%d"%restFrame[i][0])
            for i in xrange(nObs):
                names.append(obsBands[i][0])
        elif obsFrame:
            names = (1. + z)*waves
        else:
            names = waves
       
        DataFrame(output, index = galIndices, columns = names).\
        to_hdf(get_output_name(prefix, ".hdf5", snap, path), "w")
       
        if len(snapList) == 1:
            mvMags = np.zeros(nGal*nFilter, dtype = 'f4')
            mvMags[...] = mvOutput
            mags = np.asarray(mvMags, dtype = 'f4').reshape(nGal, -1)

        free_int_spectra()
        free(ageList)
        free(dustArgs)
        free(filters)
        free(absorption)
        free(cOutput)

    free_raw_spectra()
    if sfhPath is None:
        free_meraxes(snapMin, snapMax)

    if len(snapList) == 1:
        return mags


def UV_slope(fname, snapList, idxList, h,
             sfhPath = None,
             dustParams = None,
             prefix = "slope", path = "./"):

    cdef:
        int i, iG
        int snap, nSnap
        int sanpMin, snapMax

    if isscalar(snapList):
        snapMax = snapList
        nSnap = 1
        snapList = [snapList]
        idxList = [idxList]
        if sfhPath is not None:
            sfhPath = [sfhPath]
    else:
        snapMax = max(snapList)
        nSnap = len(snapList)

    if sfhPath is None:
        snapMin = read_meraxes(fname, snapMax, h)
    else:
        snapMin = 1

    waves = get_wavelength()
    centreWaves, betaFilters = beta_filters()
    cdef:
        int nWaves = len(waves)
        int nGal
        int *indices

        prop_set *galProps

        int nAgeList
        double *ageList
        
        double z

        double *logWaves = init_1d_double(np.log(centreWaves))
        double *filters = init_1d_double(betaFilters)
        int nFilter = len(centreWaves)

        dust_params *dustArgs = NULL

        int nR = 3
   
        float *cOutput 
        float[:] mvOutput
        float[:] mvMags

    for i in xrange(nSnap):
        snap = snapList[i]

        if sfhPath is None:
            galIndices = idxList[i]
            nGal = len(galIndices)
            indices = init_1d_int(np.asarray(galIndices, dtype = 'i4'))
            galProps = read_properties_by_progenitors(g_firstProgenitor, g_nextProgenitor, 
                                                      g_metals, g_sfr,
                                                      snap, indices, nGal)
            free(indices)
        else:
            galIndices = read_galaxy_indices(sfhPath[i])
            nGal = len(galIndices)
            galProps = read_properties_by_file(sfhPath[i])

        nAgeList = snap - snapMin + 1
        ageList= init_1d_double(get_age_list(fname, snap, nAgeList, h))
        z = meraxes.io.grab_redshift(fname, snap)

        if dustParams is not None:
            dustArgs = dust_parameters(dustParams[i])

        cOutput = UV_slope_cext(galProps, nGal,
                                z, ageList, nAgeList,
                                logWaves, filters, nFilter,
                                dustArgs)

        mvOutput = <float[:nGal*(nFilter + nR)]>cOutput
        output = np.hstack([np.asarray(mvOutput[nGal*nFilter:], 
                                       dtype = 'f4').reshape(nGal, -1),
                            np.asarray(mvOutput[:nGal*nFilter], 
                                       dtype = 'f4').reshape(nGal, -1)])
        
        columns = np.append(["beta", "norm", "R"], centreWaves)
        columns[-1] = "M1600"
        DataFrame(output, index = galIndices, columns = columns). \
        to_hdf(get_output_name(prefix, ".hdf5", snap, path), "w")
        
        if len(snapList) == 1:
            mvMags = np.zeros(nGal*(nFilter + nR), dtype = 'f4')
            mvMags[...] = mvOutput
            mags = np.asarray(mvMags, dtype = 'f4').reshape(nGal, -1)

        free_int_spectra()
        free(dustArgs)
        free(ageList)
        free(cOutput)
       
    free_raw_spectra()
    free(filters)
    if sfhPath is None:
        free_meraxes(snapMin, snapMax)

    if len(snapList) == 1:
        return mags


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


@boundscheck(False)
@wraparound(False)
def dust_extinction(M1600, double z, double scatter):
    """
    Calculate the dust extinction at rest frame 1600 angstrom

    M1600: rest frame 1600 angstrom magnitudes. It can be an array.
    z: redshift

    Returns: dust extinction at rest frame 1600 angstrom
             M1600_obs = M1600 + A1600,
             where M1600_obs is the dust attenuated magnitudes
    """
    # Reference Mason et al. 2015, equation 4
    #           Bouwens 2014 et al. 2014, Table 3

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
    """
    Function of the reddening curve of Calzetti et al. 2000
    
    lam: wavelengths in a unit of angstrom
    """
    # Reference Calzetti et al. 2000, Liu et al. 2016
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
    Function to add reddening

    waves: wavelength in a unit of angstrom
    M1600: rest frame 1600 angstrom magnitudes.
    z: redshift

    Returns: the output can be directly added to intrinsic magnitudes
    """
    # waves must be in a unit of angstrom
    # The reddening curve is normalised by the value at 1600 A
    A1600 = dust_extinction(M1600, z, scatter)
    if isscalar(waves):
        return reddening_curve(waves)/reddening_curve(1600.)*A1600
    else:
        waves = np.asarray(waves)
        return reddening_curve(waves)/reddening_curve(1600.)*A1600.reshape(-1, 1)


