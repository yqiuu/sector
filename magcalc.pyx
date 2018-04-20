import os, sys
from warnings import warn
from time import time
from copy import deepcopy

from libc.stdlib cimport malloc, free
from libc.stdio cimport *
from libc.string cimport memcpy
from libc.math cimport exp, log

import numpy as np
from numpy import isnan, isscalar, vectorize
from pandas import DataFrame

from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from dragons import meraxes

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Basic functions                                                               #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
global sTime

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
    print "#***********************************************************"
    print text
 

def timing_end():
    global sTime
    elapsedTime = time() - sTime
    minute = int(elapsedTime)/60
    print "# Done!"
    print "# Elapsed time: %i min %.6f sec"%(minute, elapsedTime - minute*60)
    print "#***********************************************************\n"


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Functions to load galaxy properties                                           #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
cdef extern from "mag_calc_cext.h":
    struct ssp:
        short index
        float metals
        float sfr

    struct csp:
        ssp *bursts
        int nBurst

    struct gal_params:
        double z
        int nAgeStep
        double *ageStep
        int nGal
        int *indices
        csp *histories


cdef void free_gal_params(gal_params *galParams):
    cdef:
        int iG
        int nGal = galParams.nGal
    for iG in xrange(nGal):
        free(galParams.histories[iG].bursts)
    free(galParams.histories)
    free(galParams.ageStep)
    free(galParams.indices)


cdef void save_gal_params(gal_params *galParams, char *fname):
    cdef:
        int iA, iG
        FILE *fp

        double z = galParams.z
        int nAgeStep = galParams.nAgeStep
        double *ageStep = galParams.ageStep
        int nGal = galParams.nGal
        int *indices = galParams.indices
        csp *histories = galParams.histories

        int nBurst

    fp = fopen(fname, 'wb')
    # Write redshift
    fwrite(&z, sizeof(double), 1, fp)
    # Write ageStep
    fwrite(&nAgeStep, sizeof(int), 1, fp)
    fwrite(ageStep, sizeof(double), nAgeStep, fp)
    # Write indices
    fwrite(&nGal, sizeof(int), 1, fp)
    fwrite(indices, sizeof(int), nGal, fp)
    # Write histories
    for iG in xrange(nGal):
        nBurst = histories[iG].nBurst
        fwrite(&nBurst, sizeof(int), 1, fp)
        fwrite(histories[iG].bursts, sizeof(ssp), nBurst, fp)
    fclose(fp)


cdef void read_gal_params(gal_params *galParams, char *fname):
    cdef:
        int iG
        FILE *fp

        double z
        int nAgeStep
        double *ageStep
        int nGal
        int *indices
        csp *histories

        int nBurst

    timing_start("# Read galaxy properties")
    fp = fopen(fname, 'rb')
    if fp == NULL:
        raise IOError("Fail to open the input file")
    # Read redshift
    fread(&z, sizeof(double), 1, fp)
    # Read ageStep
    fread(&nAgeStep, sizeof(int), 1, fp)
    ageStep = <double*>malloc(nAgeStep*sizeof(double))
    fread(ageStep, sizeof(double), nAgeStep, fp)
    # Read indices
    fread(&nGal, sizeof(int), 1, fp)
    indices = <int*>malloc(nGal*sizeof(int))
    fread(indices, sizeof(int), nGal, fp)
    # Read histories
    histories = <csp*>malloc(nGal*sizeof(csp))
    pHistories = histories
    for iG in xrange(nGal):
        fread(&nBurst, sizeof(int), 1, fp)
        histories[iG].nBurst = nBurst
        histories[iG].bursts = <ssp*>malloc(nBurst*sizeof(ssp))
        fread(histories[iG].bursts, sizeof(ssp), nBurst, fp)
        pHistories += 1
    fclose(fp)

    galParams.z = z
    galParams.nAgeStep = nAgeStep
    galParams.ageStep = ageStep
    galParams.nGal = nGal
    galParams.indices = indices
    galParams.histories = histories

    timing_end()


cdef:
    int **g_firstProgenitor = NULL
    int **g_nextProgenitor = NULL
    float **g_metals = NULL
    float **g_sfr = NULL
    # >>>>> New metallicity tracer
    #float *g_dTime 
    # <<<<<


def read_meraxes(fname, int snapMax, h):
    #=====================================================================
    # This function reads meraxes output. It is called by galaxy_mags(...).
    # Meraxes output is stored by g_firstProgenitor, g_nextProgenitor, g_metals
    # and g_sfr. They are external variables of mag_calc_cext.c
    #
    # fname: path of the meraxes output
    # snapMax: start snapshot
    # h: liitle h
    #
    # Return: the smallest snapshot number that contains a galaxy
    #=====================================================================
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
    # >>>>> New metallicity tracer
    #global g_dTime
    # <<<<<
    timing_start("# Read meraxes output")
    g_firstProgenitor = <int**>malloc(snapNum*sizeof(int*))
    g_nextProgenitor = <int**>malloc(snapMax*sizeof(int*))
    # Unit: 1e10 M_sun (New metallicity tracer)
    g_metals = <float**>malloc(snapNum*sizeof(float*))
    # Unit: M_sun/yr
    g_sfr = <float**>malloc(snapNum*sizeof(float*))
    # >>>>>  New metallicity tracer
    # Unit: Myr
    #g_dTime = init_1d_float(np.append([0], -np.diff(meraxes.io.read_snaplist(fname, h)[2])) \
    #                        .astype('f4'))
    # <<<<<
    meraxes.set_little_h(h = h)
    for snap in xrange(snapMax, -1, -1):
        try:
            # Copy metallicity and star formation rate to the pointers
            gals = meraxes.io.read_gals(fname, snap, 
                                        props = ["ColdGas", "MetalsColdGas", "Sfr"])
            print ''
            # <<<<< Old Metallicity tracer
            metals = gals["MetalsColdGas"]/gals["ColdGas"]
            metals[isnan(metals)] = 0.001
            g_metals[snap] = init_1d_float(metals)
            # >>>>> New metallicity tracer
            #g_metals[snap] = init_1d_float(gals["MetalsStellarMass"])
            # <<<<<
            g_sfr[snap] = init_1d_float(gals["Sfr"])
            snapMin = snap
            gals = None
        except IndexError:
            print "# No galaxies in snapshot %d"%snap
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
    #=====================================================================
    # Function to free g_firstProgenitor, g_nextProgenitor, 
    # g_metals, and g_sfr
    #=====================================================================
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
    # >>>>>  New metallicity tracer
    #free(g_dTime) 
    # <<<<<


cdef struct trace_params:
    int **firstProgenitor
    int **nextProgenitor
    # Unit: 1e10 M_sum (New metallcitiy tracer)
    float **metals
    # Unit: 1 M_sun/yr
    float **sfr
    # Unit: 1 Myr
    float *dTime
    int tSnap
    ssp *bursts
    int nBurst


DEF MAX_NODE = 100000

cdef void trace_progenitors(int snap, int galIdx, trace_params *args):
    cdef:
        float sfr
        ssp *pBursts
        int nProg
    if galIdx >= 0:
        sfr = args.sfr[snap][galIdx]
        if sfr > 0.:
            args.nBurst += 1
            nProg = args.nBurst
            if (nProg >= MAX_NODE):
                raise MemoryError("Number of progenitors exceeds MAX_NODE")
            pBursts = args.bursts + nProg
            pBursts.index = args.tSnap - snap
            # <<<<< Old metallicity tracer
            pBursts.metals = args.metals[snap][galIdx]
            # >>>>> New metallicity tracer
            #pBursts.metals = trace_metallicity(snap, galIdx, args)
            # <<<<<
            pBursts.sfr = sfr
            #print "snap %d, galIdx %d, metals %.3f sfr %.3f\n"%(snap, galIdx, 
            #                                                    args.metals[snap][galIdx], 
            #                                                    sfr)
        
        trace_progenitors(snap - 1, args.firstProgenitor[snap][galIdx], args)
        trace_progenitors(snap, args.nextProgenitor[snap][galIdx], args)


cdef inline float trace_metallicity(int snap, int galIdx, trace_params *args):
    cdef:
        float progMetalsMass = 0
        int progSnap = snap - 1
        int progIdx = args.firstProgenitor[snap][galIdx]
    if progIdx < 0:
        return args.metals[snap][galIdx]/args.sfr[snap][galIdx]/args.dTime[snap]*1e4
        # The factor 1e4 is from the unit conversion
    else:
        progMetalsMass += args.metals[progSnap][progIdx]
        progIdx = args.nextProgenitor[progSnap][progIdx]
        while(progIdx > 0):
            progMetalsMass += args.metals[progSnap][progIdx]
            progIdx = args.nextProgenitor[progSnap][progIdx]
        return (args.metals[snap][galIdx] - progMetalsMass) \
               /args.sfr[snap][galIdx]/args.dTime[snap]*1e4


cdef csp *trace_properties(int tSnap, int[:] indices):
    cdef:
        int iG
        int nGal = indices.shape[0]

        size_t memSize
        size_t totalMemSize = 0

        csp *histories = <csp*>malloc(nGal*sizeof(csp))
        csp *pHistories
        ssp bursts[MAX_NODE]
        trace_params args

        int galIdx
        int nProg
        float sfr

    args.firstProgenitor = g_firstProgenitor
    args.nextProgenitor = g_nextProgenitor
    args.metals = g_metals
    args.sfr = g_sfr
    # >>>>>  New metallicity tracer
    #args.dTime = g_dTime 
    # <<<<<
    args.tSnap = tSnap
    args.bursts = bursts

    timing_start("# Read galaxies properties")
    for iG in xrange(nGal):
        galIdx = indices[iG]
        nProg = -1
        sfr = args.sfr[tSnap][galIdx]
        if sfr > 0.:
            nProg += 1
            bursts[nProg].index = 0
            # <<<<< Old metallicity tracer
            bursts[nProg].metals = args.metals[tSnap][galIdx]
            # >>>>> New metallicity tracer
            #bursts[nProg].metals = trace_metallicity(tSnap, galIdx, &args)
            # <<<<<
            bursts[nProg].sfr = sfr
        args.nBurst = nProg
        trace_progenitors(tSnap - 1, args.firstProgenitor[tSnap][galIdx], &args)
        nProg = args.nBurst + 1
        pHistories = histories + iG
        pHistories.nBurst = nProg
        if nProg == 0:
            pHistories.bursts = NULL
            print "Warning: snapshot %d, index %d"%(tSnap, galIdx)
            print "         the star formation rate is zero throughout the histroy"
        else:
            memSize = nProg*sizeof(ssp)
            pHistories.bursts = <ssp*>malloc(memSize)
            memcpy(pHistories.bursts, bursts, memSize)
            totalMemSize += memSize

    print "# %.1f MB memory has been allocted"%(totalMemSize/1024./1024.)
    timing_end()
    return histories


cdef void *read_star_formation_history(gal_params *galParams, fname, snapshot, gals, h):
    if type(gals) is str:
        read_gal_params(galParams, gals)
    else:
        # Read redshift
        galParams.z = meraxes.io.grab_redshift(fname, snapshot)
        # Read lookback time
        galParams.nAgeStep = snapshot
        timeStep = meraxes.io.read_snaplist(fname, h)[2]*1e6 # Convert Myr to yr
        ageStep = np.zeros(snapshot, dtype = 'f8')
        for iA in xrange(snapshot):
            ageStep[iA] = timeStep[snapshot - iA - 1] - timeStep[snapshot]
        galParams.ageStep = init_1d_double(ageStep)
        gals = np.asarray(gals, dtype = 'i4')
        galParams.nGal = len(gals)
        galParams.indices = init_1d_int(gals)
        galParams.histories = trace_properties(snapshot, gals)
    return galParams


#def trace_star_formation_history(fname, snap, galIndices, h):
#    #=====================================================================
#    # Read galaxy properties from Meraxes outputs
#    #=====================================================================
#    cdef int snapMin = read_meraxes(fname, snap, h)
#    # Trace galaxy merge trees
#    cdef:
#        int iG
#        int nGal = len(galIndices)
#        int *indices = init_1d_int(np.asarray(galIndices, dtype = 'i4'))
#        csp *histories = \
#        read_gal_params_by_progenitors(snap, indices, nGal)
#    free(indices)
#    free_meraxes(snapMin, snap)
#    # Convert output to numpy array
#    cdef:
#        int iN
#        int nBurst
#        ssp *bursts
#        double[:, ::1] mvBursts
#    output = np.empty(nGal, dtype = object)
#    for iG in xrange(nGal):
#        nBurst = histories[iG].nBurst
#        bursts = histories[iG].bursts
#        mvBursts = np.zeros([nBurst, 3])
#        for iN in xrange(nBurst):
#            mvBursts[iN][0] = bursts[iN].index
#            mvBursts[iN][1] = bursts[iN].metals
#            mvBursts[iN][2] = bursts[iN].sfr
#        output[iG] = np.asarray(mvBursts)
#    return output


def save_star_formation_history(fname, snapList, idxList, h, 
                                prefix = 'sfh', outPath = './'):
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
        The name of the output file is 'prefix_XXX.hdf5', where XXX is
        number of the snapshot.
    outPath: str
        Path to the output.
    """
    cdef:
        int iS, nSnap
        int snap, snapMax, snapMin
        gal_params galParams

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
    for iS in xrange(nSnap):
        read_star_formation_history(&galParams, fname, snapList[iS], idxList[iS], h)
        save_gal_params(&galParams, get_output_name(prefix, '.bin', snapList[iS], outPath))
        free_gal_params(&galParams)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Functions to compute the IGM absorption                                       #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def Lyman_absorption_Fan(double[:] obsWaves, double z):
    #=====================================================================
    # Depreciate function. It is original to calculate the optical depth of
    # Fan et al. 2006 
    #=====================================================================
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
    #=====================================================================
    # Function to calculate the optical depth of Inoue et al. 2014
    # It is called by galaxy_mags(...).
    #
    # obsWaves: wavelength in unit of angstrom
    # z: redshift
    #
    # Return: transmission (dimensionless)
    # Reference Inoue et al. 2014
    #=====================================================================
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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Functions about ISM absorptionb                                               #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
cdef extern from "mag_calc_cext.h":
    struct dust_params:
        double tauUV_ISM
        double nISM
        double tauUV_BC
        double nBC
        double tBC


cdef dust_params *init_dust_parameters(dust):
    cdef:
        int iG
        int nGal = len(dust)
        double[:, ::1] mvDustParams = np.array(dust)
        dust_params *dustParams = <dust_params*>malloc(nGal*sizeof(dust_params))
        dust_params *pDustParams = dustParams

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
cdef extern from "mag_calc_cext.h":
    struct sed_params:
        # Raw templates
        int minZ
        int maxZ
        int nZ
        double *Z
        int nWaves
        double *waves
        int nAge
        double *age
        double *raw
        # Redshift
        double z
        # Filters
        int nFlux
        int nObs
        int *nFilterWaves
        double *filterWaves
        double *filters
        double *logWaves
        # IGM absorption
        double *LyAbsorption
        # Working templates
        int nAgeStep
        double *ageStep
        double *integrated
        double *ready
        double *working


cdef void init_templates_raw(sed_params *spectra, path):
    #=====================================================================
    # The dictionary define by *path* should contain:                               
    #                                                                               
    # "sed_Z.npy" Metallicity of SED templates in a 1-D array                         
    #                                                                               
    # "sed_waves.npy" Wavelength of SED templates in a unit of angstrom in 
    # a 1-D array                                                                         
    #                                                                               
    # "sed_age.npy" Stellar age of SED templates in a unit of yr in a 1-D 
    # array     
    #                                                                               
    # "sed_flux.npy" Flux density of SED templates in a unit of erg/s/A/cm^2 
    # in a 3-D array. The flux density should be normlised by the surface 
    # area of a 10 pc sphere. The first, second and third dimensions should 
    # be metallicity, wavelength and stellar age respectively.
    #=====================================================================
    timing_start("# Read SED templates")
    # Read metallicity range
    Z = np.load(os.path.join(path, "sed_Z.npy"))
    spectra.Z = init_1d_double(Z)
    spectra.nZ = len(Z)
    spectra.minZ = <short>(Z.min()*1000 - 0.5)
    spectra.maxZ = <short>(Z.max()*1000 - 0.5)
    print "# Metallicity range:\n#\t%.3f to %.3f"%(Z[0], Z[-1])
    # Read wavelength
    waves = np.load(os.path.join(path, "sed_waves.npy"))
    print "# Wavelength range:\n#\t%.1f AA to %.1f AA"%(waves[0], waves[-1])
    spectra.waves = init_1d_double(waves)
    spectra.nWaves = len(waves)
    # Read stellar age
    age = np.load(os.path.join(path, "sed_age.npy"))
    print "# Stellar age range:\n#\t%.2f Myr to %.2f Myr"%(age[0]*1e-6, age[-1]*1e-6)
    spectra.age = init_1d_double(age)
    spectra.nAge = len(age)
    # Read flux
    flux = np.load(os.path.join(path, "sed_flux.npy")).flatten()
    spectra.raw = init_1d_double(flux)
    timing_end()


def get_wavelength(path):
    #=====================================================================
    # Return wavelengths of SED templates in a unit of angstrom
    #=====================================================================
    return np.load(os.path.join(path, "sed_waves.npy"))


cdef void free_raw_spectra(sed_params *spectra):
    free(spectra.age)
    free(spectra.waves)
    free(spectra.raw)


cdef extern from "mag_calc_cext.h":
    void shrink_templates_raw(sed_params *spectra, double maxAge)


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
                        [2400., 2580.],
                        [1550., 1650.]])
    betaBands = np.zeros_like(windows)
    betaBands[:, 0] = (windows[:, 1] + windows[:, 0])/2.
    betaBands[:, 1] = windows[:, 1] - windows[:, 0]
    return betaBands


cdef void init_filters(sed_params *spectra, waves, int nBeta, restBands, obsBands, double z):
    #=====================================================================
    # This function is to generate transmission curves that has the
    # same wavelengths with SED templates. It is called by
    # galaxy_mags(...). The input format refer to galaxy_mags(...).
    #
    # Before integration over the filters, the fluxes must be a function
    # of wavelength.
    # After integration over the filters, the fluxex becomes a function
    # of frequency.
    #=====================================================================
    cdef:
        int iF
        int nRest = len(restBands)
        int nObs = len(obsBands)
        int nFlux = nRest + nObs
    nFilterWaves = np.zeros(nFlux, dtype = 'i4')
    filterWaves = np.empty(nFlux, dtype = object)
    filters = np.empty(nFlux, dtype = object)
    for iF in xrange(nRest):
        centre, bandWidth = restBands[iF]
        lower = centre - bandWidth/2.
        upper = centre + bandWidth/2.
        if lower < min(waves) or upper > max(waves):
            raise ValueError("Filter wavelength ranges are beyond SED templates")
        newWaves = np.array([lower])
        for w in waves:
            if w > lower and w < upper:
                newWaves = np.append(newWaves, w)
        newWaves = np.append(newWaves, upper)
        if iF < nBeta:
            newFilters = np.full(len(newWaves), 1.)
            newFilters /= np.trapz(newFilters, newWaves)
        else:
            # Transform the filter such that it can give AB magnitude
            newFilters = np.full(len(newWaves), 1./np.log(upper/lower))
            newFilters *= 3.34e4*newWaves
        nFilterWaves[iF] = len(newWaves)
        filterWaves[iF] = newWaves
        filters[iF] = newFilters
    for iF in xrange(nRest, nFlux):
        newWaves, newFilters = deepcopy(obsBands[iF - nRest][1])
        if min(newWaves) < min(waves)*(1. + z) or max(newWaves) > max(waves)*(1. + z):
            raise ValueError("Filter wavelength ranges are beyond SED templates")
        newFilters /= np.trapz(newFilters/newWaves, newWaves)
        newFilters *= 3.34e4*newWaves
        #newWaves /= (1. + z)
        nFilterWaves[iF] = len(newWaves)
        filterWaves[iF] = newWaves
        filters[iF] = newFilters
    spectra.nFlux = nFlux
    spectra.nObs = nObs
    spectra.nFilterWaves = init_1d_int(nFilterWaves)
    spectra.filterWaves = init_1d_double(np.concatenate(filterWaves).astype('f8'))
    spectra.filters = init_1d_double(np.concatenate(filters).astype('f8'))
    if nBeta > 0:
        spectra.logWaves = init_1d_double(np.log(restBands[:, 0]))
    else:
        spectra.logWaves = NULL


cdef void free_filters(sed_params *spectra):
    free(spectra.nFilterWaves)
    free(spectra.filterWaves)
    free(spectra.filters)


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


cdef extern from "mag_calc_cext.h" nogil:
    double *composite_spectra_cext(sed_params *spectra,
                                   gal_params *galParams, dust_params *dustParams,
                                   short outType, short nThread)


def composite_spectra(fname, snapList, gals, h, Om0, sedPath,
                      dust = None, IGM = 'I2014', 
                      outType = 'ph', 
                      restBands = [[1600, 100],], obsBands = [], obsFrame = False,
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

    if isscalar(snapList):
        snapMax = snapList
        nSnap = 1
        snapList = [snapList]
        gals = [gals]
    else:
        snapMax = max(snapList)
        nSnap = len(snapList)

    if type(gals[0]) is not str:
        snapMin = read_meraxes(fname, snapMax, h)

    waves = get_wavelength(sedPath)
    cdef:
        sed_params spectra
        gal_params galParams
        double z
        int nGal

        int nWaves = len(waves)
        int nRest = 0
        int nObs = 0
        int nFlux = 0
        double *logWaves= NULL
        int cOutType = 0

        int nR = 3

        dust_params *dustParams = NULL

        double *cOutput 
        double[:] mvOutput

    for iS in xrange(nSnap):
        # Read star formation rates and metallcities form galaxy merger trees
        read_star_formation_history(&galParams, fname, snapList[iS], gals[iS], h)
        z = galParams.z
        nGal = galParams.nGal
        # Set redshift
        spectra.z = z
        # Convert the format of dust parameters 
        if dust is not None:
            dustParams = init_dust_parameters(dust[iS])
        # Compute the transmission of the IGM
        if IGM == 'I2014':
            spectra.LyAbsorption = init_1d_double(Lyman_absorption_Inoue((1. + z)*waves, z))
        else:
            spectra.LyAbsorption = NULL
        # Generate Filters
        if outType == 'ph':
            nRest = len(restBands)
            nObs = len(obsBands)
            nFlux = nRest + nObs
            init_filters(&spectra, waves, 0, restBands, obsBands, z)
            cOutType = 0
        elif outType == 'sp':
            nFlux = nWaves
            spectra.nFlux = nWaves
            if obsFrame:
                nObs = nWaves
                spectra.nObs = nWaves
            else:
                spectra.nObs = 0
            spectra.nFilterWaves = NULL
            spectra.filterWaves = NULL
            spectra.filters = NULL
            spectra.logWaves = NULL
            cOutType = 1
        elif outType == 'UV slope':
            betaBands = beta_filters()
            centreWaves = betaBands[:, 0]
            nRest = len(centreWaves)
            nFlux = nRest
            init_filters(&spectra, waves, nRest - 1, betaBands, [], z)
            cOutType = 2
        else:
            raise KeyError("outType can only be 'ph', 'sp' and 'UV Slope'")
        # Read raw SED templates
        init_templates_raw(&spectra, sedPath)
        shrink_templates_raw(&spectra, galParams.ageStep[galParams.nAgeStep - 1])
        # Compute spectra
        cOutput = composite_spectra_cext(&spectra, &galParams, dustParams,
                                         cOutType, nThread)
        # Save the output to a numpy array
        if outType == 'UV slope':
            mvOutput = <double[:nGal*(nFlux + nR)]>cOutput
            output = np.hstack([np.asarray(mvOutput[nGal*nFlux:], 
                                           dtype = 'f4').reshape(nGal, -1),
                                np.asarray(mvOutput[:nGal*nFlux], 
                                           dtype = 'f4').reshape(nGal, -1)])
        else:
            mvOutput = <double[:nGal*nFlux]>cOutput
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

        free_gal_params(&galParams)
        free(dustParams)
        free(spectra.LyAbsorption)
        free_filters(&spectra)
        free(cOutput)
        free(logWaves)

    free_raw_spectra(&spectra)
    if type(gals[0]) is not str:
        free_meraxes(snapMin, snapMax)

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
        gal_params *galParams
        sed_params *spectra
        short nThread

    def __cinit__(self, sfhList, sedPath, nThread = 1):
        cdef:
            int iS
            int nSnap = len(sfhList)
            gal_params *pGalParams
            sed_params *pSpectra
            int nBeta

        self.nSnap = nSnap
        self.galParams = <gal_params*>malloc(nSnap*sizeof(gal_params))
        self.spectra = <sed_params*>malloc(nSnap*sizeof(sed_params))
        self.nThread = <short>nThread

        waves = get_wavelength(sedPath)
        betaBands = beta_filters()
        nBeta = len(betaBands) - 1

        pGalParams = self.galParams
        pSpectra = self.spectra
        for iS in xrange(nSnap):
            # Read star formation rates and metallcities form galaxy merger trees
            read_gal_params(pGalParams, sfhList[iS])
            # Set redshift
            pSpectra.z = pGalParams.z
            #
            pSpectra.LyAbsorption = NULL
            # Generate filters
            init_filters(pSpectra, waves, nBeta, betaBands, [], 0.)
            # Read raw SED templates
            init_templates_raw(pSpectra, sedPath)
            shrink_templates_raw(pSpectra, pGalParams.ageStep[pGalParams.nAgeStep - 1])
            #
            pGalParams += 1
            pSpectra += 1


    def __dealloc__(self):
        cdef:
            int iS
            int nSnap = self.nSnap
            gal_params *pGalParams = self.galParams
            sed_params *pSpectra = self.spectra
            
        for iS in xrange(nSnap):
            free_gal_params(pGalParams)
            free_filters(pSpectra)
            free_raw_spectra(pSpectra)
            #
            pGalParams += 1
            pSpectra += 1

        free(self.galParams)
        free(self.spectra)


    def run(self, dust, prefix = 'mags', outPath = './'):
        cdef:
            int iS
            int nSnap = self.nSnap
            gal_params *pGalParams = self.galParams
            sed_params *pSpectra = self.spectra

            int nGal
            dust_params *dustParams = NULL

            int nFlux = pSpectra.nFlux
            int nR = 3
            int[:] mvIndices
            double *cOutput
            double[:] mvOutput

        betaBands = beta_filters()
        centreWaves = betaBands[:, 0]

        snapList = [100, 78, 63]
        for iS in xrange(nSnap):
            # Compute spectra
            nGal = pGalParams.nGal
            dustParams = init_dust_parameters(dust[iS])
            cOutput = composite_spectra_cext(pSpectra, pGalParams, dustParams, 2, self.nThread)
            free(dustParams)
            mvOutput = <double[:nGal*(nFlux + nR)]>cOutput
            output = np.hstack([np.asarray(mvOutput[nGal*nFlux:], 
                                           dtype = 'f4').reshape(nGal, -1),
                                np.asarray(mvOutput[:nGal*nFlux], 
                                           dtype = 'f4').reshape(nGal, -1)])
            columns = []
            columns = np.append(["beta", "norm", "R"], centreWaves)
            columns[-1] = "M1600-100"
            indices = np.asarray(<int[:nGal]>pGalParams.indices, dtype = 'i4')
            DataFrame(output, index = indices, columns = columns).\
            to_hdf(get_output_name(prefix, ".hdf5", snapList[iS], outPath), "w")

            pGalParams += 1
            pSpectra += 1


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Dust model of Mason et al . 2015                                              #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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


