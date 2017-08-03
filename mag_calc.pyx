from time import time
from struct import unpack

import numpy as np
from pandas import DataFrame
from numpy import isnan, isscalar
from astropy.cosmology import FlatLambdaCDM
from dragons import meraxes

from libc.stdlib cimport malloc, free

global sTime

cdef:
    int **pFirstProgenitor
    int **pNextProgenitor
    float **pMetals
    float **pSFR


def timing_start(text):
    global sTime
    sTime = time()
    print "#*******************************************************************************"
    print text
 

def timing_end():
    global sTime
    elapsedTime = int(time() - sTime)
    print "# Done!"
    print "# Elapsed time: %i min %i sec"%(elapsedTime/60, elapsedTime%60)
    print "#*******************************************************************************\n"


def get_wavelength():
    fp = open("Input/sed_waves.bin", "rb")
    nWaves = unpack('i', fp.read(sizeof(int)))[0]
    waves = np.array(unpack('%dd'%nWaves, fp.read(nWaves*sizeof(double))))
    fp.close()
    return waves


def read_filters(restFrame, obsBands, z):
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
        filters[i] /= waves*np.log(upper/lower)
    for i in xrange(nObs):
        fWaves, trans = obsBands[i][1]
        filters[nRest + i] = np.interp(obsWaves, fWaves, trans, left = 0., right = 0.)
        filters[nRest + i] /= obsWaves*np.trapz(trans/fWaves, fWaves)
    return filters.flatten()


def HST_filters(filterNames):
    filterList = {"B435":"Input/HST_ACS_F435W.npy", 
                  "V606":"Input/HST_ACS_F606W.npy", 
                  "i775":"Input/HST_ACS_F775W.npy", 
                  "I814":"Input/HST_ACS_F814W.npy", 
                  "z850":"Input/HST_ACS_F850LP.npy", 
                  "Y098":"Input/HST_IR_F098M.npy", 
                  "Y105":"Input/HST_IR_F105W.npy", 
                  "J125":"Input/HST_IR_F125W.npy", 
                  "H160":"Input/HST_IR_F160W.npy", 
                  "3.6":"Input/HST_IRAC_3.6.npy"}
    obsBands = []
    for name in filterNames:
        obsBands.append([name, np.load(filterList[name])])
    return obsBands


def read_meraxes(fname, int snapMax, h):
    timing_start("# Read meraxes output")
    cdef:
        int snapNum = snapMax+ 1
        int snapMin = snapMax
        int snap, N
        int[:] intMemview1, intMemview2
        float[:] floatMemview1, floatMemview2
    global pFirstProgenitor 
    global pNextProgenitor
    global pMetals
    global pSFR
    pFirstProgenitor = <int**>malloc(snapNum*sizeof(int*))
    pNextProgenitor = <int**>malloc(snapMax*sizeof(int*))
    pMetals = <float**>malloc(snapNum*sizeof(float*))
    pSFR = <float**>malloc(snapNum*sizeof(float*))

    meraxes.set_little_h(h = h)
    for snap in xrange(snapMax, -1, -1):
        try:
            # Copy metallicity and star formation rate to the pointers
            gals = meraxes.io.read_gals(fname, snap, ["StellarMass", "MetalsStellarMass", "Sfr"])
            N = len(gals)
            metals = gals["MetalsStellarMass"]/gals["StellarMass"]
            metals[isnan(metals)] = 0.001
            pMetals[snap] = <float*>malloc(N*sizeof(float))
            floatMemview1 = <float[:N]>pMetals[snap]
            floatMemview2 = metals
            floatMemview1[...] = floatMemview2
            pSFR[snap] = <float*>malloc(N*sizeof(float))
            floatMemview1 = <float[:N]>pSFR[snap]
            floatMemview2 = gals["Sfr"]
            floatMemview1[...] = floatMemview2
            snapMin = snap
            gals = None
        except IndexError:
            print "# There is no galaxies in snapshot %d"%snap
            break;
    print "# snapMin = %d"%snapMin
    for snap in xrange(snapMin, snapNum):
        # Copy first progenitor indices to the pointer
        firstProgenitor = meraxes.io.read_firstprogenitor_indices(fname, snap)
        N = len(firstProgenitor)
        pFirstProgenitor[snap] = <int*>malloc(N*sizeof(int))
        intMemview1 = <int[:N]>pFirstProgenitor[snap]
        intMemview2 = firstProgenitor
        intMemview1[...] = intMemview2
        firstProgenitor = None
        # Copy next progenitor indices to the pointer
        if snap < snapMax:
            nextProgenitor = meraxes.io.read_nextprogenitor_indices(fname, snap)       
            pNextProgenitor[snap] = <int*>malloc(N*sizeof(int))
            intMemview1 = <int[:N]>pNextProgenitor[snap]
            intMemview2 = nextProgenitor
            intMemview1[...] = intMemview2
            nextProgenitor = None

    timing_end()    
    return snapMin


cdef void free_meraxes():
    free(pFirstProgenitor)
    free(pNextProgenitor)
    free(pMetals)
    free(pSFR)


cdef extern from "mag_calc_cext.h":
    void galaxy_spectra_cext(double *pOutput, 
                             double z, int snap,
                             int *indices, int nGal,
                             double *ageList, int nAgeList,
                             int **pFP, int **pNP, float **pM, float **pS) 

    void galaxy_mags_cext(float *pOutput, 
                          double z, int snap,
                          int *indices, int nGal,
                          double *ageList, int nAgeList,
                          double *filters, int nRest, int nObs,
                          int **pFP, int **pNP, float **pM, float **pS) 


def get_age_list(fname, snap, nAgeList, h):
    travelTime = meraxes.io.read_snaplist(fname, h)[2]*1e6 # Convert Myr to yr
    ageList = np.zeros(nAgeList)
    for i in xrange(nAgeList):
        ageList[i] = travelTime[snap - i - 1] - travelTime[snap]
    return ageList


def galaxy_spectra(fname, snap, indices, h):
    cdef:
        int nAgeList = snap - read_meraxes(fname, snap, h) + 1
        double *cAgeList = <double*>malloc(nAgeList*sizeof(double))
        double[:] mvCAgeList = <double[:nAgeList]>cAgeList
        double[:] mvAgeList = get_age_list(fname, snap, nAgeList, h)
    mvCAgeList[...]= mvAgeList
 
    cdef: 
        int nGal = len(indices)
        int *cIndices = <int*>malloc(nGal*sizeof(int)) 
        int[:] mvCIndices = <int[:nGal]>cIndices
        int[:] mvIndices = np.array(indices, dtype = 'i4')
    mvCIndices[...] = mvIndices

    waves = get_wavelength()
    nWaves = len(waves)
    cdef:
        double *pOutput = <double*>malloc(nGal*nWaves*sizeof(double))
        double z = meraxes.io.grab_redshift(fname, snap)

    galaxy_spectra_cext(pOutput, 
                        z, snap,
                        cIndices, nGal,
                        cAgeList, nAgeList,
                        pFirstProgenitor, pNextProgenitor, pMetals, pSFR)

    cdef:
        double[:] mvOutput = <double[:nGal*nWaves]>pOutput
        double[:] mvSpectra = np.zeros(nGal*nWaves, dtype = 'f8')
    mvSpectra[...] = mvOutput
    
    free(cAgeList)
    free(cIndices)
    free(pOutput)
    free_meraxes()

    return np.vstack([waves, np.asarray(mvSpectra).reshape(-1, nWaves)])


def galaxy_mags(fname, snapList, idxList, h, Om0, 
                restFrame = [[1600., 100.]], obsBands = [], 
                path = ""):
    cosmo = FlatLambdaCDM(H0 = 100.*h, Om0 = Om0)
    waves = get_wavelength()
    nWaves = len(waves)
   
    cdef:
        int i
        int snap, nSnap
        int sanpMin, snapMax

    if isscalar(snapList):
        snapMax = snapList
        nSnap = 1
        snapList = [snapList]
        idxList = [idxList]
    else:
        snapMax = max(snapList)
        nSnap = len(snapList)

    snapMin = read_meraxes(fname, snapMax, h)

    cdef:
        int nRest = len(restFrame)
        int nObs = len(obsBands)
        double z
        double *cFilters = <double*>malloc((nRest + nObs)*nWaves*sizeof(double))
        double[:] mvCFilter = <double[:(nRest + nObs)*nWaves]>cFilters
        double[:] mvFilter

        int nGal = 0
        int *cIndices
        int[:] mvCIndices
        int[:] mvIndices

        int nAgeList
        double *cAgeList
        double[:] mvCAgeList
        double[:] mvAgeList
        
        float *pOutput 
        float[:] mvOutput
        float[:] mvMags

    names = []
    for i in xrange(nRest):
        names.append("M%d"%restFrame[i][0])
    for i in xrange(nObs):
        names.append(obsBands[i][0])

    for i in xrange(nSnap):
        snap = snapList[i]
        indices = idxList[i]

        z = meraxes.io.grab_redshift(fname, snap)
        mvFilter = read_filters(restFrame, obsBands, z)
        mvCFilter[...] = mvFilter

        nGal = len(indices)
        cIndices = <int*>malloc(nGal*sizeof(int)) 
        mvCIndices = <int[:nGal]>cIndices
        mvIndices = np.array(indices, dtype = 'i4')
        mvCIndices[...] = mvIndices

        nAgeList = snap - snapMin + 1
        cAgeList = <double*>malloc(nAgeList*sizeof(double))
        mvCAgeList = <double[:nAgeList]>cAgeList
        mvAgeList = get_age_list(fname, snap, nAgeList, h)
        mvCAgeList[...]= mvAgeList

        pOutput = <float*>malloc(nGal*(nRest + nObs)*sizeof(float))
        mvOutput = <float[:nGal*(nRest + nObs)]>pOutput

        galaxy_mags_cext(pOutput, 
                         z, snap,
                         cIndices, nGal,
                         cAgeList, nAgeList,
                         cFilters, nRest, nObs,
                         pFirstProgenitor, pNextProgenitor, pMetals, pSFR) 

        output = np.asarray(mvOutput, dtype = 'f4').reshape(nGal, -1)
        output[:, nRest:] += cosmo.distmod(z).value
        DataFrame(output, index = indices, columns = names).to_hdf(path + "mags%03d.hdf5"%snap, "w")

        if len(snapList) == 1:
            mvMags = np.zeros(nGal*(nRest + nObs), dtype = 'f4')
            mvMags[...] = mvOutput
            mags = np.asarray(mvMags, dtype = 'f4').reshape(nGal, -1)

        free(cIndices)       
        free(cAgeList)
        free(pOutput)

    free(cFilters)
    free_meraxes()

    if len(snapList) == 1:
        return mags

  
    
