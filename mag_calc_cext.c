#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

//#define SURFACE_AREA 1.1965e40 // 4*pi*(10 pc)**2 unit cm^2
//#define JANSKY(x) (3.34e4*(x)*(x))
#define M_AB(x) (-2.5*log10(x) + 8.9) // Convert Jansky to AB magnitude
#define TOL 1e-30 // Minimum Flux
short g_nThread = 1;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Profiling functions                                                         *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifdef TIMING
    #define MAX_BLOCK 100
    #define INTEGRATION 0
    #define DUST 1
    #define WORKING1 2
    #define WORKING2 3

    static double timer[MAX_BLOCK];
    static double counter[MAX_BLOCK];
    static char blockNames[MAX_BLOCK][128];
    static struct timespec g_sTime;
    static struct timespec g_eTime;
    static struct timespec g_sTime2;
    static struct timespec g_eTime2;


    void timing_start(char* text) {
        clock_gettime(CLOCK_REALTIME, &g_sTime);
        printf("#***********************************************************\n");
        printf("# %s\n", text);
    }


    void timing_end(void) {
        clock_gettime(CLOCK_REALTIME, &g_eTime);
        double elapsedTime = g_eTime.tv_sec - g_sTime.tv_sec \
                             + (g_eTime.tv_nsec - g_sTime.tv_nsec)/1e9;
        int minute = (int)elapsedTime/60;
        printf("# 100.0%% complete!\n");
        printf("# Done!\n");
        printf("# Elapsed time: %d min %.6f sec\n", minute, elapsedTime - minute*60);
        printf("#***********************************************************\n\n");
    }


    void init_profiler(void) {
        int iB;
        for(iB = 0; iB < MAX_BLOCK; ++iB) {
            timer[iB] = 0.;
            counter[iB] = 0;
        }
    }


    void profiler_start(char* name, int blockIdx) {
        clock_gettime(CLOCK_REALTIME, &g_sTime2);
        strcpy(blockNames[blockIdx], name);
    }


    void profiler_end(int blockIdx) {
        clock_gettime(CLOCK_REALTIME, &g_eTime2);
        timer[blockIdx] += g_eTime2.tv_sec - g_sTime2.tv_sec \
                           + (g_eTime2.tv_nsec - g_sTime2.tv_nsec)/1e9;
        counter[blockIdx] += 1;
    }


    void profiler_summary(void) {
        int iB, ncall;
        printf("#***********************************************************\n");
        for(iB = 0; iB < MAX_BLOCK; ++iB) {
            ncall = counter[iB];
            if (ncall == 0)
                continue;
            printf("# %s\n", blockNames[iB]);
            printf("#\tncall: %d\taverage time: %.6f ms\n", ncall, timer[iB]/ncall*1e3);
        }
        printf("#***********************************************************\n\n");
    }


    inline void report(int i, int tot) {
        int n = tot > 10 ? tot/10 : 1;
        if (i%n == 0) {
            printf("# %5.1f%% complete!\r", 100.*(i + 1)/tot);      
            fflush(stdout);
        }
    }
#endif


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Basic functions                                                             *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
FILE *open_file(char *fName, char *mode) {
/* Open the file with specific mode */
    FILE *fp;
    if ((fp = fopen(fName, mode)) == NULL) {
        printf("File open error: \"%s\"!\n", fName);
        exit(0);
    }
    printf("# File opened: \"%s\"!\n", fName);
    return fp;
}


double **malloc_2d_double(int nRow, int nCol) {
    int i;
    double **target = malloc(nRow*sizeof(double*));
    for(i = 0; i < nRow; ++i) 
        target[i] = (double*)malloc(nCol*sizeof(double));
    return target;
}


double **memcpy_2d_double(double **source, int nRow, int nCol) {
    int i;
    double **target = malloc_2d_double(nRow, nCol);
    for(i = 0; i < nRow; ++i) 
        memcpy(target[i], source[i], nCol*sizeof(double));
    return target;
}


void free_2d_double(double **p, int nRow) {
    int i;
    for(i = 0; i < nRow; ++i)
        free(p[i]);
    free(p);
}


inline unsigned int bisection_search(double a, double *x, int nX) {
    /* return idx such x[idx] <= a < x[idx + 1] 
     * a must be x[0] <= a < x[nX - 1]
     */
    unsigned int idx0 = 0;
    unsigned int idx1 = nX - 1;
    unsigned int idxMid;
    while(idx1 - idx0 > 1) {
        idxMid = (idx0 + idx1)/2;
        if(a >= x[idxMid])
            idx0 = idxMid;
        else if(a < x[idxMid])
            idx1 = idxMid;
    }
    return idx0;
}


inline double interp(double xp, double *x, double *y, int nPts) {
    /* Interpolate a given points */
    unsigned int idx0, idx1;
    if((xp < x[0]) || (xp > x[nPts - 1])) {
        printf("Error: Point %10.5e is beyond the interpolation region\n", xp);
        exit(0);
    }
    if (xp == x[nPts - 1])
        return y[nPts - 1];
    else {
        idx0 = bisection_search(xp, x, nPts);
        if (x[idx0] == xp)
            return y[idx0];
        idx1 = idx0 + 1;
        return y[idx0] + (y[idx1] - y[idx0])*(xp - x[idx0])/(x[idx1] - x[idx0]);
    }
}


inline double trapz_table(double *y, double *x, int nPts, double a, double b) {
    /* Integrate tabular data from a to b */
    int i;
    int idx0, idx1;
    double ya, yb;
    double I;
    if (x[0] > a) {
        printf("Error: Integration range %10.5e is beyond the tabular data\n", a);
        exit(0);
    }
    if (x[nPts - 1] < b) {
        printf("Error: Integration range %10.5e is beyond the tabular data\n", b); 
        exit(0);
    }
    if (a > b) {
        printf("Error: a must be smaller than b\n");
        exit(0);
    }
    idx0 = bisection_search(a, x, nPts);
    idx1 = idx0 + 1;

    ya = y[idx0] + (y[idx1] - y[idx0])*(a - x[idx0])/(x[idx1] - x[idx0]);
    if(b <= x[idx1]) {
        yb = y[idx0] + (y[idx1] - y[idx0])*(b - x[idx0])/(x[idx1] - x[idx0]);
        return (b - a)*(yb + ya)/2.;
    }
    else 
        I = (x[idx1] - a)*(y[idx1] + ya)/2.;

    for(i = idx1; i < nPts - 1; ++i) {
        if (x[i + 1] < b)
            I += (x[i + 1] - x[i])*(y[i + 1] + y[i])/2.;
        else if (x[i] < b) {
            yb = y[i] + (y[i + 1] - y[i])*(b - x[i])/(x[i + 1] - x[i]);
            I += (b - x[i])*(yb + y[i])/2.;
        }
        else
            break;
    }
    return I;
}


inline double trapz_filter(double *filter, double *flux, double *waves, int nWaves) {
    /* integrate the flux in a filter */
    int i;
    double y0 = filter[0]*flux[0];
    double y1;
    double I = 0.;
    for(i = 1; i < nWaves; ++i) {
        y1 = filter[i]*flux[i];
        I += (waves[i] - waves[i - 1])*(y0 + y1);
        y0 = y1;
    }
    return I/2.;
}


struct linResult {
    double slope;
    double intercept;
    double R;
};


inline struct linResult linregress(double *x, double *y, int nPts) {
    int i;
   
    double xSum = 0.;
    for(i = 0; i < nPts; ++i)
        xSum += x[i];

    double ySum = 0.;
    for(i = 0; i < nPts; ++i)
        ySum += y[i];

    double xxSum = 0.;
    for(i = 0; i < nPts; ++i)
        xxSum += x[i]*x[i];

    double xySum = 0.;
    for(i = 0; i < nPts; ++i)
        xySum += x[i]*y[i];

    double denominator = nPts*xxSum - xSum*xSum;
    
    double slope = (nPts*xySum - xSum*ySum)/denominator;
    double intercept = (xxSum*ySum - xSum*xySum)/denominator;

    double yReg;
    double delta;
    double ssRes = 0.;
    for(i = 0; i < nPts; ++i) {
        yReg = slope*x[i] + intercept;
        delta = yReg - y[i];
        ssRes += delta*delta;
    }

    double yMean = ySum/nPts;
    double ssTot = 0.;
    for(i = 0; i < nPts; ++i) {
        delta = yMean - y[i];
        ssTot += delta*delta;
    }

    double R = sqrt(1. - ssRes/ssTot);

    struct linResult result;
    result.slope = slope;
    result.intercept = intercept;
    result.R = R;
    return result;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Struct to store galaxy properites                                           *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
struct ssp {
    short index;
    float metals;
    float sfr;
};

struct csp {
    struct ssp *bursts;
    int nBurst;
};

struct gal_params {
    double z;
    int nAgeStep;
    double *ageStep;
    int nGal;
    int *indices;
    struct csp *histories;
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Functions to process SEDs                                                   *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
// Struct for SED templates
struct sed_params {
    // Raw templates
    int minZ;
    int maxZ;
    int nZ;
    double *Z;
    int nWaves;
    double *waves;
    int nAge;
    double *age;
    double *raw;
    // Filters
    int nFlux;
    int nObs;
    int *nFilterWaves;
    double *filterWaves;
    double *filters;
    double *logWaves;
    // IGM absoprtion
    double *LyAbsorption;
    // Working templates
    int nAgeStep;
    double *ageStep;
    double *integrated;
    double *ready;
    double *working;
};


void init_templates_working(struct sed_params *spectra, 
                            double *ageStep, int nAgeStep, int nFlux) {
    size_t intFluxSize = spectra->nZ*nAgeStep*spectra->nWaves*sizeof(double);
    size_t workingSize = (spectra->maxZ + 1)*nAgeStep*nFlux*sizeof(double);
    spectra->ageStep = ageStep;
    spectra->nAgeStep = nAgeStep;
    spectra->integrated = (double*)malloc(intFluxSize);
    spectra->ready = (double*)malloc(intFluxSize);
    spectra->working = (double*)malloc(workingSize);
}


void free_templates_working(struct sed_params *spectra) {
    free(spectra->integrated);
    free(spectra->ready);
    free(spectra->working);
}


void shrink_templates_raw(struct sed_params *spectra, double maxAge, double z) {
    if (spectra->filters == NULL)
        return;

    int iA, iF, iW, iZ;

    int nZ = spectra->nZ;
    int nWaves = spectra->nWaves;
    double *waves = spectra->waves;
    int nAge = spectra->nAge;
    double *age = spectra->age;
    double *raw = spectra->raw;
    int nFlux = spectra->nFlux;
    int nRest = nFlux - spectra->nObs;
    int *nFilterWaves = spectra->nFilterWaves;
    double *filterWaves = spectra->filterWaves;
    double *LyAbsorption = spectra->LyAbsorption;
    int nFW;
    double *pFilterWaves;

    int inFlag = 0;
    int outFlag;
    int nNewWaves = 0;
    int *wavesIndices = malloc(nWaves*sizeof(int));
    double w;
    int nNewAge = 0;
    double *newWaves;
    double *newRaw;
    double *pNewRaw;
    double *newAbsorption;

    // Find shrinked wavelength ranges
    printf("#***********************************************************\n");
    printf("# Shrinked wavelength ranges:\n");
    for(iW = 1; iW < nWaves; ++iW) {
        w = waves[iW];
        pFilterWaves = filterWaves;
        outFlag = 1;
        for(iF = 0; iF < nFlux; ++iF) {
            nFW = nFilterWaves[iF];
            if (iF == nRest)
                w *= 1. + z;
            if (w > pFilterWaves[0] && w < pFilterWaves[nFW - 1]) {
                outFlag = 0;
                if (!inFlag) {
                    printf("#\t%.1f AA to ", waves[iW - 1]);
                    wavesIndices[nNewWaves++] = iW - 1;
                    wavesIndices[nNewWaves++] = iW;
                    inFlag = 1;
                }
                else {
                    wavesIndices[nNewWaves++] = iW;
                }
                break;
            }
            pFilterWaves += nFW;
        }
        if (inFlag && outFlag) {
            printf("%.1f AA\n", waves[iW]);
            wavesIndices[nNewWaves++] = iW;
            inFlag = 0;
        }
    }
    printf("# Original nWaves: %d\n", nWaves);
    printf("# Shrinked nWaves: %d\n", nNewWaves);
    // Find nNewAge
    printf("#\n# Shrinked age range:\n");
    for(iA = 0; age[iA] < maxAge; ++iA);
    nNewAge = iA + 1;
    spectra->nAge = nNewAge;
    printf("#\t%.2f Myr to %.2f Myr\n", age[0]*1e-6, age[iA]*1e-6);
    printf("# Original nAge: %d\n", nAge);
    printf("# Shrinked nAge: %d\n", nNewAge);
    printf("#***********************************************************\n\n");
    // Construct new wavelengths
    newWaves = (double*)malloc(nNewWaves*sizeof(double));
    for(iW = 0; iW < nNewWaves; ++iW) 
        newWaves[iW] = waves[wavesIndices[iW]];
    spectra->nWaves = nNewWaves;
    spectra->waves = newWaves;
    free(waves);
    // Construct new raw templates
    newRaw = (double*)malloc(nZ*nNewWaves*nNewAge*sizeof(double));
    pNewRaw = newRaw;
    for(iZ = 0; iZ < nZ; ++iZ)
        for(iW = 0; iW < nNewWaves; ++iW)
            for(iA = 0; iA < nNewAge; ++iA)
                *pNewRaw++ = raw[(iZ*nWaves + wavesIndices[iW])*nAge + iA];
    spectra->raw = newRaw;
    free(raw);
    // Construct new IGM absorption
    if (LyAbsorption != NULL) {
        newAbsorption = (double*)malloc(nNewWaves*sizeof(double));
        for(iW = 0; iW < nNewWaves; ++iW)
            newAbsorption[iW] = LyAbsorption[wavesIndices[iW]];
        spectra->LyAbsorption = newAbsorption;
        free(LyAbsorption);
    }
    free(wavesIndices);
}


void integrate_templates_raw(struct sed_params *spectra) {
    #ifdef TIMING
        profiler_start("Integration over time", INTEGRATION);
    #endif
    int iA, iW, iZ;
    double *pData;

    int nAgeStep = spectra->nAgeStep;
    double *ageStep = spectra->ageStep;
    int nAge;
    double *age;
    int nWaves; 
    int nZ;
    double *data;
    // Spectra after integration over time
    // The first dimension refers to metallicites and ages
    // The last dimension refers to wavelengths
    double *intData;

    nAge = spectra->nAge;
    age = spectra->age;
    nWaves = spectra->nWaves; 
    nZ = spectra->nZ;
    data = spectra->raw;
    intData = spectra->integrated;
    for(iZ = 0; iZ < nZ; ++iZ) 
        for(iA = 0; iA < nAgeStep; ++iA) {
            pData = intData + (iZ*nAgeStep + iA)*nWaves;
            for(iW = 0; iW < nWaves; ++iW) {
                if (iA == 0) 
                    // The first time step of SED templates is typicall not zero
                    // Here assumes that the templates is zero beween zero
                    // and the first time step
                    pData[iW] = trapz_table(data + (iZ*nWaves + iW)*nAge, age, nAge, 
                                            age[0], ageStep[iA]);
                else
                    pData[iW] = trapz_table(data + (iZ*nWaves + iW)*nAge, age, nAge, 
                                            ageStep[iA - 1], ageStep[iA]);
            }
        }
    memcpy(spectra->ready, spectra->integrated, nZ*nAgeStep*nWaves*sizeof(double));
    #ifdef TIMING
        profiler_end(INTEGRATION);
    #endif
}
 

struct dust_params {
    double tauUV_ISM;
    double nISM;
    double tauUV_BC;
    double nBC;
    double tBC;
};


inline double *dust_absorption(struct sed_params *spectra, struct dust_params *dustParams) {
    /* tBC: life time of the birth clound
     * nu: fraction of ISM dust absorption
     * tauUV: V-band absorption optical depth
     * nBC: power law index of tauBC
     * nISM: power law index of tauISM
     * 
     * Reference: da Cunha et al. 2008
     */
    #ifdef TIMING
        profiler_start("Dust absorption", DUST);
    #endif
    int nAge = spectra->nAge;
    double *age = spectra->age;
    double *rawData = spectra->raw;
    int nWaves = spectra->nWaves;
    double *waves = spectra->waves;
    int nZ = spectra->nZ;

    int nAgeStep = spectra->nAgeStep;
    double *ageStep = spectra->ageStep;
    memcpy(spectra->ready, spectra->integrated, nZ*nAgeStep*nWaves*sizeof(double));
    double *data = spectra->ready;

    int iAgeBC;
    double t0, t1;

    double tauUV_ISM = dustParams->tauUV_ISM;
    double nISM = dustParams->nISM;
    double tauUV_BC = dustParams->tauUV_BC;
    double nBC = dustParams->nBC;
    double tBC = dustParams->tBC;

    double *transISM = malloc(nWaves*sizeof(double));
    double *transBC = malloc(nWaves*sizeof(double));

    // Find the time inverval containning the birth cloud
    if (tBC >= ageStep[nAgeStep - 1]) {
        iAgeBC = nAgeStep;
        t0 = 0.;
        t1 = 0.;
    }
    else if(tBC < ageStep[0]) {
        iAgeBC = 0;
        t0 = age[0];
        t1 = ageStep[0];
    }
    else {
        iAgeBC = bisection_search(tBC, ageStep, nAgeStep) + 1;
        t0 = ageStep[iAgeBC - 1];
        t1 = ageStep[iAgeBC];
    } 
    
    // Compute the optical depth of both the birth cloud and the ISM
    #pragma omp parallel \
    default(none) \
    firstprivate(nAge, age, rawData, \
                 nWaves, waves, nZ, nAgeStep, ageStep, data, \
                 iAgeBC, t0, t1, \
                 tauUV_ISM, nISM, tauUV_BC, nBC, tBC, \
                 transISM, transBC) \
    num_threads(g_nThread)
    {
        int iW, i, n;
        double *pData; 
        double ratio;        
        
        #pragma omp for schedule(static,1)
        for(iW = 0; iW < nWaves; ++iW) {
            ratio = waves[iW]/1600.;
            transISM[iW] = exp(-tauUV_ISM*pow(ratio, nISM));
            transBC[iW] = exp(-tauUV_BC*pow(ratio, nBC));
        }
        
        // t_s < tBC < t_s + dt
        if (iAgeBC != nAgeStep) {
            n = nZ*nWaves;
            #pragma omp for schedule(static,1) 
            for(i = 0; i < n; ++i) {
                iW = i%nWaves;
                pData = data + (i/nWaves*nAgeStep + iAgeBC)*nWaves;
                pData[iW] = transBC[iW] \
                    *trapz_table(rawData + i*nAge, age, nAge, t0, tBC) \
                    + trapz_table(rawData + i*nAge, age, nAge, tBC, t1);
            }     
        }
        
        // tBC > t_s       
        n = iAgeBC*nZ;
        #pragma omp for schedule(static,1) 
        for(i = 0; i < n; ++i) {
            pData = data + (i/iAgeBC*nAgeStep + i%iAgeBC)*nWaves;
            for(iW = 0; iW < nWaves; ++iW) 
                pData[iW] *= transBC[iW];
        }
        
        n = nAgeStep*nZ;
        #pragma omp for schedule(static,1) 
        for(i = 0; i < n; ++i) {
            pData = data + i*nWaves;
            for(iW = 0; iW < nWaves; ++iW) 
                pData[iW] *= transISM[iW];
        }
    }

    free(transISM);
    free(transBC);
    #ifdef TIMING
        profiler_end(DUST);
    #endif   
    return data;
}


inline void templates_working(struct sed_params *spectra, double z) {
    int nWaves = spectra->nWaves;
    double *waves = spectra->waves;
    int nZ = spectra->nZ;

    int nFlux = spectra->nFlux;
    int nObs = spectra->nObs;
    int *nFilterWaves = spectra->nFilterWaves;
    double *filterWaves = spectra->filterWaves;
    double *filters = spectra->filters;
    double *LyAbsorption = spectra->LyAbsorption;

    int nAge = spectra->nAgeStep;
    double *readyData = spectra->ready;
    double *workingData = spectra->working;

    double *obsWaves = NULL;
    double *obsData = NULL;
    if (nObs > 0) {
        obsWaves = (double*)malloc(nWaves*sizeof(double));
        obsData = (double*)malloc(nZ*nAge*nWaves*sizeof(double));
    }
    // Spectra to be interploated along metallicities
    // The first dimension refers to filters/wavelengths and ages
    // Thw last dimension refers to metallicites
    double *refSpectra = malloc(nFlux*nAge*nZ*sizeof(double));
    
    int minZ = spectra->minZ;
    int maxZ = spectra->maxZ;   
    double *Z = spectra->Z;

    #pragma omp parallel \
    default(none)  \
    firstprivate(z, nWaves, waves, nZ, \
                 nFlux, nObs, nFilterWaves, filterWaves, filters, LyAbsorption, \
                 nAge, readyData, workingData, obsWaves, obsData, refSpectra, \
                 minZ, maxZ, Z) \
    num_threads(g_nThread) 
    {
        int iA, iW, iZ, iAZ, iF, iFW, i, n;
        int nRest = nFlux - nObs;
        double *pData;
        double *pObsData;
        double *pRefData;
        int nFW;
        double *pFilterWaves = filterWaves;
        double *pFilters = filters;
        double *filterData;
        double I;
        double interpZ;
        
        #pragma omp single
        if (nObs > 0) {
            // Transform everything to observer frame
            // Note the fluxes in this case is a function of wavelength
            // Therefore the fluxes has a factor of 1/(1 + z)
            for(iW = 0; iW < nWaves; ++iW)
                obsWaves[iW] = waves[iW]*(1. + z);
            for(iAZ = 0; iAZ < nAge*nZ; ++iAZ) {
                pData = readyData + iAZ*nWaves;
                pObsData = obsData + iAZ*nWaves;
                for(iW = 0; iW < nWaves; ++iW)
                    pObsData[iW] = pData[iW]/(1. + z);           
            }
            if (LyAbsorption != NULL)
                // Add IGM absorption
                 for(iAZ = 0; iAZ < nAge*nZ; ++iAZ) {
                    pObsData = obsData + iAZ*nWaves;
                    for(iW = 0; iW < nWaves; ++iW)
                        pObsData[iW] *= LyAbsorption[iW];
                    }       
        }
        if (filters == NULL) {
            // Tranpose the templates such that the last dimension is the metallicity
            #pragma omp single
            if (nObs > 0) {
                for(iZ = 0; iZ < nZ; ++iZ) 
                    for(iA = 0; iA < nAge; ++iA) {
                        // Use fluxes in the observer frame
                        pObsData = obsData + (iZ*nAge + iA)*nWaves;
                        for(iW = 0; iW < nWaves; ++iW)
                            refSpectra[(iW*nAge + iA)*nZ + iZ] = pObsData[iW];
                    }
            }
            else {
                for(iZ = 0; iZ < nZ; ++iZ) 
                    for(iA = 0; iA < nAge; ++iA) {
                        pData = readyData + (iZ*nAge + iA)*nWaves;
                        for(iW = 0; iW < nWaves; ++iW)
                            refSpectra[(iW*nAge + iA)*nZ + iZ] = pData[iW];
                    }
            }
        }
        else {
            // Intgrate SED templates over filters
            #ifdef TIMING
                #pragma omp single
                profiler_start("Filter fluxes", WORKING1);
            #endif
            pRefData = refSpectra;
            for(iF = 0; iF < nRest; ++iF) {
                nFW = nFilterWaves[iF];
                filterData = (double*)malloc(nFW*sizeof(double));
                n = nAge*nZ;
                #pragma omp for schedule(static,1) 
                for(i = 0; i < n; ++i) {
                    pData = readyData + (i%nZ*nAge + i/nZ)*nWaves;
                    for(iFW= 0; iFW < nFW; ++iFW)
                        filterData[iFW] = interp(pFilterWaves[iFW], 
                                                 waves, pData, nWaves);
                    for(iFW = 0; iFW < nFW; ++iFW)
                        filterData[iFW] *= pFilters[iFW];
                    I = 0.;
                    for(iFW = 1; iFW < nFW; ++iFW)
                        I += (pFilterWaves[iFW] - pFilterWaves[iFW - 1]) \
                             *(filterData[iFW] + filterData[iFW - 1]);
                    pRefData[iF*n + i] = I/2.;
                }
                free(filterData);
                pFilterWaves += nFW;
                pFilters += nFW;
            }
            for(iF = nRest; iF < nFlux; ++iF) {
                nFW = nFilterWaves[iF];
                filterData = (double*)malloc(nFW*sizeof(double));
                n = nAge*nZ;
                #pragma omp for schedule(static,1) 
                for(i= 0; i< nAge*nZ; ++i) {
                    pObsData = obsData + (i%nZ*nAge + i/nZ)*nWaves;
                    for(iFW = 0; iFW < nFW; ++iFW)
                        filterData[iFW] = interp(pFilterWaves[iFW], 
                                                 obsWaves, pObsData, nWaves);
                    for(iFW = 0; iFW < nFW; ++iFW)
                        filterData[iFW] *= pFilters[iFW];
                    I = 0.;
                    for(iFW = 1; iFW < nFW; ++iFW)
                        I += (pFilterWaves[iFW] - pFilterWaves[iFW - 1]) \
                             *(filterData[iFW] + filterData[iFW - 1]);
                    pRefData[iF*n + i] = I/2.;
                }
                free(filterData);
                pFilterWaves += nFW;
                pFilters += nFW;
            }
            #ifdef TIMING
                #pragma omp single
                profiler_end(WORKING1);
            #endif
        }

        // Interploate SED templates along metallicities
        #ifdef TIMING
            #pragma omp single
            profiler_start("Interploation over metallicites", WORKING2);
        #endif
        n = (maxZ - minZ + 1)*nAge;
        #pragma omp for schedule(static,1)
        for(i = 0; i < n; ++i) {
            interpZ = (minZ + i/nAge + 1.)/1000.;
            pData = workingData + i*nFlux;
            for(iF = 0; iF < nFlux; ++iF) 
                pData[iF] = interp(interpZ, Z, refSpectra + (iF*nAge+ i%nAge)*nZ, nZ);
        }
        #ifdef TIMING
            #pragma omp single
            profiler_end(WORKING2);
        #endif
    }
    
    if (nObs > 0) {
        free(obsWaves);
        free(obsData);
    }
    free(refSpectra);

}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Primary Functions                                                           *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
double *composite_spectra_cext(struct sed_params *spectra,
                               struct gal_params *galParams, struct dust_params *dustParams,
                               short outType, short nThread) {
    #ifdef TIMING
        init_profiler();
        timing_start("Compute magnitudes");
    #endif
    g_nThread = nThread;

    int iF, iG, iP, iFG;
    // Initialise galaxies parameters
    double z = galParams->z;
    int nAgeStep= galParams->nAgeStep;
    double *ageStep = galParams->ageStep;
    int nGal = galParams->nGal;
    struct csp *pHistories = galParams->histories;
    struct ssp *pBursts;
    int nProg;
    double sfr;
    int metals;

    // Generate templates
    int nFlux = spectra->nFlux;
    double *logWaves = spectra->logWaves;
    init_templates_working(spectra, ageStep, nAgeStep, nFlux);
    integrate_templates_raw(spectra);
    double *workingData = spectra->working;
    double *pWorkingData;

    // Initialise outputs
    double *output = malloc(nGal*nFlux*sizeof(double));
    double *pOutput = output;
    for(iFG = 0; iFG < nGal*nFlux; ++iFG)
        *pOutput++ = TOL;
    pOutput = output;

    int minZ = spectra->minZ;
    int maxZ = spectra->maxZ;

    if (dustParams == NULL)
        templates_working(spectra, z);
    for(iG = 0; iG < nGal; ++iG) {
        // Add dust absorption to SED templates
        if (dustParams != NULL) {
            dust_absorption(spectra, dustParams + iG);
            templates_working(spectra, z);
        }
        // Sum contributions from all progenitors
        nProg = pHistories->nBurst;
        for(iP = 0; iP < nProg; ++iP) {
            pBursts = pHistories->bursts + iP;
            sfr = pBursts->sfr;
            metals = (int)(pBursts->metals*1000 - .5);
            if (metals < minZ)
                metals = minZ;
            else if (metals > maxZ)
                metals = maxZ;
            pWorkingData = workingData + (metals*nAgeStep + pBursts->index)*nFlux;
            for(iF = 0 ; iF < nFlux; ++iF)
                pOutput[iF] += sfr*pWorkingData[iF];
        }
        ++pHistories;
        pOutput += nFlux;

        #ifdef TIMING
            report(iG, nGal);
        #endif
    }
    free_templates_working(spectra);

    if (outType == 0) {
        pOutput = output;
        for(iFG = 0; iFG < nFlux*nGal; ++iFG) {
            *pOutput = M_AB(*pOutput);
            ++pOutput;
        }
        #ifdef TIMING
            timing_end();
            profiler_summary();
        #endif
        return output;
    }
    else if (outType == 1) {
        #ifdef TIMING
            timing_end();
            profiler_summary();
        #endif
        return output;
    }
    
    // Fit UV slopes
    int nR = 3;
    struct linResult result;

    output = (double*)realloc(output, (nFlux + nR)*nGal*sizeof(double));
    pOutput = output + nFlux*nGal;
    double *pFit = output;

    int nFit = nFlux - 1;
    double *logf = malloc(nFit*sizeof(double));

    for(iG = 0; iG < nGal; ++iG) {
        for(iF = 0; iF < nFit; ++iF) 
            logf[iF] = log(pFit[iF]);
        pFit += nFlux;
        //printf("waves = %.1f, logf = %.1f\n", logWaves[1], logf[1]);
        result = linregress(logWaves, logf, nFit);
        pOutput[0] = (double)result.slope;
        pOutput[1] = (double)result.intercept;
        pOutput[2] = (double)result.R;
        pOutput += nR;
        //printf("Slope = %.1f\n", result.slope);
    }
    // Convert to AB magnitude
    pOutput = output + nFit;
    for(iG = 0; iG < nGal; ++iG) {
        *pOutput = M_AB(*pOutput);
        pOutput += nFlux;
    }

    #ifdef TIMING
        timing_end();
        profiler_summary();
    #endif
    return output;
}

