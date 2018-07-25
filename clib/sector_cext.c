#define _SECTOR_
#include"tools.h"
#include"sector_cext.h"
#include"hdf5.h"
#include"hdf5_hl.h"

//#define SURFACE_AREA 1.1965e40 // 4*pi*(10 pc)**2 unit cm^2
//#define JANSKY(x) (3.34e4*(x)*(x))
#define M_AB(x) (-2.5*log10(x) + 8.9) // Convert Jansky to AB magnitude
#define TOL 1e-30 // Minimum Flux

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
    #define SUM 4
    #define FIT 5

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
            printf("#  call: %6d  total: %2.3f sec  mean: %2.3f ms\n",
                   ncall, timer[iB], timer[iB]/ncall*1e3);
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
 * Struct to store galaxy properites                                           *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void trim_gal_params(struct gal_params *galParams, int minZ, int maxZ) {
    /* Set the metallicity of each SSP to the given range */
    int iB, iG;

    float f_minZ = (minZ + 1.)/1000.;
    float f_maxZ = (maxZ + 1.)/1000.;
    int metals;
    int nGal = galParams->nGal;
    struct csp *pHistories = galParams->histories;
    int nBurst;
    struct ssp *pBursts;

    for(iG = 0; iG < nGal; ++iG) {
        nBurst = pHistories->nBurst;
        pBursts = pHistories->bursts;
        for(iB = 0; iB < nBurst; ++iB) {
            metals = (int)(pBursts->metals*1000 - .5);
            if (metals < minZ)
                pBursts->metals = f_minZ;
            else if (metals > maxZ)
                pBursts->metals = f_maxZ;
            ++pBursts;
        }
        ++pHistories;
    }
}


int *age_flag(struct csp *histories, int nAgeStep) {
    int iA, iB;

    int *ageFlag = malloc(nAgeStep*sizeof(int));
    int nB = histories->nBurst;
    struct ssp *bursts = histories->bursts;

    for(iA = 0; iA < nAgeStep; ++iA)
        ageFlag[iA] = 1;
    for(iB = 0; iB < nB; ++iB)
        ageFlag[bursts[iB].index] = 0;

    return ageFlag;
}


int *Z_flag(struct csp *histories, int nMaxZ) {
    int iZ, iB;

    int *ZFlag = malloc(nMaxZ*sizeof(int));
    int nB = histories->nBurst;
    struct ssp *bursts = histories->bursts;

    for(iZ = 0; iZ < nMaxZ; ++iZ)
        ZFlag[iZ] = 1;
    for(iB = 0; iB < nB; ++iB)
        ZFlag[(int)(1000*bursts[iB].metals - 0.5)] = 0;

    return ZFlag;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Functions to process SEDs                                                   *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void init_templates_raw(struct sed_params *spectra, char *fName) {
    /* File name should be "sed_library.hdf5" and contain:
     * "/metals" Metallicity grid (1-D dataset)
     * "/waves" Wavelength grid [AA] (1-D dataset)
     * "/age" Stellar age grid [yr] (1-D dataset)
     * "/flux" Flux density [erg/s/AA/cm^2] (3-D dataset)
     *
     * Flux densities should be normlised by the surface area of 10 pc
     * sphere. The first, second, and third dimensions should be metallicty,
     * wavelength and stellar age respectively.
     */
    hid_t file_id;
    hsize_t dims[3];

    printf("#***********************************************************\n");
    printf("# Read SED templates\n");
    file_id = H5Fopen(fName, H5F_ACC_RDONLY, H5P_DEFAULT);
    // Read dimensions
    H5LTget_dataset_info(file_id, "/flux", dims, NULL, NULL);
    int nZ = dims[0];
    int nWaves = dims[1];
    int nAge = dims[2];
    // Read metallicity range
    spectra->nZ = nZ;
    spectra->Z = (double*)malloc(nZ*sizeof(double));
    H5LTread_dataset_double(file_id, "/metals", spectra->Z);
    spectra->minZ = (int)(spectra->Z[0]*1000 - .5);
    spectra->maxZ = (int)(spectra->Z[nZ - 1]*1000 - .5);
    spectra->nMaxZ = spectra->maxZ - spectra->minZ + 1;
    printf("# Metallicity range:\n#\t%.3f to %.3f\n", spectra->Z[0], spectra->Z[nZ - 1]);
    // Read wavelength
    spectra->nWaves = nWaves;
    spectra->waves = (double*)malloc(nWaves*sizeof(double));
    H5LTread_dataset_double(file_id, "/waves", spectra->waves);
    printf("# Wavelength range:\n#\t%.1f AA to %.1f AA\n",
           spectra->waves[0], spectra->waves[nWaves - 1]);
    // Read stellar age
    spectra->nAge = nAge;
    spectra->age = (double*)malloc(nAge*sizeof(double));
    H5LTread_dataset_double(file_id, "/age", spectra->age);
    printf("# Stellar age range:\n#\t%.2f Myr to %.2f Myr\n",
           spectra->age[0]*1e-6, spectra->age[nAge - 1]*1e-6);
    // Read flux
    spectra->raw = (double*)malloc(nZ*nWaves*nAge*sizeof(double));
    H5LTread_dataset_double(file_id, "/flux", spectra->raw);
    H5Fclose(file_id);
    printf("#***********************************************************\n\n");
}


void shrink_templates_raw(struct sed_params *spectra, double maxAge) {
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
    int nFW;
    int *nFilterWaves = spectra->nFilterWaves;
    double *filterWaves = spectra->filterWaves;
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

    // Find shrinked wavelength ranges
    printf("#***********************************************************\n");
    printf("# Shrinked wavelength ranges:\n");
    for(iW = 1; iW < nWaves; ++iW) {
        w = waves[iW];
        pFilterWaves = filterWaves;
        outFlag = 1;
        for(iF = 0; iF < nFlux; ++iF) {
            nFW = nFilterWaves[iF];
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

    free(wavesIndices);
}


void init_filters(struct sed_params *spectra,
                  double *betaBands, int nBeta, double *restBands, int nRest,
                  double *obsTrans, double *obsWaves, int *nObsWaves, int nObs, double z) {
    // Initialise filter parameters
    //   -This function should be called after ``init_templates_raw``
    int iF, iW;
    double *pData = betaBands;

    int nFlux = nBeta + nRest + nObs;
    int nWaves = spectra->nWaves;
    double *waves = spectra->waves;
    int nNewWaves;
    double *newWaves = malloc(nWaves*sizeof(double));
    double w;
    double lower, upper;
    double norm;

    int *nFilterWaves = malloc(nFlux*sizeof(int));
    double *filterWaves = NULL;
    double *filters = NULL;
    int offset = 0;

    double *centreWaves = malloc(nFlux*sizeof(double));
    double *logWaves = malloc(nFlux*sizeof(double));

    // Initialise rest-frame filters and those that are related to the UV slope
    nRest += nBeta;
    for(iF = 0; iF < nRest; ++iF) {
        if (iF == nBeta)
            pData = restBands;
        //   -Construct wavelengths for a new filter
        lower = *pData++;
        upper = *pData++;
        if (lower < waves[0] || upper > waves[nWaves - 1]) {
            printf("Filter wavelength range are beyond SED templates\n");
            exit(0);
        }
        nNewWaves = 1;
        newWaves[0] = lower;
        for(iW = 0; iW < nWaves; ++iW) {
            w = waves[iW];
            if (w > lower && w < upper)
                newWaves[nNewWaves++] = w;
        }
        newWaves[nNewWaves++] = upper;
        nFilterWaves[iF] = nNewWaves;
        filterWaves = realloc(filterWaves, (offset + nNewWaves)*sizeof(double));
        memcpy(filterWaves + offset, newWaves, nNewWaves*sizeof(double));
        //  -Construct transmissions
        filters = realloc(filters, (offset + nNewWaves)*sizeof(double));
        if (iF < nBeta) {
            norm = 1./(upper - lower);
            for(iW = 0; iW < nNewWaves; ++iW)
                *(filters + offset + iW) = norm;
        }
        else {
            norm = 1./log(upper/lower);
            for(iW = 0; iW < nNewWaves; ++iW)
                *(filters + offset + iW) = 3.34e4*newWaves[iW]*norm;
        }
        offset += nNewWaves;
        //
        w = (lower + upper)/2.;
        centreWaves[iF] = w;
        logWaves[iF] = log(w);
    }

    // Initialise observer-frame filters
    if (nObs > 0) {
        memcpy(nFilterWaves + nRest, nObsWaves, nObs*sizeof(int));
        nNewWaves = 0;
        for(iF = 0; iF < nObs; ++iF)
            nNewWaves += nObsWaves[iF];
        filterWaves = realloc(filterWaves, (offset + nNewWaves)*sizeof(double));
        memcpy(filterWaves + offset, obsWaves, nNewWaves*sizeof(double));
        filters = realloc(filters, (offset + nNewWaves)*sizeof(double));
        memcpy(filters + offset, obsTrans, nNewWaves*sizeof(double));
        for(iF = 0; iF < nObs; ++iF) {
            nNewWaves = nObsWaves[iF];
            lower = filterWaves[offset];
            upper = filterWaves[offset + nNewWaves - 1];
            if (lower < waves[0]*(1. + z) || upper > waves[nWaves - 1]*(1. + z)) {
                printf("Filter wavelength range are beyond SED templates\n");
                exit(0);
            }
            //   -Divide the transmission by the wavelength to calculate the normalisation
            for(iW = 0; iW < nNewWaves; ++iW)
                filters[offset + iW] /= filterWaves[offset + iW];
            norm = 1./trapz_table(filters + offset, filterWaves + offset, nNewWaves, lower, upper);
            //   -Compute the pivot wavelength and transform it to the rest-frame
            for(iW = 0; iW < nNewWaves; ++iW)
                filters[offset + iW] *= filterWaves[offset + iW]*filterWaves[offset + iW];
            w = sqrt(norm*trapz_table(filters + offset, filterWaves + offset,
                                      nNewWaves, lower, upper))/(1. + z);
            centreWaves[nRest + iF] = w;
            logWaves[nRest + iF] = log(w);
            //  -Transform the filter such that it can give AB magnitude
            for(iW = 0; iW < nNewWaves; ++iW)
                filters[offset + iW] *= 3.34e4*norm;
            //  -Transform the filter to the rest-frame
            for(iW = 0; iW < nNewWaves; ++iW)
                filterWaves[offset + iW] /= 1. + z;
            offset += nNewWaves;
        }
    }
    else
        //   -Disable IGM absorption if there is no observer-frame filter
        spectra->igm = 0;

    spectra->nFlux = nFlux;
    spectra->nObs = nObs;
    spectra->nFilterWaves = nFilterWaves;
    spectra->filterWaves = filterWaves;
    spectra->filters = filters;
    spectra->centreWaves = centreWaves;
    spectra->logWaves = logWaves;

    free(newWaves);
}



void init_templates_integrated(struct sed_params *spectra) {
    #ifdef TIMING
        profiler_start("Integration over time", INTEGRATION);
    #endif
    int iA, iW, iZ;
    double *pData;

    int nAgeStep = spectra->nAgeStep;
    double *ageStep = spectra->ageStep;
    int nAge = spectra->nAge;
    double *age = spectra->age;
    int nWaves = spectra->nWaves;
    int nZ = spectra->nZ;
    double *data = spectra->raw;
    // Spectra after integration over time
    // The first dimension refers to metallicites and ages
    // The last dimension refers to wavelengths
    double *intData = malloc(nZ*nAgeStep*nWaves*sizeof(double));

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
    spectra->integrated = intData;
    #ifdef TIMING
        profiler_end(INTEGRATION);
    #endif
}


void init_templates_working(struct sed_params *spectra, struct csp *pHistories,
                                   struct dust_params *dustParams, int iG) {
    int *ageFlag = NULL;
    int *ZFlag = NULL;

    int minZ = spectra->minZ;
    int nMaxZ = spectra->nMaxZ;
    int nZ = spectra->nZ;
    int nWaves = spectra->nWaves;
    int nAgeStep = spectra->nAgeStep;
    int nFlux = spectra->nFlux;
    double *readyData = spectra->ready;

    if (dustParams != NULL) {
        memcpy(readyData, spectra->integrated, nZ*nAgeStep*nWaves*sizeof(double));
        ageFlag = age_flag(pHistories, nAgeStep);
        ZFlag = Z_flag(pHistories, nMaxZ);
        dust_absorption(spectra, dustParams + iG);
    }
    else if (iG == -1) {
        memcpy(readyData, spectra->integrated, nZ*nAgeStep*nWaves*sizeof(double));
        ageFlag = calloc(nAgeStep, sizeof(int));
        ZFlag = calloc(nMaxZ, sizeof(int));
    }
    else
        return;

    int iA, iF, iW, iZ, i, n;
    double *pData;
    // Spectra to be interploated along metallicities
    // Dimension: filter/wavelength × age × metallicity
    double *refSpectra = NULL;

    if (spectra->filters == NULL) {
        double z = spectra->z;

        n = nZ*nAgeStep;
        refSpectra = malloc(nWaves*n*sizeof(double));
        if (spectra->nObs > 0) {
            // Transform everything to observer frame
            // Note the fluxes in this case is a function of wavelength
            // Therefore the fluxes has a factor of 1/(1 + z)
            for(i = 0; i < n; ++i) {
                pData = readyData + i*nWaves;
                for(iW = 0; iW < nWaves; ++iW)
                   pData[iW] /= 1. + z;
            }
        }
        // Tranpose the templates such that the last dimension is the metallicity
        for(i = 0; i < n; ++i) {
            pData = readyData + i*nWaves;
            for(iW = 0; iW < nWaves; ++iW)
                refSpectra[(iW*nAgeStep + i%nAgeStep)*nZ + i/nAgeStep] = pData[iW];
        }

    }
    else {
        // Intgrate SED templates over filters
        int nFW;
        int *nFilterWaves = spectra->nFilterWaves;
        double *pFilterWaves = spectra->filterWaves;
        double *pFilters = spectra->filters;
        double *waves = spectra->waves;

        n = nZ*nAgeStep;
        refSpectra = malloc(nFlux*n*sizeof(double));
        for(iF = 0; iF < nFlux; ++iF) {
            nFW = nFilterWaves[iF];
            for(i = 0; i < n; ++i) {
                iA = i/nZ;
                if (ageFlag[iA])
                    continue;
                refSpectra[iF*n + i] = trapz_filter(pFilters, pFilterWaves, nFW,
                                                    readyData + (i%nZ*nAgeStep + iA)*nWaves,
                                                    waves, nWaves);
            }
            pFilterWaves += nFW;
            pFilters += nFW;
        }
    }
    // Interploate SED templates along metallicities
    double *Z = spectra->Z;
    double *workingData = spectra->working;
    double interpZ;

    n = nMaxZ*nAgeStep;
    for(i = 0; i < n; ++i) {
        iZ = i/nAgeStep;
        if (ZFlag[iZ])
            continue;
        iA = i%nAgeStep;
        if (ageFlag[iA])
            continue;
        interpZ = (minZ + iZ + 1.)/1000.;
        pData = workingData + i*nFlux;
        for(iF = 0; iF < nFlux; ++iF)
            pData[iF] = interp(interpZ, Z, refSpectra + (iF*nAgeStep+ iA)*nZ, nZ);
    }

    free(ageFlag);
    free(ZFlag);
    free(refSpectra);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Primary Functions                                                           *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void fit_UV_slope(double *pTarget, double *pFit, int nGal, int nFlux,
                  double *logWaves, int nFit, int nR) {
    #ifdef TIMING
        profiler_start("Slope fit", FIT);
    #endif
    int iF, iG;
    struct linResult result;
    double *logf = malloc(nFit*sizeof(double));

    for(iG = 0; iG < nGal; ++iG) {
        for(iF = 0; iF < nFit; ++iF)
            logf[iF] = log(pFit[iF]);
        pFit += nFlux;
        result = linregress(logWaves, logf, nFit);
        *pTarget++ = (double)result.slope;
        if (nR > 1)
            *pTarget++ = (double)result.intercept;
        if (nR > 2)
            *pTarget++ = (double)result.R;
    }
    #ifdef TIMING
        profiler_end(FIT);
    #endif
}


void compute_spectra(double *target, struct sed_params *spectra,
                     struct gal_params *galParams, struct dust_params *dustParams,
                     int approx, short nThread) {
    // Initialise SED templates
    spectra->ready = NULL;
    spectra->working = NULL;

    // Integrate SED templates over given time steps
    spectra->nAgeStep = galParams->nAgeStep;
    spectra->ageStep = galParams->ageStep;
    init_templates_integrated(spectra);
    init_IGM_absorption(spectra);

    // Initialise templates for birth cloud if necessary
    if (dustParams == NULL) {
        spectra->inBC = NULL;
        spectra->outBC = NULL;
    }

    #ifdef TIMING
        profiler_start("Summation over progenitors", SUM);
    #endif
    #pragma omp parallel \
    default(none) \
    firstprivate(spectra, galParams, dustParams, target, approx) \
    num_threads(nThread)
    {
        int iF, iG;
        int offset;
        double *pTarget;

        int iP, nProg;
        int nAgeStep = galParams->nAgeStep;
        int nGal = galParams->nGal;
        struct csp *histories = galParams->histories;
        struct csp *pHistories;
        struct ssp *pBursts;
        double sfr;
        int metals;

        int nFlux = spectra->nFlux;
        struct sed_params omp_spectra;
        memcpy(&omp_spectra, spectra, sizeof(struct sed_params));
        double *readyData = malloc(spectra->nZ*nAgeStep*spectra->nWaves*sizeof(double));
        double *workingData = malloc(spectra->nMaxZ*nAgeStep*nFlux*sizeof(double));
        omp_spectra.ready = readyData;
        omp_spectra.working = workingData;


        if (approx) {
            init_templates_working(&omp_spectra, histories, NULL, -1);
            init_templates_special(&omp_spectra, dustParams->tBC, approx);

            int iAge;
            // Assume that tBC is the same for every galaxy
            int iAgeBC = birth_cloud_interval(dustParams->tBC, galParams->ageStep, nAgeStep);
            double *inBC = omp_spectra.inBC;
            double *outBC = omp_spectra.outBC;
            double *inBCFlux = malloc(nFlux*sizeof(double));
            double *outBCFlux = malloc(nFlux*sizeof(double));

            #pragma omp for schedule(static, 1)
            for(iG = 0; iG < nGal; ++iG) {
                pHistories = histories + iG;
                // Sum contributions from all progenitors
                nProg = pHistories->nBurst;
                pTarget = target + iG*nFlux;
                for(iF = 0; iF < nFlux; ++iF) {
                    inBCFlux[iF] = 0.;
                    outBCFlux[iF] = 0.;
                }
                for(iP = 0; iP < nProg; ++iP) {
                    pBursts = pHistories->bursts + iP;
                    iAge = pBursts->index;
                    sfr = pBursts->sfr;
                    metals = (int)(pBursts->metals*1000 - .5);
                    if (iAge > iAgeBC) {
                        offset = (metals*nAgeStep + iAge)*nFlux;
                        for(iF = 0; iF < nFlux; ++iF)
                            outBCFlux[iF] += sfr*workingData[offset + iF];
                    }
                    else if (iAge == iAgeBC) {
                        offset = metals*nFlux;
                        for(iF = 0; iF < nFlux; ++iF) {
                            inBCFlux[iF] += sfr*inBC[offset + iF];
                            outBCFlux[iF] += sfr*outBC[offset + iF];
                        }
                    }
                    else {
                        offset = (metals*nAgeStep + iAge)*nFlux;
                        for(iF = 0; iF < nFlux; ++iF)
                            inBCFlux[iF] += sfr*workingData[offset +iF];
                    }
                }

                // Apply dust absorption
                dust_absorption_approx(inBCFlux, outBCFlux, omp_spectra.centreWaves, nFlux,
                                       dustParams + iG);

                for(iF = 0; iF < nFlux; ++iF)
                    pTarget[iF] += inBCFlux[iF] + outBCFlux[iF];
                #ifdef TIMING
                    report(iG, nGal);
                #endif
            }
            free(inBC);
            free(outBC);
            free(inBCFlux);
            free(outBCFlux);
        }
        else {
            if (dustParams == NULL)
                init_templates_working(&omp_spectra, histories, NULL, -1);
            else
                //  -Assume that tBC is the same for every galaxy
                init_templates_special(&omp_spectra, dustParams->tBC, approx);

            #pragma omp for schedule(static, 1)
            for(iG = 0; iG < nGal; ++iG) {
                pHistories = histories + iG;
                // Add dust absorption to SED templates
                init_templates_working(&omp_spectra, pHistories, dustParams, iG);
                // Sum contributions from all progenitors
                nProg = pHistories->nBurst;
                pTarget = target + iG*nFlux;
                for(iP = 0; iP < nProg; ++iP) {
                    pBursts = pHistories->bursts + iP;
                    sfr = pBursts->sfr;
                    metals = (int)(pBursts->metals*1000 - .5);
                    offset = (metals*nAgeStep + pBursts->index)*nFlux;
                    for(iF = 0 ; iF < nFlux; ++iF)
                        pTarget[iF] += sfr*workingData[offset + iF];
                }
                #ifdef TIMING
                    report(iG, nGal);
                #endif
            }
        }
        free(readyData);
        free(workingData);
    }
    free(spectra->integrated);
    #ifdef TIMING
        profiler_end(SUM);
    #endif
}


double *composite_spectra_cext(struct sed_params *spectra,
                               struct gal_params *galParams, struct dust_params *dustParams,
                               short outType, short approx, short nThread) {
    // Trim the metallicity of each SSP such that it is within the range of
    // input SED templates
    trim_gal_params(galParams, spectra->minZ, spectra->maxZ);

    #ifdef TIMING
        init_profiler();
        timing_start("Compute magnitudes");
    #endif

    // Initialise outputs
    int iGF;
    int nGal = galParams->nGal;
    int nFlux = spectra->nFlux;
    int nGF = nGal*nFlux;
    double *output = malloc(nGF*sizeof(double));
    double *pOutput = output;
    for(iGF = 0; iGF < nGF; ++iGF)
        *pOutput++ = TOL;

    // Compute spectra for given inputs
    //   -Only use approximation when considering dust
    compute_spectra(output, spectra, galParams, dustParams,
                    (short)(approx && dustParams != NULL), nThread);

    if (outType == 0) {
        // Convert to AB magnitude
        pOutput = output;
        for(iGF = 0; iGF < nGF; ++iGF) {
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
    int iG;
    int nR = 3;
    output = (double*)realloc(output, (nFlux + nR)*nGal*sizeof(double));
    fit_UV_slope(output + nFlux*nGal, output, nGal, nFlux, spectra->logWaves, nFlux - 1, nR);

    // Convert to AB magnitude
    pOutput = output + nFlux - 1;
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