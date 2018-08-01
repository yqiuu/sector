#define _SECTOR_
#include"tools.h"
#include"sector.h"

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
 * Galaxy formation model interface                                            *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifdef CALC_MAGS
void init_templates_mini(mini_sed_params_t *miniSpectra, char *fName,
                         double *LTTime, int *targetSnap, double *redshifts,
                         double *restBands, int nRest, int nBeta, double tBC) {
    // Initialise full templates
    int iA, iS;
    struct sed_params spectra[MAGS_N_SNAPS];
    int nAgeStep;
    double *ageStep;

    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        nAgeStep = targetSnap[iS];
        //// Initialise raw templates
        init_templates_raw(spectra + iS, fName);
        //// Initialise filters
        init_filters(spectra + iS, NULL, 0, restBands, nRest,
                     NULL, NULL, NULL, 0, 1. + redshifts[iS]);
        if (spectra[iS].nFlux != MAGS_N_BANDS) {
            printf("MAGS_N_BANDS does not match!\n");
            exit(EXIT_FAILURE);
        }
        //// Initialise time step
        spectra[iS].nAgeStep = nAgeStep;
        ageStep = (double*)malloc(nAgeStep*sizeof(double));
        ////   -Should be in a unit of yr
        for(int iA = 0; iA < nAgeStep; ++iA)
            ageStep[iA] = LTTime[nAgeStep - iA - 1] - LTTime[nAgeStep];
        spectra[iS].ageStep = ageStep;
        ////   -This function may be omitted
        shrink_templates_raw(spectra + iS, ageStep[nAgeStep - 1]);
        ////   -Disable IGM absorption
        spectra[iS].igm = 0;
        //// Integrate templates over given time steps
        init_templates_integrated(spectra + iS);
        //// Initialise working templates
        spectra[iS].ready = \
        (double*)malloc(spectra[iS].nZ*nAgeStep*spectra[iS].nWaves*sizeof(double));
        spectra[iS].working = \
        (double*)malloc(spectra[iS].nMaxZ*nAgeStep*spectra[iS].nFlux*sizeof(double));
        init_templates_working(spectra + iS, NULL, NULL, -1);
        // Initialise special templates for birth cloud
        init_templates_special(spectra + iS, tBC, 1);
    }

    // Initialise mini templates
    int nSize = 0;
    int nMaxZ = spectra->nMaxZ;
    double *working;
    size_t totalSize = 0;
    int offsetWorking = 0;
    int offsetInBC = 0;
    int offsetOutBC = 0;
    int offsetWaves = 0;

    //// Compute size of working templates
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS)
        totalSize += targetSnap[iS];
    totalSize *= nMaxZ*MAGS_N_BANDS;
    //// Compute size of special templates
    totalSize += 2*MAGS_N_SNAPS*nMaxZ*MAGS_N_BANDS;
    ///  Compute size of wavelengths
    totalSize += 2*MAGS_N_BANDS;
    totalSize *= sizeof(double);
    ////
    working = (double*)malloc(totalSize);
    //// Copy working templates
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        nSize = targetSnap[iS]*nMaxZ*MAGS_N_BANDS;
        memcpy(working + offsetWorking, spectra[iS].working, nSize*sizeof(double));
        offsetWorking += nSize;
    }
    //// Copy special templates
    offsetInBC = offsetWorking;
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        nSize = nMaxZ*MAGS_N_BANDS;
        memcpy(working + offsetInBC, spectra[iS].inBC, nSize*sizeof(double));
        offsetInBC += nSize;
    }
    offsetOutBC = offsetInBC;
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        nSize = nMaxZ*MAGS_N_BANDS;
        memcpy(working + offsetOutBC, spectra[iS].outBC, nSize*sizeof(double));
        offsetOutBC += nSize;
    }
    //// Copy wavelengths (same at each target snapshot)
    offsetWaves = offsetOutBC;
    memcpy(working + offsetWaves, spectra->centreWaves, MAGS_N_BANDS*sizeof(double));
    offsetWaves += MAGS_N_BANDS;
    memcpy(working + offsetWaves, spectra->logWaves, MAGS_N_BANDS*sizeof(double));
    ////
    memcpy(miniSpectra->targetSnap, targetSnap, MAGS_N_SNAPS*sizeof(int));
    miniSpectra->nBeta = 0;
    miniSpectra->minZ = spectra->minZ;
    miniSpectra->maxZ = spectra->maxZ;
    miniSpectra->nMaxZ = nMaxZ;
    ////   -Find the interval for birth cloud
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS)
        miniSpectra->iAgeBC[iS] = \
        birth_cloud_interval(tBC, spectra[iS].ageStep, spectra[iS].nAgeStep);
    miniSpectra->totalSize = totalSize;
    miniSpectra->working = working;
    miniSpectra->inBC = working + offsetWorking;
    miniSpectra->outBC = working + offsetInBC;
    miniSpectra->centreWaves = working + offsetOutBC;
    miniSpectra->logWaves = working + offsetWaves;

    // Free full templates
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        free(spectra[iS].Z);
        free(spectra[iS].waves);
        free(spectra[iS].age);
        free(spectra[iS].raw);
        free(spectra[iS].nFilterWaves);
        free(spectra[iS].filterWaves);
        free(spectra[iS].filters);
        free(spectra[iS].integrated);
        free(spectra[iS].ready);
        free(spectra[iS].working);
        free(spectra[iS].inBC);
        free(spectra[iS].outBC);
        free(spectra[iS].centreWaves);
        free(spectra[iS].logWaves);
    }
}


void init_luminosities(double *inBCFlux, double *outBCFlux) {
    int iSF;
    int nSF = MAGS_N_SNAPS*MAGS_N_BANDS;

    for(iSF = 0; iSF < nSF; ++iSF) {
        inBCFlux[iSF] = TOL;
        outBCFlux[iSF] = TOL;
    }
}


void add_luminosities(double *pInBCFlux, double *pOutBCFlux, mini_sed_params_t *spectra,
                      int snapshot, double metals, double sfr) {
    /* Add luminosities when there is a burst
     *   -SFRs should be in a unit of M_solar/yr. However, one can convert the unit on
     *    final results rather than here in order to achieve better performance */

    // Compute integer metallicity
    int Z = (int)(metals*1000 - .5);
    if (Z < spectra->minZ)
        Z = spectra->minZ;
    else if (Z > spectra->maxZ)
        Z = spectra->maxZ;

    // Add luminosities
    int iA, iF, iS, iAgeBC;
    int offset;
    int nAgeStep;
    int nZF = spectra->nMaxZ*MAGS_N_BANDS;
    double *pWorking = spectra->working;
    double *pInBC = spectra->inBC;
    double *pOutBC = spectra->outBC;

    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        nAgeStep = spectra->targetSnap[iS];
        iA = nAgeStep - snapshot;
        if(iA < 0)
            continue;
        iAgeBC = spectra->iAgeBC[iS];
        if (iA > iAgeBC) {
            offset = (Z*nAgeStep + iA)*MAGS_N_BANDS;
            for(iF = 0; iF < MAGS_N_BANDS; ++iF)
                pOutBCFlux[iF] += sfr*pWorking[offset + iF];
        }
        else if (iA == iAgeBC) {
            offset = Z*MAGS_N_BANDS;
            for(iF = 0; iF < MAGS_N_BANDS; ++iF) {
                pInBCFlux[iF] += sfr*pInBC[offset + iF];
                pOutBCFlux[iF] += sfr*pOutBC[offset + iF];
            }
        }
        else {
            offset = (Z*nAgeStep + iA)*MAGS_N_BANDS;
            for(iF = 0; iF < MAGS_N_BANDS; ++iF)
                pInBCFlux[iF] += sfr*pWorking[offset + iF];
        }
        pWorking += nAgeStep*nZF;
        pInBC += nZF;
        pOutBC += nZF;
        pInBCFlux += MAGS_N_BANDS;
        pOutBCFlux += MAGS_N_BANDS;
    }
}


void merge_luminosities(double *inBCFluxTgt, double *outBCFluxTgt,
                        double *inBCFlux, double *outBCFlux) {
    int iSF;
    int nSF = MAGS_N_SNAPS*MAGS_N_BANDS;

    for(iSF = 0; iSF < nSF; ++iSF) {
        inBCFluxTgt[iSF] += inBCFlux[iSF];
        outBCFluxTgt[iSF] += outBCFlux[iSF];
    }
}


void get_magnitudes(double *mags, double *inBCFlux, double *outBCFlux) {
    /* Get magnitudes for only one snapshot */
    int iF;
    for(iF = 0; iF < MAGS_N_BANDS; ++iF)
        mags[iF] = M_AB(inBCFlux[iF] + outBCFlux[iF]);
}
#endif


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Galaxy properites related                                                   *
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
