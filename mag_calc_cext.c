#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

//#define SURFACE_AREA 1.1965e40 // 4*pi*(10 pc)**2 unit cm^2
//#define JANSKY(x) (3.34e4*(x)*(x))
#define INTGER_Z(Z, minZ, maxZ) (int)(1000*Z - 0.5)
#define M_AB(x) (-2.5*log10(x) + 8.9) // Convert Jansky to AB magnitude
#define TOL 1e-30 // Minimum Flux


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Basic Functions                                                             *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
clock_t g_sTime;


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


inline void report(int i, int tot) {
    int n = tot > 10 ? tot/10 : 1;
    if (i%n == 0) {
        printf("# %5.1f%% complete!\r", 100.*(i + 1)/tot);      
        fflush(stdout);
    }
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


void timing_start(char* text) {
    g_sTime = clock();
    printf("#**********************************************************\n");
    printf(text);
}

void timing_end(void) {
    float elapsedTime = (float)(clock() - g_sTime)/CLOCKS_PER_SEC;
    int min = (int)elapsedTime/60;
    printf("# 100%% complete!\n");
    printf("# Done!\n");
    printf("# Elapsed time: %d min %.5f sec\n", min, elapsedTime - min*60);
    printf("#**********************************************************\n\n");
}


inline int bisection_search(double a, double *x, double nX) {
    /* return idx such x[idx] <= a < x[idx + 1] 
     * a must be x[0] <= a < x[nX - 1]
     */
    int idx0 = 0;
    int idx1 = nX - 1;
    int idxMid;
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
    int idx0, idx1;
    if((xp < x[0]) || (xp > x[nPts - 1])) {
        printf("Error: Point %10.5e is beyond the interpolation region\n", xp);
        exit(0);
    }
    if (xp == x[nPts - 1])
        return y[nPts - 1];
    else {
        idx0 = bisection_search(xp, x, nPts);
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
        I += (waves[i] - waves[i - 1])*(y0 + y1)/2.;
        y0 = y1;
    }
    return I;
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
struct props {
    short index;
    float metals;
    float sfr;
};

struct prop_set {
   struct props *nodes;
   int nNode;
};


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Functions to process SEDs                                                   *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
// Struct for SED templates
struct sed_params {
    double *Z;
    int nZ;
    int minZ;
    int maxZ;
    double *waves;
    int nWaves;
    double *age;
    int nAge;
    double *data;
};


struct sed_params *g_rawSpectra = NULL;
struct sed_params *g_intSpectra = NULL;

struct sed_params *init_template(double *age, int nAge, 
                               double *waves, int nWaves, 
                               double *data, int nZ) {
    struct sed_params *spectra = malloc(sizeof(struct sed_params));
    spectra->age = age;
    spectra->nAge = nAge;
    spectra->waves = waves;
    spectra->nWaves = nWaves;
    spectra->data = data;
    spectra->nZ = nZ;
    return spectra;
}


void free_template(struct sed_params *spectra, int nRows) {
    free(spectra->age);
    free(spectra->waves);
    free(spectra->data);
}


void free_raw_spectra(void) {
    free(g_rawSpectra->age);
    free(g_rawSpectra->waves);
    free(g_rawSpectra->data);
    g_rawSpectra = NULL;
}


void free_int_spectra(void) {
    free(g_intSpectra->data);
    g_intSpectra = NULL;
}


void templates_time_integration(double *ageList, int nAgeList) {
    int iA, iW, iZ;
    double *pData;

    int nAge;
    double *age;
    int nWaves; 
    double *waves;
    int nZ;
    double *data;
    // Spectra after integration over time
    // The first dimension refers to metallicites and ages
    // The last dimension refers to wavelengths
    double *intData;
    
    if(g_intSpectra == NULL) {
        timing_start("# Integrate SED templates over time\n");
        nAge = g_rawSpectra->nAge;
        age = g_rawSpectra->age;
        nWaves = g_rawSpectra->nWaves; 
        waves = g_rawSpectra->waves;
        nZ = g_rawSpectra->nZ;
        data = g_rawSpectra->data;
        intData = (double*)malloc(nZ*nAgeList*nWaves*sizeof(double));
        for(iZ = 0; iZ < nZ; report(iZ++, nZ)) 
            for(iA = 0; iA < nAgeList; ++iA) {
                pData = intData + (iZ*nAgeList + iA)*nWaves;
                for(iW = 0; iW < nWaves; ++iW) {
                    if (iA == 0) 
                        // The first time step of SED templates is typicall not zero
                        // Here assumes that the templates is zero beween zero
                        // and the first time step
                        pData[iW] = trapz_table(data + (iZ*nWaves + iW)*nAge, age, nAge, 
                                                age[0], ageList[iA]);
                    else
                        pData[iW] = trapz_table(data + (iZ*nWaves + iW)*nAge, age, nAge, 
                                                ageList[iA - 1], ageList[iA]);
                }
            }
        g_intSpectra = init_template(ageList, nAgeList, waves, nWaves, intData, nZ);
        timing_end();
    }

}
 

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Functions of the dust model                                                 *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
struct dust_params {
    double tauUV_ISM;
    double nISM;
    double tauUV_BC;
    double nBC;
    double tBC;
};


inline double *dust_absorption(struct dust_params *dustArgs) {
    /* tBC: life time of the birth clound
     * nu: fraction of ISM dust absorption
     * tauUV: V-band absorption optical depth
     * nBC: power law index of tauBC
     * nISM: power law index of tauISM
     * 
     * Reference: da Cunha et al. 2008
     */
    int iA, iW, iZ, iAZ;
    double *pData;

    int nAge = g_intSpectra->nAge;
    double *age = g_intSpectra->age;
    int nWaves = g_intSpectra->nWaves;
    double *waves = g_intSpectra->waves;
    int nZ = g_intSpectra->nZ;
    double *data = malloc(nZ*nAge*nWaves*sizeof(double));
    memcpy(data, g_intSpectra->data, nZ*nAge*nWaves*sizeof(double));

    int nRawAge = g_rawSpectra->nAge;
    double *rawAge = g_rawSpectra->age;
    double *rawData = g_rawSpectra->data;

    int iAgeBC;
    double t0, t1;
    double ratio;

    double tauUV_ISM = dustArgs->tauUV_ISM;
    double nISM = dustArgs->nISM;
    double tauUV_BC = dustArgs->tauUV_BC;
    double nBC = dustArgs->nBC;
    double tBC = dustArgs->tBC;

    double *tauISM = malloc(nWaves*sizeof(double));
    double *tauBC = malloc(nWaves*sizeof(double));

    // Compute the optical depth of both the birth cloud and the ISM
    for(iW = 0; iW < nWaves; ++iW) {
        ratio = waves[iW]/1600.;
        tauISM[iW] = tauUV_ISM*pow(ratio, nISM);
        tauBC[iW] = tauUV_BC*pow(ratio, nBC);
    }   

    // Birth cloud part 
    // Find the time inverval containning the birth cloud
    if (tBC >= age[nAge - 1])
        iAgeBC = nAge;
    else if(tBC < age[0]) {
        iAgeBC = 0;
        t0 = rawAge[0];
        t1 = age[0];
    }
    else {
        iAgeBC = bisection_search(tBC, age, nAge) + 1;
        t0 = age[iAgeBC - 1];
        t1 = age[iAgeBC];
    }
    // t_s < tBC < t_s + dt
    if (iAgeBC < nAge)
        for(iZ = 0; iZ < nZ; ++iZ) {
            pData = data + (iZ*nAge + iAgeBC)*nWaves;
            for(iW = 0; iW < nWaves; ++iW) 
                pData[iW] = exp(-tauBC[iW]) \
                            *trapz_table(rawData + (iZ*nWaves + iW)*nRawAge, 
                                         rawAge, nRawAge, t0, tBC) \
                            + trapz_table(rawData + (iZ*nWaves + iW)*nRawAge, 
                                         rawAge, nRawAge, tBC, t1);
        }
    
    // tBC > t_s
    for(iA = 0; iA < iAgeBC; ++iA)
        for(iZ = 0; iZ < nZ; ++iZ) {
            pData = data + (iZ*nAge + iA)*nWaves;
            for(iW = 0; iW < nWaves; ++iW) 
                pData[iW] *= exp(-tauBC[iW]);
        }

    // ISM part
    for(iAZ = 0; iAZ < nAge*nZ; ++iAZ) {
        pData = data + iAZ*nWaves;
        for(iW = 0; iW < nWaves; ++iW) 
            pData[iW] *= exp(-tauISM[iW]);
    }
    return data;
}


void templates_working(double *fluxTmp, double z, double *filters, int nRest, int nObs, 
                       double *LyAbsorption, struct dust_params *dustArgs) {
    int iA, iW, iZ, iAZ, iF;
    double *pData;

    int nZ = g_intSpectra->nZ;
    int nAge = g_intSpectra->nAge;
    int nWaves = g_intSpectra->nWaves;
    double *waves = g_intSpectra->waves;

    double *data;

    //timing_start("# Process working SED templates\n");
    if (dustArgs != NULL)
        data = dust_absorption(dustArgs);
    else
        data = g_intSpectra->data;

    double *obsWaves = NULL;
    double *obsData = NULL;
    double *pObsData;

    if (nObs > 0) {
        // Transform everything to observer frame
        // Note the fluxes in this case is a function of wavelength
        // Therefore the fluxes has a factor of 1/(1 + z)
        obsWaves = (double*)malloc(nWaves*sizeof(double));
        obsData = (double*)malloc(nZ*nAge*nWaves*sizeof(double));
        for(iW = 0; iW < nWaves; ++iW)
            obsWaves[iW] = waves[iW]*(1. + z);
        for(iAZ = 0; iAZ < nZ*nAge; ++iAZ) {
            pData = data + iAZ*nWaves;
            pObsData = obsData + iAZ*nWaves;
            for(iW = 0; iW < nWaves; ++iW)
                pObsData[iW] = pData[iW]/(1. + z);           
        }
        if (LyAbsorption != NULL)
            // Add IGM absorption
             for(iAZ = 0; iAZ < nZ*nAge; ++iAZ) {
                pObsData = obsData + iAZ*nWaves;
                for(iW = 0; iW < nWaves; ++iW)
                    pObsData[iW] *= LyAbsorption[iW];
                }       
    }

    int nFilter = filters == NULL ? nWaves : nRest + nObs;
    double *pFilter;   

    // Spectra to be interploated along metallicities
    // The first dimension refers to filters/wavelengths and ages
    // Thw last dimension refers to metallicites
    double *refSpectra = malloc(nFilter*nAge*nZ*sizeof(double));
    // Output
    // The first dimension refers to metallicites
    // The second dimension refers to ages
    // The last dimension refers to filters/wavelengths

    if (filters != NULL) {
        // Intgrate SED templates over filters
        // Compute fluxes in rest frame filters
        //printf("# Compute fluxes in rest frame filters\n");
        for(iF = 0; iF < nRest; ++iF) {
            pFilter = filters + iF*nWaves;
            for(iA = 0; iA < nAge; ++iA) {
                pData = refSpectra + (iF*nAge+ iA)*nZ;
                for(iZ = 0; iZ < nZ; ++iZ)
                    pData[iZ] = trapz_filter(pFilter, data + (iZ*nAge + iA)*nWaves,
                                             waves, nWaves);
            }
        }
        //printf("# Compute fluxes in observer frame filters\n");
        for(iF = nRest; iF < nFilter; ++iF) {
            pFilter = filters + iF*nWaves;
            for(iA = 0; iA < nAge; ++iA) {
                pData = refSpectra + (iF*nAge+ iA)*nZ;
                for(iZ = 0; iZ < nZ; ++iZ)
                    pData[iZ] = trapz_filter(pFilter, obsData + (iZ*nAge+ iA)*nWaves, 
                                             obsWaves, nWaves);
            }
        }
    }
    else {
        if (nObs > 0) {
            // Tranpose the templates such that the last dimension is the metallicity
            for(iZ = 0; iZ < nZ; ++iZ) 
                for(iA = 0; iA < nAge; ++iA) {
                    pData = obsData + (iZ*nAge + iA)*nWaves;
                    for(iW = 0; iW < nWaves; ++iW)
                        refSpectra[(iW*nAge + iA)*nZ + iZ] = pData[iW];
                }
        }
        else {
            for(iZ = 0; iZ < nZ; ++iZ) 
                for(iA = 0; iA < nAge; ++iA) {
                    pData = data + (iZ*nAge + iA)*nWaves;
                    for(iW = 0; iW < nWaves; ++iW)
                        refSpectra[(iW*nAge + iA)*nZ + iZ] = pData[iW];
                }
        }
    }

    // Interploate SED templates along metallicities
    //printf("# Interpolate SED templates along metallicities\n");
    double *Z = g_rawSpectra->Z;
    int minZ = g_rawSpectra->minZ;
    int maxZ = g_rawSpectra->maxZ;
    pData = fluxTmp;
    for(iZ = minZ; iZ < maxZ + 1; ++iZ)
        for(iA = 0; iA < nAge; ++iA) 
            for(iF = 0; iF < nFilter; ++iF) 
                *pData++ = interp((iZ + 1.)/1000., Z, refSpectra + (iF*nAge+ iA)*nZ, nZ);

    if (nObs > 0) {
        free(obsWaves);
        free(obsData);
    }
    if (dustArgs != NULL)
        free(data);
    free(refSpectra); 
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Primary Functions                                                           *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
float *composite_spectra_cext(struct sed_params *rawSpectra,
                              struct prop_set *galProps, int nGal,
                              double z, double *ageList, int nAgeList,
                              double *filters, double* logWaves, int nRest, int nObs,
                              double *absorption, struct dust_params *dustArgs,
                              int outType) {
    int iF, iG, iP, iFG;
    double *pData;

    //Generate templates
    g_rawSpectra = rawSpectra;
    templates_time_integration(ageList, nAgeList);

    int nFilter = filters == NULL ? g_intSpectra->nWaves : nObs + nRest;
    double *fluxTmp = (double*)malloc((g_rawSpectra->maxZ + 1)*nAgeList*nFilter*sizeof(double));
    double *flux = malloc(nFilter*sizeof(double));
    float *output = malloc(nGal*nFilter*sizeof(float));
    float *pOutput = output;

    struct prop_set *pGalProps;
    struct props *pNodes;
    int nProg;

    double sfr;
    int metals;
    int minZ = rawSpectra->minZ;
    int maxZ = rawSpectra->maxZ;

    timing_start("# Compute magnitudes\n");
    if (dustArgs == NULL)
        templates_working(fluxTmp, z, filters, nRest, nObs, absorption, dustArgs);
    for(iG = 0; iG < nGal; report(iG++, nGal)) {
        // Initialise fluxes
        for(iF = 0; iF < nFilter; ++iF)
            flux[iF] = TOL;
        // Sum contributions from all progentiors
        pGalProps = galProps + iG;
        nProg = pGalProps->nNode;
        if (dustArgs != NULL)
            templates_working(fluxTmp, z, filters, nRest, nObs, absorption, dustArgs + iG);
        for(iP = 0; iP < nProg; ++iP) {
            pNodes = pGalProps->nodes + iP;
            sfr = pNodes->sfr;
            metals = (int)(pNodes->metals*1000 - .5);
            if (metals < minZ)
                metals = minZ;
            else if (metals > maxZ)
                metals = maxZ;
            pData = fluxTmp + (metals*nAgeList + pNodes->index)*nFilter;
            for(iF = 0 ; iF < nFilter; ++iF) {
                flux[iF] += sfr*pData[iF];
            }
        }
        // Store output
        for(iF = 0; iF < nFilter; ++iF) 
            *pOutput++ = (float)flux[iF];
    }
    free(fluxTmp);
    free(flux);

    if (outType == 0) {
        pOutput = output;
        for(iFG = 0; iFG < nFilter*nGal; iFG++) {
            *pOutput = M_AB(*pOutput);
            ++pOutput;
        }
        timing_end();
        return output;
    }
    else if (outType == 1) {
        timing_end();   
        return output;
    }
    
    // Fit UV slopes
    int nR = 3;
    struct linResult result;

    output = (float*)realloc(output, (nFilter + nR)*nGal*sizeof(float));
    pOutput = output + nFilter*nGal;
    float *pFit = output;

    int nFit = nFilter - 1;
    double *logf = malloc(nFit*sizeof(double));

    for(iG = 0; iG < nGal; ++iG) {
        for(iF = 0; iF < nFit; ++iF) 
            logf[iF] = log(pFit[iF]);
        pFit += nFilter;

        //printf("waves = %.1f, logf = %.1f\n", logWaves[1], logf[1]);
        result = linregress(logWaves, logf, nFit);
        pOutput[0] = (float)result.slope;
        pOutput[1] = (float)result.intercept;
        pOutput[2] = (float)result.R;
        pOutput += nR;
        //printf("Slope = %.1f\n", result.slope);
    }
        
    // Convert to AB magnitude
    pOutput = output + nFit;
    for(iG = 0; iG < nGal; ++iG) {
        *pOutput = M_AB(*pOutput);
        pOutput += nFilter;
    }

    timing_end();
    return output;
}

