#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define MAX_NODE 100000 // Max length of galaxy merger tree

#define NMETALS 5 // Number of input metallicity
#define NUM_Z 40 // Number of interpolated metallicity
#define MAX_Z 39 // Maximum metallicity index
#define MIN_Z 0 // Minimum metallicity index

#define SURFACE_AREA 1.1965e40 // 4*pi*(10 pc)**2 unit cm^2
#define JANSKY(x) (3.34e4*(x)*(x))
#define M_AB(x) (-2.5*log10(x) + 8.9) // Convert Jansky to AB magnitude
#define TOL 1e-50 // Minimum Flux

clock_t sTime;
// Meraxes output
int **firstProgenitor = NULL;
int **nextProgenitor = NULL;
float **galMetals = NULL;
float **galSFR = NULL;
// Variable for galaxy merger tree
struct node {
    short snap;
    short metals;
    float sfr;
};
//
char sedAge[] = "/lustre/projects/p113_astro/yqiu/magcalc/input/sed_age.bin";
char sedWaves[] = "/lustre/projects/p113_astro/yqiu/magcalc/input/sed_waves.bin";
char *sedTemplates[NMETALS] = {"/lustre/projects/p113_astro/yqiu/magcalc/input/sed_0.001.bin",
                               "/lustre/projects/p113_astro/yqiu/magcalc/input/sed_0.004.bin",
                               "/lustre/projects/p113_astro/yqiu/magcalc/input/sed_0.008.bin",
                               "/lustre/projects/p113_astro/yqiu/magcalc/input/sed_0.020.bin",
                               "/lustre/projects/p113_astro/yqiu/magcalc/input/sed_0.040.bin"};



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
    int n = tot > 10 ? tot/10:1;
    if (i%n == 0) {
        printf("# %5.1f%% complete!\r", 100.*(i + 1)/tot);      
        fflush(stdout);
    }
}


double **malloc_2d_double(int nRow, int nCol) {
    int i;
    double **target;
    target = (double**)malloc(nRow*sizeof(double*));
    for(i = 0; i < nRow; ++i) 
        target[i] = (double*)malloc(nCol*sizeof(double));
    return target;
}


void free_2d_double(double **p, int nRow) {
    int i;
    for(i = 0; i < nRow; ++i)
        free(p[i]);
    free(p);
}


void timing_start(char* text) {
    sTime = clock();
    printf("#**********************************************************\n");
    printf(text);
}

void timing_end(void) {
    float elapsedTime = (float)(clock() - sTime)/CLOCKS_PER_SEC;
    int min = (int)elapsedTime/60;
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


// Struct for SED templates
struct template {
    double *age;
    int nAge;
    double *waves;
    int nWaves;
    // The first dimension refers to metallicities and ages
    // The list dimension refers to wavelengths
    double **data;
    int nZ;
};


void init_template(struct template *spectra, double *age, int nAge, 
                   double *waves, int nWaves, double **data, int nZ) {
    spectra->age = age;
    spectra->nAge = nAge;
    spectra->waves = waves;
    spectra->nWaves = nWaves;
    spectra->data = data;
    spectra->nZ = nZ;
}


void read_sed_templates(struct template *spectra) {
    /* SED templates must be normalised to 1 M_sun with unit erg/s/A 
     * Wavelengths must be in a unit of angstrom
     * Ages must be in a unit of year 
     * The first dimension of templates must be wavelengths
     * The last dimension of templates must be ages
     */
    FILE *fp;
    int iA, iW, iZ, iWZ;
    double *pData;

    int nAge;
    double *age;
    int nWaves;
    double *waves;
    int nZ;
    double **data;

    timing_start("# Read SED templates\n");
    fp = open_file(sedAge, "r");
    fread(&nAge, sizeof(int), 1, fp);
    age = (double*)malloc(nAge*sizeof(double));
    fread(age, sizeof(double), nAge, fp);
    fclose(fp);

    fp = open_file(sedWaves, "r");
    fread(&nWaves, sizeof(int), 1, fp);
    waves = (double*)malloc(nWaves*sizeof(double));
    fread(waves, sizeof(double), nWaves, fp);
    fclose(fp);

    nZ = NMETALS;
    data = malloc_2d_double(nZ*nWaves, nAge);
    for(iZ = 0; iZ < nZ; ++iZ) {
        fp = open_file(sedTemplates[iZ], "r");
        for(iW = 0; iW < nWaves; ++iW) 
            fread(data[iZ*nWaves + iW], sizeof(double), nAge, fp);
        fclose(fp);
    }
    
    // Convert flux to flux density at 10 pc
    for(iWZ = 0; iWZ < nZ*nWaves; ++iWZ) {
        pData = data[iWZ];
        for(iA = 0; iA < nAge; ++iA)
            pData[iA] /= SURFACE_AREA;
    }


    init_template(spectra, age, nAge, waves, nWaves, data, nZ);
    timing_end();
}


void free_template(struct template *spectra) {
    free(spectra->age);
    free(spectra->waves);
    free_2d_double(spectra->data, spectra->nZ*spectra->nWaves);
}

/*
double *init_templates_sp(double *ageList, int nAgeList) {
    int i, j, k;
    double refMetals[NMETALS] = {0., 3., 7., 19., 39.};
    struct template rawSpectra;
    int nWaves = 1221;
    double *waves;
    // Spectra to be interploated along metallicities
    // The first dimension refers to ages and wavelengths
    // The last dimension refers to metallicites
    double **refSpectra;
    // Output
    // The first dimension refers to metallicites
    // The second dimension refers to ages
    // The last dimension refers to wavelengths
    double *fluxTmp;
    double *pData;

    read_sed_templates(&rawSpectra);
    timing_start("# Process SED templates\n");
    nWaves = rawSpectra.nWaves;
    waves = (double*)malloc(nWaves*sizeof(double));
    memcpy(waves, rawSpectra.waves, nWaves*sizeof(double));
    // Integrate raw SED templates over time
    printf("# Integrate SED templates over time\n");
    refSpectra = malloc_2d_double(nAgeList*nWaves, NMETALS);
    
    for(i = 0; i < nAgeList; report(i++, nAgeList)) 
        for(j = 0; j < nWaves; ++j)  {
            pData = refSpectra[i*nWaves + j];
            for(k = 0; k < NMETALS; ++k) {
                if (i == 0) {
                    // The first time step of SED templates is typicall not zero
                    // Here assumes that the templates is constant beween zero
                    // and the first time step
                    //pData[k] = rawSpectra.data[k*nWaves + j][0]*rawSpectra.age[0];
                    pData[k] = trapz_table(rawSpectra.data[k*nWaves + j],
                                           rawSpectra.age, rawSpectra.nAge, 
                                           rawSpectra.age[0], ageList[i]);
                }
                else
                    pData[k] = trapz_table(rawSpectra.data[k*nWaves + j], 
                                           rawSpectra.age, rawSpectra.nAge, 
                                           ageList[i - 1], ageList[i]);
                // Convert to rest frame flux                       
                // pData[k] *= JANSKY(waves[j]);
            }
        }
    free_template(&rawSpectra);
    // Interpolate SED templates along metallicites
    printf("# Interpolate SED templates along metallicities\n");
    fluxTmp = (double*)malloc(NUM_Z*nAgeList*nWaves*sizeof(double));
    pData = fluxTmp;
    for(i = 0; i < NUM_Z; ++i)
        for(j = 0; j < nAgeList; ++j) 
            for(k = 0; k < nWaves; ++k)
                *pData++ = interp((double)i, refMetals, refSpectra[j*nWaves + k], NMETALS);
        
    free(refSpectra);
    timing_end();
    return fluxTmp;
}
*/

void templates_time_integration(struct template *spectra, struct template *rawSpectra,
                                double *ageList, int nAgeList) {
    int i, j, k;
    double *pData;

    int nAge = rawSpectra->nAge;
    double *age = rawSpectra->age;
    int nWaves = rawSpectra->nWaves; 
    double *waves = rawSpectra->waves;
    int nZ = rawSpectra->nZ;
    double **data = rawSpectra->data;
    // Spectra after integration over time
    // The first dimension refers to metallicites and ages
    // The last dimension refers to wavelengths
    double **intData = malloc_2d_double(NMETALS*nAgeList, nWaves);

    timing_start("# Process SED templates\n");
    // Integrate raw SED templates over time
    printf("# Integrate SED templates over time\n");
    for(i = 0; i < NMETALS; report(i++, NMETALS)) 
        for(j = 0; j < nAgeList; ++j) {
            pData = intData[i*nAgeList + j];
            for(k = 0; k < nWaves; ++k) {
                if (j == 0) {
                    // The first time step of SED templates is typicall not zero
                    // Here assumes that the templates is constant beween zero
                    // and the first time step
                    //pData[k] = data[i*nWaves + k][0]*age[0];
                    pData[k] = trapz_table(data[i*nWaves + k], age, nAge, 
                                           age[0], ageList[j]);

                }
                else
                    pData[k] = trapz_table(data[i*nWaves + k], age, nAge, 
                                           ageList[j - 1], ageList[j]);
            }
        }

    // Update templates
    double *waves2 = malloc(nWaves*sizeof(double));
    memcpy(waves2, waves, nWaves*sizeof(double));
    init_template(spectra, ageList, nAgeList, waves2, nWaves, intData, nZ);
}
 

void dust_absorption(struct template *spectra, struct template *rawSpectra,
                     double tBC, double mu, double tauV, double nBC, double nISM) {
    /* tBC: life time of the birth clound
     * nu: fraction of ISM dust absorption
     * tauV: V-band absorption optical depth
     * nBC: power law index of tauBC
     * nISM: power law index of tauISM
     * 
     * Reference: da Cunha et al. 2008
     */
    int iA, iW, iZ, iAZ;
    double *pData;
    int nAge = spectra->nAge;
    double *age = spectra->age;
    int nWaves = spectra->nWaves;
    double *waves = spectra->waves;
    int nZ = spectra->nZ;
    double **data = spectra->data;

    int nRawAge = rawSpectra->nAge;
    double *rawAge = rawSpectra->age;
    double **rawData = rawSpectra->data;

    double t0, t1;

    double ratio;
    int iAgeBC;
    double tauV_BC = (1 - mu)*tauV;
    double tauV_ISM = mu*tauV;
    double *tauISM = malloc(nWaves*sizeof(double));
    double *tauBC = malloc(nWaves*sizeof(double));

    // Compute the optical depth of both the birth cloud and the ISM
    for(iW = 0; iW < nWaves; ++iW) {
        ratio = waves[iW]/5500.;
        tauISM[iW] = tauV_ISM*pow(ratio, nISM);
        tauBC[iW] = tauV_BC*pow(ratio, nBC);
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
            pData = data[iZ*nAge + iAgeBC];
            for(iW = 0; iW < nWaves; ++iW) 
                pData[iW] = exp(-tauBC[iW]) \
                            *trapz_table(rawData[iZ*nWaves + iW], rawAge, nRawAge, t0, tBC) \
                            + trapz_table(rawData[iZ*nWaves + iW], rawAge, nRawAge, tBC, t1);
        }
    
    // tBC > t_s
    for(iA = 0; iA < iAgeBC; ++iA)
        for(iZ = 0; iZ < nZ; ++iZ) {
            pData = data[iZ*nAge + iA];
            for(iW = 0; iW < nWaves; ++iW) 
                pData[iW] *= exp(-tauBC[iW]);
        }

    // ISM part
    for(iAZ = 0; iAZ < nAge*nZ; ++iAZ) {
        pData = data[iAZ];
        for(iW = 0; iW < nWaves; ++iW) 
            pData[iW] *= exp(-tauISM[iW]);
    }

}


double *templates_working(struct template *spectra, double z,
                          double *filters, int nRest, int nObs, double *absorption) {
    int iA, iW, iZ, iAZ, iF;
    double *pData;

    int nAge = spectra->nAge;
    int nWaves = spectra->nWaves;
    double *waves = spectra->waves;
    int nZ = spectra->nZ;
    double **data = spectra->data;

    double refMetals[NMETALS] = {0., 3., 7., 19., 39.};
    int nFilter = filters == NULL ? nWaves:nRest + nObs;
    double *pFilter;   

    // Spectra to be interploated along metallicities
    // The first dimension refers to filters/wavelengths and ages
    // Thw last dimension refers to metallicites
    double **refSpectra = malloc_2d_double(nFilter*nAge, nZ);
    // Output
    // The first dimension refers to metallicites
    // The second dimension refers to ages
    // The last dimension refers to filters/wavelengths
    double *fluxTmp;

    if (filters != NULL) {
        // Intgrate SED templates over filters
        // Compute fluxes in rest frame filters
        printf("# Compute fluxes in rest frame filters\n");
        for(iF = 0; iF < nRest; ++iF) {
            pFilter = filters + iF*nWaves;
            for(iA = 0; iA < nAge; ++iA) {
                pData = refSpectra[iF*nAge+ iA];
                for(iZ = 0; iZ < nZ; ++iZ)
                    pData[iZ] = trapz_filter(pFilter, data[iZ*nAge + iA], waves, nWaves);
            }
        }
    }

    if (nObs > 0) {
        // Transform everything to observer frame
        // Note the fluxes in this case is a function of wavelength
        // Therefore the fluxes has a factor of 1/(1 + z)
        for(iW = 0; iW < nWaves; ++iW)
            waves[iW] *= 1. + z;
        for(iAZ = 0; iAZ < nZ*nAge; ++iAZ) {
            pData = data[iAZ];
            for(iW = 0; iW < nWaves; ++iW)
                pData[iW] /= 1. + z;           
        }
        if (absorption != NULL)
            // Add IGM absorption
             for(iAZ = 0; iAZ < nZ*nAge; ++iAZ) {
                pData = data[iAZ];
                for(iW = 0; iW < nWaves; ++iW)
                    pData[iW] *= absorption[iW];
                }       
    }
    
    if (filters != NULL) {
        // Compute fluxes in observer frame filters
        printf("# Compute fluxes in observer frame filters\n");
        for(iF = nRest; iF < nFilter; ++iF) {
            pFilter = filters + iF*nWaves;
            for(iA = 0; iA < nAge; ++iA) {
                pData = refSpectra[iF*nAge+ iA];
                for(iZ = 0; iZ < nZ; ++iZ)
                    pData[iZ] = trapz_filter(pFilter, data[iZ*nAge+ iA], waves, nWaves);
            }
        }
    }

    if (filters == NULL)
        // Tranpose the templates such that the last dimension is the metallicity
        for(iZ = 0; iZ < nZ; ++iZ) 
            for(iA = 0; iA < nAge; ++iA) {
                pData = data[iZ*nAge + iA];
                for(iW = 0; iW < nWaves; ++iW)
                    refSpectra[iW*nAge + iA][iZ] = pData[iW];
            }
    
    // Interploate SED templates along metallicities
    printf("# Interpolate SED templates along metallicities\n");
    fluxTmp = (double*)malloc(NUM_Z*nAge*nFilter*sizeof(double));
    pData = fluxTmp;
    for(iZ = 0; iZ < NUM_Z; ++iZ)
        for(iA = 0; iA < nAge; ++iA) 
            for(iF = 0; iF < nFilter; ++iF) 
                *pData++ = interp((double)iZ, refMetals, refSpectra[iF*nAge+ iA], nZ);

    free_2d_double(refSpectra, nFilter*nAge); 
    timing_end();
    return fluxTmp;   
}


void trace_progenitors(int snap, int galIdx, struct node *branch, int *pNProg) {
    int metals;
    float sfr;
    if (galIdx >= 0) {
        sfr = galSFR[snap][galIdx];
        if (sfr > 0.) {
            *pNProg += 1;
            if (*pNProg >= MAX_NODE) {
                printf("Error: Number of progenitors exceeds MAX_NODE\n");
                exit(0);
            }
            metals = (int)(galMetals[snap][galIdx]*1000 - .5);
            if (metals < MIN_Z)
                metals = MIN_Z;
            else if (metals > MAX_Z)
                metals = MAX_Z;
            branch[*pNProg].snap = snap;
            branch[*pNProg].metals = metals;
            branch[*pNProg].sfr = sfr;
            //printf("snap %d, metals %d, sfr %.3f\n", snap, metals, sfr);
        }
        trace_progenitors(snap - 1, firstProgenitor[snap][galIdx], branch, pNProg);
        trace_progenitors(snap, nextProgenitor[snap][galIdx], branch, pNProg);
    }
}


inline int trace_merger_tree(int snap, int galIdx, struct node *branch) {
    int nProg = -1;
    int metals;
    float sfr = galSFR[snap][galIdx];
    
    if (sfr > 0.) {
        ++nProg;
        metals = (int)(galMetals[snap][galIdx]*1000 - .5);
        if (metals < MIN_Z)
            metals = MIN_Z;
        else if (metals > MAX_Z)
            metals = MAX_Z;
        branch[nProg].snap = snap;
        branch[nProg].metals = metals;
        branch[nProg].sfr = sfr;
    }
    trace_progenitors(snap - 1, firstProgenitor[snap][galIdx], branch, &nProg);
    ++nProg;
    if (nProg == 0) {
        printf("Warning: snapshot %d, index %d\n", snap, galIdx);
        printf("         the star formation rate is zero throughout the histroy\n");
    }
    return nProg;
}


inline void compute_spectrum(double *spectrum, int cSnap,
                             struct node *branch, int nProg, struct template *spectra) {
    int i, j;
    int snap;
    double *pData;
    double sfr;
    for(i = 0; i < spectra->nWaves; ++i)
        *(spectrum + i) = 0.;

    for(i = 0; i < nProg; ++i) {
        snap = branch[i].snap;
        pData = spectra->data[branch[i].metals*spectra->nAge + cSnap - branch[i].snap];
        sfr = branch[i].sfr;
        for(j = 0; j < spectra->nWaves; ++j) 
            *(spectrum + j) += sfr*pData[j];    
        
    }    
}


float *composite_spectra_cext(double z, int tSnap,
                              int *indices, int nGal,
                              double *ageList, int nAgeList,
                              double *filters, int nRest, int nObs,
                              double *absorption,
                              int dust, double tBC, double mu, 
                              double tauV, double nBC, double nISM) {
    int iF, iG, iP, iFG;
    double *pData;
    
    //Generate templates
    struct template rawSpectra;
    struct template spectra;
    read_sed_templates(&rawSpectra);
    templates_time_integration(&spectra, &rawSpectra, ageList, nAgeList);
    if (dust)
        dust_absorption(&spectra, &rawSpectra, tBC, mu, tauV, nBC, nISM);
    
    double *fluxTmp = templates_working(&spectra, z, filters, nRest, nObs, absorption);
    
    int nFilter = filters == NULL ? spectra.nWaves:nObs + nRest;
    double *flux = malloc(nFilter*sizeof(double));
    float *output = malloc(nGal*nFilter*sizeof(float));
    float *pOutput = output;

    int nProg;
    struct node branch[MAX_NODE];
    double sfr;
    
    timing_start("# Compute magnitudes\n");
    for(iG = 0; iG < nGal; report(iG++, nGal)) {
        nProg = trace_merger_tree(tSnap, indices[iG], branch);
        // Initialise fluxes
        for(iF = 0; iF < nFilter; ++iF)
            flux[iF] = TOL;
        // Sum contributions from all progentiors
        for(iP = 0; iP < nProg; ++iP) {
            sfr = branch[iP].sfr;    
            pData = fluxTmp + (branch[iP].metals*nAgeList + tSnap - branch[iP].snap)*nFilter;
            for(iF = 0 ; iF < nFilter; ++iF) {
                flux[iF] += sfr*pData[iF];
            }
        }
        // Store output
        for(iF = 0; iF < nFilter; ++iF) 
            *pOutput++ = (float)flux[iF];
        
    }
    pOutput = output;
    if (filters != NULL)
        for(iFG = 0; iFG < nFilter*nGal; iFG++) {
            *pOutput = M_AB(*pOutput);
            ++pOutput;
        }

    free(flux);
    free(fluxTmp);
    timing_end();
    return output;
}



