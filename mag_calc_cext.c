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
#define NFILTER 10 // Number of input filters
#define MAX_FILTER 100 // Max Number of filters

#define SURFACE_AREA 1.1965e40 // 4*pi*(10 pc)**2 unit cm^2
#define TOL 1e-50

time_t sTime;
// Meraxes output
int **pFirstProgenitor = NULL;
int **pNextProgenitor = NULL;
float **pMetals = NULL;
float **pSFR = NULL;
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
// Struct for SED templates
struct template {
    int nAge;
    double *age;
    int nWaves;
    double *waves;
    int nZ;
    // The first dimension refers to metallicities and ages
    // The list dimension refers to wavelengths
    double **data;
};


void free_template(struct template *spectra) {
    int i;
    int nRow = spectra->nZ*spectra->nAge;
    free(spectra->age);
    free(spectra->waves);
    for(i = 0; i < nRow; ++i)
        free(spectra->data[i]);   
    free(spectra->data);
}


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


void report(int i, int tot) {
    if((tot > 10) && i%(tot/10) == 0) {
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
    sTime = time(NULL);
    printf("#*******************************************************************************\n");
    printf(text);
}

void timing_end(void) {
    int elapsedTime;
    elapsedTime = (int)difftime(time(NULL), sTime);
    printf("# Done!\n");
    printf("# Elapsed time: %d min %d sec\n", elapsedTime/60, elapsedTime%60);
    printf("#*******************************************************************************\n\n");
}


inline double interp(double xp, double *x, double *y, int nPts) {
    /* Interpolate a given points */
    int idx0, idx1, idxMid;
    if((xp < x[0]) || (xp > x[nPts - 1])) {
        printf("Error: The given point %10.5e is outside of the interpolation region\n", xp);
        exit(0);
    }
    idx0 = 0;
    idx1 = nPts - 1;
    while(idx1 - idx0 > 1) {
        idxMid = (idx0 + idx1)/2;
        if(xp > x[idxMid])
            idx0 = idxMid;
        else if(xp < x[idxMid])
            idx1 = idxMid;
        else
            return y[idxMid];
    }
    return y[idx0] + (y[idx1] - y[idx0])*(xp - x[idx0])/(x[idx1] - x[idx0]);
}


inline double trapz_table(double *y, double *x, int nPts, double a, double b) {
    /* Integrate tabular data from a to b */
    int i;
    int idx0 = 0;
    int idx1 = nPts - 1;
    int idxMid;
    double ya, yb;
    double I;
    if (x[0] > a) {
        printf("Error: Integration range %10.5e is outside the tabular data\n", a);
        exit(0);
    }
    if (x[nPts - 1] < b) {
        printf("Error: Integration range %10.5e is outside the tabular data\n", b); 
        exit(0);
    }
    if (a > b) {
        printf("Error: a must be smaller than b\n");
        exit(0);
    }
    // Use bisection to search the interval that contains a
    while(idx1 - idx0 > 1) {
        idxMid = (idx0 + idx1)/2;
        if (a > x[idxMid])
            idx0 = idxMid;
        else if (a < x[idxMid])
            idx1 = idxMid;
        else
            break;
    }

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


inline double trapz(double *y, double *x, int nPts) {
    int i;
    double I = 0.;
    for(i = 1; i < nPts; ++i)
        I += (x[i] - x[i - 1])*(y[i] + y[i - 1])/2.;
    return I;
}


void read_sed_templates(struct template *spectra) {
    /* SED templates must be normalised to 1 M_sun with unit erg/s/A 
     * Wavelength must be in a unit of angstrom
     * Age must be in a unit of year 
     */
    FILE *fp;
    int i, j;

    timing_start("# Read SED templates\n");
    fp = open_file(sedAge, "r");
    fread(&spectra->nAge, sizeof(int), 1, fp);
    spectra->age = (double*)malloc(spectra->nAge*sizeof(double));
    fread(spectra->age, sizeof(double), spectra->nAge, fp);
    fclose(fp);

    fp = open_file(sedWaves, "r");
    fread(&spectra->nWaves, sizeof(int), 1, fp);
    spectra->waves = (double*)malloc(spectra->nWaves*sizeof(double));
    fread(spectra->waves, sizeof(double), spectra->nWaves, fp);
    fclose(fp);

    spectra->nZ = NMETALS;
    spectra->data = malloc_2d_double(spectra->nZ*spectra->nAge, spectra->nWaves);
    for(i = 0; i < NMETALS; ++i) {
        fp = open_file(sedTemplates[i], "r");
        for(j = 0; j < spectra->nAge; ++j) 
            fread(spectra->data[i*spectra->nAge + j], sizeof(double), spectra->nWaves, fp);
        fclose(fp);
    }

    timing_end();
}


void get_integrand(double *integrand, struct template *spectra, int idxZ, int idxW) {
    int i;
    for(i = 0; i < spectra->nAge; ++i) 
        integrand[i] = spectra->data[idxZ*spectra->nAge + i][idxW];
}


void integrate_sed_templates(double *ageList, int nAgeList, 
                             struct template *spectra, struct template *rawSpectra) {
    int i, j, k;
    double refMetals[NMETALS] = {0., 3., 7., 19., 39.};
    double **refSpectra;
    double *pData;
    double *integrand = (double*)malloc(rawSpectra->nAge*sizeof(double));

    timing_start("# Process SED templates\n");
    spectra->nAge = nAgeList;
    spectra->age = (double*)malloc(nAgeList*sizeof(double));
    memcpy(spectra->age, ageList, nAgeList*sizeof(double));

    spectra->nWaves = rawSpectra->nWaves;
    spectra->waves = (double*)malloc(spectra->nWaves*sizeof(double));
    memcpy(spectra->waves, rawSpectra->waves, rawSpectra->nWaves*sizeof(double));
    
    spectra->nZ = NUM_Z;
    refSpectra = malloc_2d_double(spectra->nAge*spectra->nWaves, NMETALS);
    
    for(i = 0; i < spectra->nAge; ++i, report(i, spectra->nAge)) 
        for(j = 0; j < spectra->nWaves; ++j) {
            // In this case the first dimension of data represents both age s
            // and wavelengths 
            // The second dimension of data represents metallicities
            // This makes it easier to interpolate metallicities
            pData = refSpectra[i*spectra->nWaves + j];
            for(k = 0; k < NMETALS; ++k) {
                get_integrand(integrand, rawSpectra, k, j);
                if (i == 0) {
                    // The first time step of SED templates is typicall not zero
                    // Here assumes that the templates is constant beween zero
                    // and the first time step
                    pData[k] = integrand[0]*rawSpectra->age[0];
                    pData[k] += trapz_table(integrand, rawSpectra->age, rawSpectra->nAge, 
                                            rawSpectra->age[0], ageList[i]);
                }
                else
                    pData[k] = trapz_table(integrand, rawSpectra->age, rawSpectra->nAge, 
                                           ageList[i - 1], ageList[i]);
            }
        }

    spectra->data = malloc_2d_double(spectra->nZ*spectra->nAge, spectra->nWaves);
    for(i = 0; i < spectra->nZ; ++i)
        for(j = 0; j < spectra->nAge; ++j) {
            pData = spectra->data[i*spectra->nAge + j];
            for(k = 0; k < spectra->nWaves; ++k) 
                pData[k] = interp((double)i, refMetals, 
                                  refSpectra[j*spectra->nWaves + k], NMETALS);
        }
    
    printf("# Interpolate SED templates in terms of metallicity\n");
    timing_end();
    
    free(refSpectra);
    free(integrand);
}


void init_templates_ph(double **fluxTmp, double z,
                       double *ageList, int nAgeList, 
                       double *filters, int nRest, int nObs,
                       double *absorption) {

    int i, j, k;
    double refMetals[NMETALS] = {0., 3., 7., 19., 39.};
 
    int nWaves;
    double *waves;
    struct template rawSpectra;
    // Spectra that are integrated over each time step
    // The first dimension refers to metallicites and ages
    // The last dimension refers to wavelengths
    double **spectra1;    
    // Spectra that are integrated over filters
    // The first dimension refers to filters and ages
    // The last dimension refers to metallicites
    double **spectra2;
    double *pData;
    double *integrand;

    read_sed_templates(&rawSpectra);
    nWaves = rawSpectra.nWaves;
    integrand = (double*)malloc(rawSpectra.nAge*sizeof(double));

    timing_start("# Process SED templates\n");
    spectra1 = malloc_2d_double(NMETALS*nAgeList, nWaves);
    
    // Integrate raw SED templates over time
    for(i = 0; i < NMETALS; ++i, report(i, NMETALS)) 
        for(j = 0; j < nAgeList; ++j) {
            pData = spectra1[i*nAgeList + j];
            for(k = 0; k < nWaves; ++k) {
                get_integrand(integrand, &rawSpectra, k, j);
                if (i == 0) {
                    // The first time step of SED templates is typicall not zero
                    // Here assumes that the templates is constant beween zero
                    // and the first time step
                    pData[k] = integrand[0]*rawSpectra.age[0];
                    pData[k] += trapz_table(integrand, rawSpectra.age, rawSpectra.nAge, 
                                            rawSpectra.age[0], ageList[i]);
                }
                else
                    pData[k] = trapz_table(integrand, rawSpectra.age, rawSpectra.nAge, 
                                           ageList[i - 1], ageList[i]);
            }
        }
    free(integrand);
    free_template(&rawSpectra);
    
    // Intgrate SED templates over filters
    printf("# Integrate SED templates over filters\n");

    waves = (double*)malloc(nWaves*sizeof(double));
    memcpy(waves, rawSpectra.waves, nWaves);
    spectra2 = malloc_2d_double((nRest + nObs)*nAgeList, NMETALS);
    integrand = (double*)malloc(nWaves*sizeof(double));

    // Compute rest frame flux
    // Unit Jy
    for(i = 0; i < NMETALS; ++i)
        for(j = 0; j < nAgeList; ++j) {
            pData = spectra1[i*nAgeList + j];
            for(k = 0; i < nWaves; ++k)
                pData[i] *= 3.34e4*waves[i]*waves[i]/SURFACE_AREA; 
        }
    for(i = 0; i < nRest; ++i) {
        pFilter = filters + i*nWaves;
        for(j = 0; j < nAgeList; ++j) {

        for(j = 0; j < nWaves; ++j) {
            fFlux[j] = flux[j]*pFilter[j];
            if (fFlux[j] < TOL)
                fFlux[j] = TOL;
        }
        mags[i] = (float)(-2.5*log10(trapz(fFlux, waves, nWaves)) + 8.9);
    }
    // Compute observed frame magnitudes
    for(i = 0; i < nWaves; ++i)
        // Unit Jy
        flux[i] = 3.34e4*obsWaves[i]*obsWaves[i]*spectrum[i]/SURFACE_AREA/(1. + z)*absorption[i]; 
    for(i = 0; i < nObs; ++i) {
        pFilter = filters + (nRest + i)*nWaves;
        for(j = 0; j < nWaves; ++j) {
            fFlux[j] = flux[j]*pFilter[j];
            if (fFlux[j] < TOL)
                fFlux[j] = TOL;
        }
        mags[nRest + i] = (float)(-2.5*log10(trapz(fFlux, obsWaves, nWaves)) + 8.9);
 
    free(waves);

    // Interploate SED templates along metallicites
    printf("# Interpolate SED templates in terms of metallicities\n");
    spectra2 = malloc_2d_double(NUM_Z*nAgeList, nWaves);
    for(i = 0; i < NUM_Z; ++i)
        for(j = 0; j < nAgeList; ++j) {
            pData = spectra2[i*nAgeList + j];
            for(k = 0; k < nWaves; ++k) 
                pData[k] = interp((double)i, refMetals, spectra1[j*nWaves + k], NMETALS);
        }
    free_2d_double(spectra1, nAgeList*nWaves);   
    timing_end();
    

    free_2d_double(spectra2, NUM_Z*nAgeList);

}


void trace_progenitors(int snap, int galIdx, struct node *branch, int *pNProg) {
    int metals;
    float sfr;
    if (galIdx >= 0) {
        sfr = pSFR[snap][galIdx];
        if (sfr > 0.) {
            *pNProg += 1;
            if (*pNProg >= MAX_NODE) {
                printf("Error: Number of progenitors exceeds MAX_NODE\n");
                exit(0);
            }
            metals = (int)(pMetals[snap][galIdx]*1000 - .5);
            if (metals < MIN_Z)
                metals = MIN_Z;
            else if (metals > MAX_Z)
                metals = MAX_Z;
            branch[*pNProg].snap = snap;
            branch[*pNProg].metals = metals;
            branch[*pNProg].sfr = sfr;
            //printf("snap %d, metals %d, sfr %.3f\n", snap, metals, sfr);
        }
        trace_progenitors(snap - 1, pFirstProgenitor[snap][galIdx], branch, pNProg);
        trace_progenitors(snap, pNextProgenitor[snap][galIdx], branch, pNProg);
    }
}


inline int trace_merger_tree(int snap, int galIdx, struct node *branch) {
    int nProg = -1;
    int metals;
    float sfr = pSFR[snap][galIdx];
    
    if (sfr > 0.) {
        ++nProg;
        metals = (int)(pMetals[snap][galIdx]*1000 - .5);
        if (metals < MIN_Z)
            metals = MIN_Z;
        else if (metals > MAX_Z)
            metals = MAX_Z;
        branch[0].snap = snap;
        branch[0].metals = metals;
        branch[0].sfr = sfr;
    }
    trace_progenitors(snap - 1, pFirstProgenitor[snap][galIdx], branch, &nProg);
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


inline void compute_mags(float *mags, double z, 
                         double *spectrum, double *absorption, 
                         double *waves, double *obsWaves, int nWaves,
                         double *filters, int nRest, int nObs) {
    int i, j;
    double *flux;
    double *fFlux;
    double *pFilter;

    flux = (double*)malloc(nWaves*sizeof(double));
    fFlux = (double*)malloc(nWaves*sizeof(double));

    // Compute rest frame magnitudes
    for(i = 0; i < nWaves; ++i)
        flux[i] = 3.34e4*waves[i]*waves[i]*spectrum[i]/SURFACE_AREA; // Unit Jy
    for(i = 0; i < nRest; ++i) {
        pFilter = filters + i*nWaves;
        for(j = 0; j < nWaves; ++j) {
            fFlux[j] = flux[j]*pFilter[j];
            if (fFlux[j] < TOL)
                fFlux[j] = TOL;
        }
        mags[i] = (float)(-2.5*log10(trapz(fFlux, waves, nWaves)) + 8.9);
    }
    // Compute observed frame magnitudes
    for(i = 0; i < nWaves; ++i)
        // Unit Jy
        flux[i] = 3.34e4*obsWaves[i]*obsWaves[i]*spectrum[i]/SURFACE_AREA/(1. + z)*absorption[i]; 
    for(i = 0; i < nObs; ++i) {
        pFilter = filters + (nRest + i)*nWaves;
        for(j = 0; j < nWaves; ++j) {
            fFlux[j] = flux[j]*pFilter[j];
            if (fFlux[j] < TOL)
                fFlux[j] = TOL;
        }
        mags[nRest + i] = (float)(-2.5*log10(trapz(fFlux, obsWaves, nWaves)) + 8.9);
    }

    free(flux);
    free(fFlux);
}

/*
void galaxy_spectra_cext(double *pOutput, 
                         double z, int snap,
                         int *indices, int nGal,
                         double *ageList, int nAgeList,
                         int **pFP, int **pNP, float **pM, float **pS) {
    int i;
    int nProg;
    struct template rawSpectra;
    struct template spectra;
    struct node branch[MAX_NODE];

    pFirstProgenitor = pFP;
    pNextProgenitor = pNP;
    pMetals = pM;
    pSFR = pS;

    read_sed_templates(&rawSpectra);   
    integrate_sed_templates(ageList, nAgeList, &spectra, &rawSpectra);

    timing_start("# Compute SEDs\n");
    
    for(i = 0; i < nGal; ++i, report(i, nGal)) {
        nProg = trace_merger_tree(snap, indices[i], branch);
        compute_spectrum(pOutput + i*spectra.nWaves, snap, branch, nProg, &spectra);
    }

    timing_end();
    
    free_template(&rawSpectra);
    free_template(&spectra);
}
*/

void galaxy_mags_cext(float *pOutput, 
                      double z, int snap,
                      int *indices, int nGal,
                      double *ageList, int nAgeList,
                      double *filters, int nRest, int nObs,
                      double *absorption) {
    int i;
    int nProg;
    int nWaves;
    int nTotal = nObs + nRest;
    double *spectrum;
    double *obsWaves;
    struct template rawSpectra;
    struct template spectra;
    struct node branch[MAX_NODE];

    read_sed_templates(&rawSpectra);   
    integrate_sed_templates(ageList, nAgeList, &spectra, &rawSpectra);

    timing_start("# Compute magnitudes\n");
    printf("# Snapshot %d\n", snap);
    nWaves = spectra.nWaves;
    spectrum = (double*)malloc(nWaves*sizeof(double));
    obsWaves = (double*)malloc(nWaves*sizeof(double));
    for(i = 0; i < nWaves; ++i)
        obsWaves[i] = (1. + z)*spectra.waves[i];

    for(i = 0; i < nGal; ++i, report(i, nGal)) {
        nProg = trace_merger_tree(snap, indices[i], branch);
        compute_spectrum(spectrum, snap, branch, nProg, &spectra);
        compute_mags(pOutput + i*nTotal, z, 
                     spectrum, absorption, 
                     spectra.waves, obsWaves, nWaves,
                     filters, nRest, nObs);
    }
    timing_end();

    free(spectrum);
    free(obsWaves);
    free_template(&rawSpectra);
    free_template(&spectra);
}



