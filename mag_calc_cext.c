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
#define JANSKY(x) (3.34e4*(x)*(x)/SURFACE_AREA) // convert erg/s/A to Jy at 10 pc
#define M_AB(x) (-2.5*log10(x) + 8.9) // convert Jansky to AB magnitude
#define TOL 1e-50

time_t sTime;
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


void free_template(struct template *spectra) {
    free(spectra->age);
    free(spectra->waves);
    free_2d_double(spectra->data, spectra->nZ*spectra->nAge);
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


inline void get_integrand(double *integrand, struct template *spectra, int idxZ, int idxW) {
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

    printf("# Interpolate SED templates in terms of metallicity\n");
    spectra->data = malloc_2d_double(spectra->nZ*spectra->nAge, spectra->nWaves);
    for(i = 0; i < spectra->nZ; ++i)
        for(j = 0; j < spectra->nAge; ++j) {
            pData = spectra->data[i*spectra->nAge + j];
            for(k = 0; k < spectra->nWaves; ++k) 
                pData[k] = interp((double)i, refMetals, 
                                  refSpectra[j*spectra->nWaves + k], NMETALS);
        }
    
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
    int nFilter = nRest + nObs;
    double *pFilter;   
    struct template rawSpectra;
    int nWaves;
    double *waves;
    // Spectra after integration over time
    // The first dimension refers to metallicites and ages
    // The last dimension refers to wavelengths
    double **spectra;
    double *integrand;
    // Spectra to be interploated along metallicities
    // The first dimension refers to filters and ages
    // Thw last dimension refers to metallicites
    double **refSpectra;
    double *pData;

    read_sed_templates(&rawSpectra);
    timing_start("# Process SED templates\n");
    nWaves = rawSpectra.nWaves;
    waves = (double*)malloc(nWaves*sizeof(double));
    memcpy(waves, rawSpectra.waves, nWaves*sizeof(double));
    // Integrate raw SED templates over time
    printf("# Integrate SED templates over time\n");
    spectra = malloc_2d_double(NMETALS*nAgeList, nWaves);
    integrand = (double*)malloc(rawSpectra.nAge*sizeof(double));
    for(i = 0; i < NMETALS; report(i++, NMETALS)) 
        for(j = 0; j < nAgeList; ++j) {
            pData = spectra[i*nAgeList + j];
            for(k = 0; k < nWaves; ++k) {
                get_integrand(integrand, &rawSpectra, i, k);
                if (j == 0) {
                    // The first time step of SED templates is typicall not zero
                    // Here assumes that the templates is constant beween zero
                    // and the first time step
                    pData[k] = integrand[0]*rawSpectra.age[0];
                    pData[k] += trapz_table(integrand, rawSpectra.age, rawSpectra.nAge, 
                                            rawSpectra.age[0], ageList[j]);

                }
                else
                    pData[k] = trapz_table(integrand, rawSpectra.age, rawSpectra.nAge, 
                                           ageList[j - 1], ageList[j]);
            }
        }
    free_template(&rawSpectra);
    free(integrand);
    
    // Intgrate SED templates over filters
    refSpectra = malloc_2d_double(nFilter*nAgeList, NMETALS);
    // Compute rest frame flux
    // Unit Jy
    printf("# Compute rest frame flux\n");
    for(i = 0; i < NMETALS*nAgeList; ++i) {
        pData = spectra[i];
        for(k = 0; k < nWaves; ++k)
            pData[k] *= JANSKY(waves[k]);
    }
    for(i = 0; i < nRest; ++i) {
        pFilter = filters + i*nWaves;
        for(j = 0; j < nAgeList; ++j) {
            pData = refSpectra[i*nAgeList + j];
            for(k = 0; k < NMETALS; ++k)
                pData[k] = trapz_filter(pFilter, spectra[k*nAgeList + j], waves, nWaves);
        }
    }
    // Compute observer frame flux
    printf("# Compute observer frame flux\n");
    // Transform everything to observer frame
    // Note the flux in this case is a function frequency
    // Therefore the flux has a factor of 1 + z
    for(i = 0; i < nWaves; ++i)
        waves[i] *= 1. + z;
    for(i = 0; i < NMETALS*nAgeList; ++i) {
        pData = spectra[i];
        for(k = 0; k < nWaves; ++k)
            pData[k] *= (1. + z)*absorption[k];
    }
    for(i = nRest; i < nFilter; ++i) {
        pFilter = filters + i*nWaves;
        for(j = 0; j < nAgeList; ++j) {
            pData = refSpectra[i*nAgeList + j];
            for(k = 0; k < NMETALS; ++k)
                pData[k] = trapz_filter(pFilter, spectra[k*nAgeList + j], waves, nWaves);
        }
    }

    // Interploate SED templates along metallicities
    printf("# Interpolate SED templates along metallicities\n");
    for(i = 0; i < NUM_Z; ++i)
        for(j = 0; j < nAgeList; ++j) {
            pData = fluxTmp[i*nAgeList + j];
            for(k = 0; k < nFilter; ++k) 
                pData[k] = interp((double)i, refMetals, 
                                  refSpectra[k*nAgeList + j], NMETALS);
        }
 
    timing_end();
    
    free_2d_double(spectra, NMETALS*nAgeList);
    free_2d_double(refSpectra, nFilter*nAgeList); 
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

    firstProgenitor = pFP;
    nextProgenitor = pNP;
    galMetals = pM;
    galSFR = pS;

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

void galaxy_mags_cext(float *mags, 
                      double z, int tSnap,
                      int *indices, int nGal,
                      double *ageList, int nAgeList,
                      double *filters, int nRest, int nObs,
                      double *absorption) {
    int i, j, k;

    int nFilter = nObs + nRest;
    double **fluxTmp = malloc_2d_double(NUM_Z*nAgeList, nFilter);
    double *flux = (double*)malloc(nFilter*sizeof(double));
    double *pData;

    int nProg;
    struct node branch[MAX_NODE];
    int snap;
    double sfr;

    init_templates_ph(fluxTmp, z,
                      ageList, nAgeList, 
                      filters, nRest, nObs,
                      absorption);
    
    timing_start("# Compute magnitudes\n");
    for(i = 0; i < nGal; report(i++, nGal)) {
        nProg = trace_merger_tree(tSnap, indices[i], branch);
        // Initialise flux
        for(j = 0; j < nFilter; ++j)
            flux[j] = TOL;
        // Sum over contributions from all progentiors
        for(j = 0; j < nProg; ++j) {
            snap = branch[j].snap;
            sfr = branch[j].sfr;    
            pData = fluxTmp[branch[j].metals*nAgeList + tSnap - snap];
            for(k = 0 ; k < nFilter; ++k) 
                flux[k] += sfr*pData[k];
        }
        // Convert fluxes to magnitudes
        for(j = 0; j < nFilter; ++j)
            *mags++ = (float)M_AB(flux[j]);
    }
    timing_end();
    
    free(flux);
    free_2d_double(fluxTmp, nFilter);
}



