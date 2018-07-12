#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Basic functions                                                             *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
FILE *open_file(char *fName, char *mode);

double **malloc_2d_double(int nRow, int nCol);
double **memcpy_2d_double(double **source, int nRow, int nCol);
void free_2d_double(double **p, int nRow);

int bisection_search(double a, double *x, int nX);
double interp(double xp, double *x, double *y, int nPts);
double trapz_table(double *y, double *x, int nPts, double a, double b);
double trapz_filter(double *filter, double *filterWaves, int nFilterwaves,
                    double *flux, double *waves, int nWaves);

struct linResult {
    double slope;
    double intercept;
    double R;
};
struct linResult linregress(double *x, double *y, int nPts);


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * SFHs related                                                                *
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
 * SEDs and dust related                                                       *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
struct sed_params {
    // Raw templates
    int minZ;
    int maxZ;
    int nMaxZ;
    int nZ;
    double *Z;
    int nWaves;
    double *waves;
    int nAge;
    double *age;
    double *raw;
    // Redshift
    double z;
    // Filters
    int nFlux;
    int nObs;
    int *nFilterWaves;
    double *filterWaves;
    double *filters;
    double *centreWaves;
    double *logWaves;
    // IGM absorption
    int igm;
    // Working templates
    int nAgeStep;
    double *ageStep;
    double *integrated;
    double *ready;
    double *working;
    double *inBC;
    double *outBC;
};

struct dust_params {
    double tauUV_ISM;
    double nISM;
    double tauUV_BC;
    double nBC;
    double tBC;
};

#ifndef _SECTOR_
void init_templates_raw(struct sed_params *spectra, char *fName);
void init_filters(struct sed_params *spectra,
                  double *betaBands, int nBeta, double *restBands, int nRest,
                  double *obsTrans, double *obsWaves, int *nObsWaves, int nObs, double z);
void shrink_templates_raw(struct sed_params *spectra, double maxAge);
void init_templates_integrated(struct sed_params *spectra);
void init_templates_working(struct sed_params *spectra, struct csp *pHistories,
                            struct dust_params *dustParams, int iG);
double *composite_spectra_cext(struct sed_params *spectra,
                               struct gal_params *galParams, struct dust_params *dustParams,
                               short outType, short approx, short nThread);
#endif

#ifndef _DUST_
int birth_cloud_interval(double tBC, double *ageStep, int nAgeStep);
void init_templates_special(struct sed_params *spectra, double tBC, int approx);
void dust_absorption(struct sed_params *spectra, struct dust_params *dustParams);
void dust_absorption_approx(double *inBCFlux, double *outBCFlux,
                            struct sed_params *spectra, struct dust_params *dustParams);
#endif

#ifndef _IGM_
void add_Lyman_absorption(double *target, double *waves, int nWaves, double z);
void init_IGM_absorption(struct sed_params *spectra);
#endif
