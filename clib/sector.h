#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

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
#ifdef CALC_MAGS
typedef struct mini_sed_params_t {
    int iS;
    int targetSnap[MAGS_N_SNAPS];
    int nBeta;
    int minZ;
    int maxZ;
    int nMaxZ;
    int iAgeBC[MAGS_N_SNAPS];
    size_t totalSize;
    double *working;
    double *inBC;
    double *outBC;
    double *centreWaves;
    double *logWaves;
} mini_sed_params_t;
#endif

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
#ifdef CALC_MAGS
void init_templates_mini(mini_sed_params_t *miniSpectra, char *fName,
                         double *LTTime, int *targetSnaps, double *redshifts,
                         double *restBands, int nRest, int nBeta, double tBC);
#endif

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

#ifndef _SPECTRA_
void init_templates_raw(struct sed_params *spectra, char *fName);
void shrink_templates_raw(struct sed_params *spectra, double maxAge);
void init_filters(struct sed_params *spectra,
                  double *betaBands, int nBeta, double *restBands, int nRest,
                  double *obsTrans, double *obsWaves, int *nObsWaves, int nObs, double z);
void init_templates_integrated(struct sed_params *spectra);
void init_templates_working(struct sed_params *spectra, struct csp *pHistories,
                                   struct dust_params *dustParams, int iG);
#endif


#ifndef _DUST_
int birth_cloud_interval(double tBC, double *ageStep, int nAgeStep);
void init_templates_special(struct sed_params *spectra, double tBC, int approx);
void dust_absorption(struct sed_params *spectra, struct dust_params *dustParams);
void dust_absorption_approx(double *inBCFlux, double *outBCFlux, double *centreWaves, int nFlux,
                            struct dust_params *dustParams);
#endif

#ifndef _IGM_
void add_Lyman_absorption(double *target, double *waves, int nWaves, double z);
void init_IGM_absorption(struct sed_params *spectra);
#endif
