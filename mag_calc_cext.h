int **firstProgenitor;
int **nextProgenitor;
float **galMetals;
float **galSFR;

float *galaxy_spectra_cext(double z, int tSnap, 
                           int *indices, int nGal,
                           double *ageList, int nAgeList,
                           int nWaves);
 

float *galaxy_mags_cext(double z, int tSnap,
                        int *indices, int nGal,
                        double *ageList, int nAgeList,
                        double *filters, int nRest, int nObs,
                        double *absorption);
