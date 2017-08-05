int **firstProgenitor;
int **nextProgenitor;
float **galMetals;
float **galSFR;

/*
void galaxy_spectra_cext(double *pOutput, 
                         double z, int snap,
                         int *indices, int nGal,
                         double *ageList, int nAgeList,
                         int **pFP, int **pNP, float **pM, float **pS);
*/

void galaxy_mags_cext(float *mags, 
                      double z, int tSnap,
                      int *indices, int nGal,
                      double *ageList, int nAgeList,
                      double *filters, int nRest, int nObs,
                      double *absorption);
