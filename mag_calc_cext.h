void galaxy_spectra_cext(double *pOutput, 
                         double z, int snap,
                         int *indices, int nGal,
                         double *ageList, int nAgeList,
                         int **pFP, int **pNP, float **pM, float **pS);

void galaxy_mags_cext(float *pOutput, 
                      double z, int snap,
                      int *indices, int nGal,
                      double *ageList, int nAgeList,
                      double *filters, int nRest, int nObs,
                      double *absorption,
                      int **pFP, int **pNP, float **pM, float **pS);


