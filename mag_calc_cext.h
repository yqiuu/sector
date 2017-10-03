int **firstProgenitor;
int **nextProgenitor;
float **galMetals;
float **galSFR;

float *composite_spectra_cext(double z, int tSnap,
                              int *indices, int nGal,
                              double *ageList, int nAgeList,
                              double *filters, int nRest, int nObs, int mAB,
                              double *absorption,
                              int dust, double tBC, double mu, 
                              double tauV, double nBC, double nISM);

float *UV_slope_cext(double z, int tSnap,
                     int *indices, int nGal,
                     double *ageList, int nAgeList,
                     double *logWaves, double *filters, int nFilter,
                     int dust, double tBC, double mu, 
                     double tauV, double nBC, double nISM);
