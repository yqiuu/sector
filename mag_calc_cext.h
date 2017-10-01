int **firstProgenitor;
int **nextProgenitor;
float **galMetals;
float **galSFR;

float *composite_spectra_cext(double z, int tSnap,
                              int *indices, int nGal,
                              double *ageList, int nAgeList,
                              double *filters, int nRest, int nObs,
                              double *absorption,
                              int dust, double tBC, double mu, 
                              double tauV, double nBC, double nISM);
