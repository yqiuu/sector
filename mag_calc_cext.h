struct sed_params {
    int minZ;
    int maxZ;
    int nZ;
    double *Z;
    int nWaves;
    double *waves;
    int nAge;
    double *age;
    int nAgeStep;
    double *ageStep;
    double *raw;
    double *integrated;
    double *ready;
    double *working;
};


struct dust_params {
    double tauUV_ISM;
    double nISM;
    double tauUV_BC;
    double nBC;
    double tBC;
};


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


float *composite_spectra_cext(struct sed_params *rawSpectra,
                              struct gal_params *galParams,
                              double *filters, double *logWaves, int nFlux, int nObs,
                              double *absorption, struct dust_params *dustArgs,
                              short outType, short nThread);
