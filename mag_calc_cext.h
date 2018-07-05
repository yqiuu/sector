struct sed_params {
    // Raw templates
    int minZ;
    int maxZ;
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
    double *logWaves;
    // IGM absorption
    double *LyAbsorption;
    // Working templates
    int nAgeStep;
    double *ageStep;
    double *integrated;
    double *ready;
    double *working;
    double *inBC;
    double *outBC;
};

void init_templates_raw(struct sed_params *spectra, char *fName);
void shrink_templates_raw(struct sed_params *spectra, double maxAge);


//void init_filters(struct sed_params *spectra, double *filters, int nFlux, int nObs,
//                  double *logWaves, double *LyAbsorption);


//void free_filters(struct sed_params *spectra);


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


double *composite_spectra_cext(struct sed_params *spectra,
                               struct gal_params *galParams, struct dust_params *dustParams,
                               short outType, short nThread);
