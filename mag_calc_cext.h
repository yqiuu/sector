int **g_firstProgenitor;
int **g_nextProgenitor;
float **g_metals;
float **g_sfr;

struct props {
    short index;
    float metals;
    float sfr;
};

struct prop_set {
    struct props *nodes;
    int nNode;
};

struct prop_set *read_properties_by_progenitors(int **firstProgenitor, int **nextProgenitor,
                                                float **galMetals, float **galSFR,
                                                int tSnap, int *indices, int nGal);


struct sed_params {
    double *Z;
    int nZ;
    int minZ;
    int maxZ;
    double *waves;
    int nWaves;
    double *age;
    int nAge;
    double *data;
};


struct dust_params {
    double tauUV_ISM;
    double nISM;
    double tauUV_BC;
    double nBC;
    double tBC;
};


float *composite_spectra_cext(struct sed_params *rawSpectra,
                              struct prop_set *galProps, int nGal,
                              double z, double *ageList, int nAgeList,
                              double *filters, double *logWaves, int nFlux, int nObs,
                              double *absorption, struct dust_params *dustArgs,
                              int outType);
