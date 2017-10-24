int **g_firstProgenitor;
int **g_nextProgenitor;
float **g_metals;
float **g_sfr;

struct props {
    short index;
    short metals;
    float sfr;
};
struct prop_set {
    struct props *nodes;
    int nNode;
};

struct prop_set *read_properties_by_progenitors(int **firstProgenitor, int **nextProgenitor,
                                                float **galMetals, float **galSFR,
                                                int tSnap, int *indices, int nGal);

void free_raw_spectra(void);


void free_int_spectra(void);


float *composite_spectra_cext(struct prop_set *galProps, int nGal,
                              double z, double *ageList, int nAgeList,
                              double *filters, int nRest, int nObs, int mAB,
                              double *absorption,
                              int dust, double tBC, double mu, 
                              double tauV, double nBC, double nISM);

float *UV_slope_cext(struct prop_set *galProps, int nGal,
                     double z, double *ageList, int nAgeList,
                     double *logWaves, double *filters, int nFilter,
                     int dust, double tBC, double mu, 
                     double tauV, double nBC, double nISM);
