cdef extern from "clib/sector.h":
    struct ssp:
        short index
        float metals
        float sfr

    struct csp:
        ssp *bursts
        int nBurst

    struct gal_params:
        double z
        int nAgeStep
        double *ageStep
        int nGal
        int *indices
        csp *histories

    struct sed_params:
        # Raw templates
        int minZ
        int maxZ
        int nMaxZ
        int nZ
        double *Z
        int nWaves
        double *waves
        int nAge
        double *age
        double *raw
        # Redshift
        double z
        # Filters
        int nFlux
        int nObs
        int *nFilterWaves
        double *filterWaves
        double *filters
        double *centreWaves
        double *logWaves
        # IGM absorption
        int igm
        # Working templates
        int nAgeStep
        double *ageStep
        double *integrated
        double *ready
        double *working
        double *inBC
        double *outBC

    struct dust_params:
        double tauUV_ISM
        double nISM
        double tauUV_BC
        double nBC
        double tBC


    void add_Lyman_absorption(double *target, double *waves, int nWaves, double z)
    void init_templates_raw(sed_params *spectra, char *fName)
    void shrink_templates_raw(sed_params *spectra, double maxAge)
    void init_filters(sed_params *spectra,
                      double *betaBands, int nBeta, double *restBands, int nRest,
                      double *obsTrans, double *obsWaves, int *nObsWaves, int nObs, double z)


cdef extern from "clib/sector.h" nogil:
    double *composite_spectra_cext(sed_params *spectra,
                                   gal_params *galParams, dust_params *dustParams,
                                   short outType, short approx, short nThread)
