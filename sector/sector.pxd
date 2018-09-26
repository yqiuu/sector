cdef extern from "clib/sector.h":
    ctypedef struct ssp_t:
        short index
        float metals
        float sfr

    ctypedef struct csp_t:
        ssp_t *bursts
        int nBurst

    ctypedef struct gal_params_t:
        double z
        int nAgeStep
        double *ageStep
        int nGal
        int *indices
        csp_t *histories

    ctypedef struct sed_params_t:
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

    ctypedef struct dust_params_t:
        double tauUV_ISM
        double nISM
        double tauUV_BC
        double nBC
        double tBC


    void add_Lyman_absorption(double *target, double *waves, int nWaves, double z)
    void init_templates_raw(sed_params_t *spectra, char *fName)
    void shrink_templates_raw(sed_params_t *spectra, double maxAge)
    void init_filters(sed_params_t *spectra,
                      double *betaBands, int nBeta, double *restBands, int nRest,
                      double *obsTrans, double *obsWaves, int *nObsWaves, int nObs, double z)


cdef extern from "clib/sector.h" nogil:
    double *composite_spectra_cext(sed_params_t *spectra,
                                   gal_params_t *galParams, dust_params_t *dustParams,
                                   short outType, short approx, short nThread)
