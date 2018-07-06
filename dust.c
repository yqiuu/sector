#define _DUST_
#include"sector_cext.h"

inline int birth_cloud_interval(double tBC, double *ageStep, int nAgeStep) {
    if (tBC >= ageStep[nAgeStep - 1])
        return nAgeStep;
    else if (tBC < ageStep[0])
        return 0;
    else
        return bisection_search(tBC, ageStep, nAgeStep) + 1;
}


void init_templates_special(struct sed_params *spectra, double tBC ) {
    // Special templates are for birth cloud 
    // Dimension: metallicity × wavelength
    int iZW;
    double *pData;
    int nZ = spectra->nZ;
    int nWaves = spectra->nWaves;
    int nAge = spectra->nAge;
    int nZW = nZ*nWaves;
    double *age = spectra->age;
    double *rawData = spectra->raw;

    // Find the time inverval containning the birth cloud
    int nAgeStep = spectra->nAgeStep;
    double *ageStep = spectra->ageStep;
    int iAgeBC = birth_cloud_interval(tBC, ageStep, nAgeStep);
    double t0, t1;

    if (iAgeBC == nAgeStep) {
        spectra->inBC = (double*)calloc(nZ*nWaves, sizeof(double));
        spectra->outBC = (double*)calloc(nZ*nWaves, sizeof(double));
        return;
    }
    else if (iAgeBC == 0) {
        t0 = age[0];
        t1 = ageStep[0];
        if (tBC < t0)
            tBC = t0;
    }
    else {
        t0 = ageStep[iAgeBC - 1];
        t1 = ageStep[iAgeBC];
    }

    spectra->inBC = (double*)malloc(nZ*nWaves*sizeof(double));
    pData = spectra->inBC;
    for(iZW = 0; iZW < nZW; ++iZW) 
        *pData++ = trapz_table(rawData + iZW*nAge, age, nAge, t0, tBC);
    spectra->outBC = (double*)malloc(nZ*nWaves*sizeof(double));
    pData = spectra->outBC;
    for(iZW = 0; iZW < nZW; ++iZW) 
        *pData++ = trapz_table(rawData + iZW*nAge, age, nAge, tBC, t1);
}


void dust_absorption_full(struct sed_params *spectra, struct dust_params *dustParams) {
    /* tBC:   life time of the birth clound
     * nu:    fraction of ISM dust absorption
     * tauUV: UV-band absorption optical depth
     * nBC:   power law index of tauBC
     * nISM:  power law index of tauISM
     *
     * Reference: da Cunha et al. 2008
     */
    int iA, iW, iZ;
    int nZ = spectra->nZ;
    int nWaves = spectra->nWaves;
    double *waves = spectra->waves;
    int nAgeStep = spectra->nAgeStep;
    double *pData = spectra->ready;
    double *pInBC = spectra->inBC;
    double *pOutBC = spectra->outBC;

    double tauUV_ISM = dustParams->tauUV_ISM;
    double nISM = dustParams->nISM;
    double tauUV_BC = dustParams->tauUV_BC;
    double nBC = dustParams->nBC;
    int iAgeBC = birth_cloud_interval(dustParams->tBC, spectra->ageStep, nAgeStep);
    double *transISM = malloc(nWaves*sizeof(double));
    double *transBC = malloc(nWaves*sizeof(double));
    double ratio;

    // Compute the optical depth of both birth cloud and ISM
    for(iW = 0; iW < nWaves; ++iW) {
        ratio = waves[iW]/1600.;
        transISM[iW] = exp(-tauUV_ISM*pow(ratio, nISM));
        transBC[iW] = exp(-tauUV_BC*pow(ratio, nBC));
    }

    for(iZ = 0; iZ < nZ; ++iZ) {
        for(iA = 0; iA < nAgeStep; ++iA) {
            // Apply optical depth of birth cloud
            if (iA < iAgeBC)
                for(iW = 0; iW < nWaves; ++iW)
                    pData[iW] *= transBC[iW];
            else if (iA == iAgeBC)
                for(iW = 0; iW < nWaves; ++iW)
                    pData[iW] = transBC[iW]*pInBC[iW] + pOutBC[iW];
            // Apply optical depth of ISM
            for(iW = 0; iW < nWaves; ++iW)
                pData[iW] *= transISM[iW];
            pData += nWaves;
        }
        pInBC += nWaves;
        pOutBC += nWaves;
    }

    free(transISM);
    free(transBC);
}

