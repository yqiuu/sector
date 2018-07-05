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


void dust_absorption(struct sed_params *spectra, struct dust_params *dustParams, int *ageFlag) {
    /* tBC: life time of the birth clound
     * nu: fraction of ISM dust absorption
     * tauUV: V-band absorption optical depth
     * nBC: power law index of tauBC
     * nISM: power law index of tauISM
     *
     * Reference: da Cunha et al. 2008
     */
    int iA, iW, i, n;
    double *pData;

    int nZ = spectra->nZ;
    int nWaves = spectra->nWaves;
    double *waves = spectra->waves;
    int nAgeStep = spectra->nAgeStep;
    double *data = spectra->ready;
    double *inBC = spectra->inBC;
    double *outBC = spectra->outBC;

    double tauUV_ISM = dustParams->tauUV_ISM;
    double nISM = dustParams->nISM;
    double tauUV_BC = dustParams->tauUV_BC;
    double nBC = dustParams->nBC;
    double tBC = dustParams->tBC;
    double *transISM = malloc(nWaves*sizeof(double));
    double *transBC = malloc(nWaves*sizeof(double));
    double ratio;
    int iAgeBC = birth_cloud_interval(tBC, spectra->ageStep, nAgeStep);

    // Compute the optical depth of both the birth cloud and the ISM
    for(iW = 0; iW < nWaves; ++iW) {
        ratio = waves[iW]/1600.;
        transISM[iW] = exp(-tauUV_ISM*pow(ratio, nISM));
        transBC[iW] = exp(-tauUV_BC*pow(ratio, nBC));
    }

    // t_s < tBC < t_s + dt
    if (iAgeBC != nAgeStep && !ageFlag[iAgeBC]) {
        // loop info: n = nZ*nWaves
        //            iZ = i/nWaves
        //            iW = i%nWaves
        n = nZ*nWaves;
        for(i = 0; i < n; ++i) {
            iW = i%nWaves;
            data[(i/nWaves*nAgeStep + iAgeBC)*nWaves + iW] = transBC[iW]*inBC[i] + outBC[i];
        }
    }

    // tBC > t_s
    n = iAgeBC*nZ;
    for(i = 0; i < n; ++i) {
        iA = i%iAgeBC;
        if (ageFlag[iA])
            continue;
        pData = data + (i/iAgeBC*nAgeStep + iA)*nWaves;
        for(iW = 0; iW < nWaves; ++iW)
            pData[iW] *= transBC[iW];
    }

    n = nAgeStep*nZ;
    for(i = 0; i < n; ++i) {
        if (ageFlag[i%nAgeStep])
            continue;
        pData = data + i*nWaves;
        for(iW = 0; iW < nWaves; ++iW)
            pData[iW] *= transISM[iW];
    }

    free(transISM);
    free(transBC);
}

