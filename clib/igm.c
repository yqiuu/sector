#define _IGM_
#include"sector.h"

#define NLYMAN 39 

void add_Lyman_absorption(double *target, double *waves, int nWaves, double z) {
    double LymanSeries[NLYMAN] = {1215.67, 1025.72, 972.537, 949.743, 937.803,
                                  930.748, 926.226, 923.150, 920.963, 919.352,
                                  918.129, 917.181, 916.429, 915.824, 915.329,
                                  914.919, 914.576, 914.286, 914.039, 913.826,
                                  913.641, 913.480, 913.339, 913.215, 913.104,
                                  913.006, 912.918, 912.839, 912.768, 912.703,
                                  912.645, 912.592, 912.543, 912.499, 912.458,
                                  912.420, 912.385, 912.353, 912.324};
    double LAF1[NLYMAN] = {1.690e-02, 4.692e-03, 2.239e-03, 1.319e-03, 8.707e-04,
                           6.178e-04, 4.609e-04, 3.569e-04, 2.843e-04, 2.318e-04,
                           1.923e-04, 1.622e-04, 1.385e-04, 1.196e-04, 1.043e-04,
                           9.174e-05, 8.128e-05, 7.251e-05, 6.505e-05, 5.868e-05,
                           5.319e-05, 4.843e-05, 4.427e-05, 4.063e-05, 3.738e-05,
                           3.454e-05, 3.199e-05, 2.971e-05, 2.766e-05, 2.582e-05,
                           2.415e-05, 2.263e-05, 2.126e-05, 2.000e-05, 1.885e-05,
                           1.779e-05, 1.682e-05, 1.593e-05, 1.510e-05};
    double LAF2[NLYMAN] = {2.354e-03, 6.536e-04, 3.119e-04, 1.837e-04, 1.213e-04,
                           8.606e-05, 6.421e-05, 4.971e-05, 3.960e-05, 3.229e-05,
                           2.679e-05, 2.259e-05, 1.929e-05, 1.666e-05, 1.453e-05,
                           1.278e-05, 1.132e-05, 1.010e-05, 9.062e-06, 8.174e-06,
                           7.409e-06, 6.746e-06, 6.167e-06, 5.660e-06, 5.207e-06,
                           4.811e-06, 4.456e-06, 4.139e-06, 3.853e-06, 3.596e-06,
                           3.364e-06, 3.153e-06, 2.961e-06, 2.785e-06, 2.625e-06,
                           2.479e-06, 2.343e-06, 2.219e-06, 2.103e-06};
    double LAF3[NLYMAN] = {1.026e-04, 2.849e-05, 1.360e-05, 8.010e-06, 5.287e-06,
                           3.752e-06, 2.799e-06, 2.167e-06, 1.726e-06, 1.407e-06,
                           1.168e-06, 9.847e-07, 8.410e-07, 7.263e-07, 6.334e-07,
                           5.571e-07, 4.936e-07, 4.403e-07, 3.950e-07, 3.563e-07,
                           3.230e-07, 2.941e-07, 2.689e-07, 2.467e-07, 2.270e-07,
                           2.097e-07, 1.943e-07, 1.804e-07, 1.680e-07, 1.568e-07,
                           1.466e-07, 1.375e-07, 1.291e-07, 1.214e-07, 1.145e-07,
                           1.080e-07, 1.022e-07, 9.673e-08, 9.169e-08};
    double DLA1[NLYMAN] = {1.617e-04, 1.545e-04, 1.498e-04, 1.460e-04, 1.429e-04,
                           1.402e-04, 1.377e-04, 1.355e-04, 1.335e-04, 1.316e-04,
                           1.298e-04, 1.281e-04, 1.265e-04, 1.250e-04, 1.236e-04,
                           1.222e-04, 1.209e-04, 1.197e-04, 1.185e-04, 1.173e-04,
                           1.162e-04, 1.151e-04, 1.140e-04, 1.130e-04, 1.120e-04,
                           1.110e-04, 1.101e-04, 1.091e-04, 1.082e-04, 1.073e-04,
                           1.065e-04, 1.056e-04, 1.048e-04, 1.040e-04, 1.032e-04,
                           1.024e-04, 1.017e-04, 1.009e-04, 1.002e-04};
    double DLA2[NLYMAN] = {5.390e-05, 5.151e-05, 4.992e-05, 4.868e-05, 4.763e-05,
                           4.672e-05, 4.590e-05, 4.516e-05, 4.448e-05, 4.385e-05,
                           4.326e-05, 4.271e-05, 4.218e-05, 4.168e-05, 4.120e-05,
                           4.075e-05, 4.031e-05, 3.989e-05, 3.949e-05, 3.910e-05,
                           3.872e-05, 3.836e-05, 3.800e-05, 3.766e-05, 3.732e-05,
                           3.700e-05, 3.668e-05, 3.637e-05, 3.607e-05, 3.578e-05,
                           3.549e-05, 3.521e-05, 3.493e-05, 3.466e-05, 3.440e-05,
                           3.414e-05, 3.389e-05, 3.364e-05, 3.339e-05};

    int iW, iL;
    double lam, tau, ratio;

    for(iW = 0; iW < nWaves; ++iW) {
        tau = 0.;
        lam = waves[iW];
        // Lyman series
        for(iL = 0; iL < NLYMAN; ++iL) {
            ratio = lam/LymanSeries[iL];
            if (ratio < 1. + z) {
                // LAF terms
                if (ratio < 2.2)
                    tau += LAF1[iL]*pow(ratio, 1.2);
                else if (ratio < 5.7)
                    tau += LAF2[iL]*pow(ratio, 3.7);
                else
                    tau += LAF3[iL]*pow(ratio, 5.5);
                /// DLA terms
                if (ratio < 3.)
                    tau += DLA1[iL]*pow(ratio, 2.);
                else
                    tau += DLA2[iL]*pow(ratio, 3.);
            }
        }
        // Lyman continuum
        ratio = lam/912.;
        // LAF terms
        if (z < 1.2) {
            if (ratio < 1. + z)
                tau += .325*(pow(ratio, 1.2) - pow(1. + z, -.9)*pow(ratio, 2.1));
        }
        else if (z < 4.7) {
            if (ratio < 2.2)
                tau += 2.55e-2*pow(1. + z, 1.6)*pow(ratio, 2.1) + .325*pow(ratio, 1.2) \
                       - .25*pow(ratio, 2.1);
            else if (ratio < 1. + z)
                tau += 2.55e-2*(pow(1. + z, 1.6)*pow(ratio, 2.1) - pow(ratio, 3.7));
        }
        else {
            if (ratio < 2.2)
                tau += 5.22e-4*pow(1. + z, 3.4)*pow(ratio, 2.1) + .325*pow(ratio, 1.2) \
                       - 3.14e-2*pow(ratio, 2.1);
            else if (ratio < 5.7)
                tau += 5.22e-4*pow(1. + z, 3.4)*pow(ratio, 2.1) + .218*pow(ratio, 2.1) \
                       - 2.55e-2*pow(ratio, 3.7);
            else if (ratio < 1. + z)
                tau += 5.22e-4*(pow(1. + z, 3.4)*pow(ratio, 2.1) - pow(ratio, 5.5));
        }
        // DLA terms
        if (z < 2.) {
            if (ratio < 1. + z)
                tau += .211*pow(1. + z, 2.) - 7.66e-2*pow(1. + z, 2.3)*pow(ratio, -.3) \
                       - .135*pow(ratio, 2.);
        }
        else {
            if (ratio < 3.)
                tau += .634 + 4.7e-2*pow(1. + z, 3.) - 1.78e-2*pow(1. + z, 3.3)*pow(ratio, -.3) \
                       -.135*pow(ratio, 2.) - .291*pow(ratio, -.3);
            else if (ratio < 1. + z)
                tau += 4.7e-2*pow(1. + z, 3.) - 1.78e-2*pow(1. + z, 3.3)*pow(ratio, -.3) \
                       -2.92e-2*pow(ratio, 3.);
        }
        target[iW] *= exp(-tau);
    }
}


void add_IGM_absorption_filters(struct sed_params *spectra) {
    if (spectra->igm > 0 && spectra->nObs > 0) {
        // Initialise observed wavelengths
        int iF, iW;
        double z = spectra->z;
        int nFlux = spectra->nFlux;
        int nRest = nFlux - spectra->nObs;
        int *nFilterWaves = spectra->nFilterWaves;
        double *filterWaves = spectra->filterWaves;
        int offset = 0;
        int nWaves = 0;
        double *obsWaves;

        for(iF = 0; iF < nFlux; ++iF) {
            if (iF < nRest)
                offset += nFilterWaves[iF];
            else
                nWaves += nFilterWaves[iF];
        }
        obsWaves = (double*)malloc(nWaves*sizeof(double));
        for(iW = 0; iW < nWaves; ++iW)
            obsWaves[iW] = filterWaves[offset + iW]*(1. + z);
        // Add Lyman absorption to filters
        add_Lyman_absorption(spectra->filters + offset, obsWaves, nWaves, z);

        free(obsWaves);
    }
}


void add_IGM_absorption_spectra(struct sed_params *spectra, double *pData, int nGal) {
    if (spectra->igm > 0 && spectra->nObs > 0 && spectra->filters == NULL) {
        // Compute the transmission of Lyman absorption
        int iW;
        int nWaves = spectra->nWaves;
        double *waves = spectra->waves;
        double *obsWaves = malloc(nWaves*sizeof(double));
        double *trans = malloc(nWaves*sizeof(double));
        double z = spectra->z;

        for(iW = 0; iW < nWaves; ++iW) {
            obsWaves[iW] = waves[iW]*(1. + z);
            trans[iW] = 1.;
        }
        add_Lyman_absorption(trans, obsWaves, nWaves, z);
        // Apply the transmission to galaxies
        int iG;
        for(iG = 0; iG < nGal; ++iG) {
            for(iW = 0; iW < nWaves; ++iW)
                pData[iW] *= trans[iW];
            pData += nWaves;
        }

        free(obsWaves);
        free(trans);
    }
}
