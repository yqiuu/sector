import logging
import numpy as np
import pandas as pd

from .dust import compute_mags_mhysa
from copy import deepcopy
from astropy.stats import biweight_location


class likelihood_UV:
    def __init__(self, obsData, volume, keys = None, mincnt = 5, nan = 10e6, blob = False):
        obsData = self._drop_bright_bins(np.atleast_1d(obsData), volume, mincnt)
        if keys is None:
            self.obsData = obsData
            self.keys = np.arange(len(obsData))
        else:
            keys = np.atleast_1d(keys)
            if len(keys) != len(obsData):
                raise ValueError("Number of keys should be equal to number of observations!")
            self.obsData = {k:d for (k, d) in zip(keys, obsData)}
            self.keys = keys
        self.volume = volume
        self.nan = nan
        self.blob = blob


    def _drop_bright_bins(self, obsData, volume, mincnt):
        obsData = deepcopy(obsData)
        for d in obsData:
            cond = np.full(len(d['obsLF']), True)
            for iB, ((lower, upper), lf) in enumerate(zip(d['binLF'], d['obsLF'])):
                if lf*(upper - lower)*volume < mincnt:
                    cond[iB] = False
            #
            d['binLF'] = d['binLF'][cond]
            d['obsLF'] = d['obsLF'][cond]
            d['obsLFerr'] = d['obsLFerr'][cond]
            # Require that the lower bound of the brighest bin of CMRs be fainter than that of LFs
            # Assume all bins are sorted, and the first bin is brightest.
            cond = d['binCMR'][:, 0] >= d['binLF'][0, 0]
            d['binCMR'] = d['binCMR'][cond]
            d['obsCMR'] = d['obsCMR'][cond]
            d['obsCMRerr'] = d['obsCMRerr'][cond]
        return obsData


    def _fix_dim(self, keys, M1600, beta):
        if keys is None:
            keys = self.keys
        else:
            keys = np.atleast_1d(keys)
        if len(keys) == 1:
            M1600 = np.atleast_2d(M1600)
            if beta is not None:
                beta = np.atleast_2d(beta)
        return keys, M1600, beta


    @staticmethod
    def _compute_lnL(model, obs, obsErr):
        delta = (model - obs)/obsErr
        return -.5*np.sum(delta*delta + np.log(2*np.pi*obsErr*obsErr))


    def _retval(self, lnL, blob):
        if np.isnan(lnL):
            lnL = self.nan
        if self.blob:
            return lnL, blob
        else:
            return lnL


    def eval_LF(self, M1600, keys = None):
        lnL = 0.
        blob = {}
        keys, M1600, _ = self._fix_dim(keys, M1600, None)
        volume = self.volume
        for k, mags in zip(keys, M1600):
            obsLF = self.obsData[k]['obsLF']
            obsLFerr = self.obsData[k]['obsLFerr']
            modelLF = np.zeros(len(obsLF))
            for iB, (lower, upper) in enumerate(self.obsData[k]['binLF']):
                width = upper - lower
                number = np.sum((mags >= lower) & (mags < upper))
                if number == 0:
                    modelLF[iB] = np.nan
                else:
                    modelLF[iB] = float(number)/width/volume
            lnL += self._compute_lnL(modelLF, obsLF, obsLFerr)
            blob[k] = modelLF
        return self._retval(lnL, blob)


    def eval_CMR(self, M1600, beta, keys = None):
        lnL = 0.
        blob = {}
        keys, M1600, beta = self._fix_dim(keys, M1600, beta)
        for k, mags, b in zip(keys, M1600, beta):
            obsCMR = self.obsData[k]['obsCMR']
            obsCMRerr = self.obsData[k]['obsCMRerr']
            modelCMR = np.zeros(len(obsCMR))
            for iB, (lower, upper) in enumerate(self.obsData[k]['binCMR']):
                sample = b[(mags >= lower) & (mags < upper)]
                if len(sample) < 2:
                    modelCMR[iB] = np.nan
                else:
                    modelCMR[iB] = biweight_location(sample)
            lnL += self._compute_lnL(modelCMR, obsCMR, obsCMRerr)
            blob[k] = modelCMR
        return self._retval(lnL, blob)


    def __call__(self, M1600, beta, keys = None):
        if self.blob:
            val1, blob1 = self.eval_LF(M1600, keys)
            val2, blob2 = self.eval_CMR(M1600, beta, keys)
            return val1 + val2, (blob1, blob2)
        else:
            return self.eval_LF(M1600, keys) + self.eval_CMR(M1600, beta, keys)


try:
    from mhysa import Constraint

    class LF_CMR(Constraint):
        def __init__(
            self, name, snapshots, obsData, dustParams, magCut = -15., label = None,
            saveMags = False, prefix = "mags", outPath = "outputs/", **lnKwargs):
            #
            super().__init__(name, snapshots, label)
            self.obsData = obsData
            self.dustParams = dustParams
            self.magCut = magCut
            self.saveMags = saveMags
            self.prefix = prefix
            self.outPath = outPath
            self.lnKwargs = lnKwargs


        def add_parameter(self, name, p0, *args, **kwargs):
            if name == 'tBC':
                raise ValueError("Varying tBC is not supported!")
            if name in self.dustParams.fixed_names():
                raise ValueError("Dust parameter %s is fixed!"%name)
            super().add_parameter(name, p0, *args, **kwargs)


        def controller_setup(self, sampler):
            # Get magnitude parameters
            snapshots, nBeta, nRest, tBC, centreWaves, logWaves = sampler.meraxes_globals.mag_params()
            # Set lifetime of birth cloud
            dustParams = self.dustParams
            dustParams.set_fixed_params(tBC = tBC)
            #
            if np.sum(np.isnan(dustParams.fixedParams)) != len(self.params):
                raise ValueError("All parameters in %s should be added!"%dustParams.paramNames.__str__())
            # Construct the map between snapshots and indices of flux arrays
            snapshots = list(snapshots)
            snapDict = {}
            for snap in self.required_snapshots:
                try:
                    snapDict[snap] = snapshots.index(snap)
                except:
                    raise Exception("Snapshots and those in the input file mismatch!")
            self.snapDict = snapDict
            # Set likelihood function
            self.lnKwargs['blob'] = True
            self.estimator = likelihood_UV(
                self.obsData, sampler.meraxes_globals.comoving_volume(), keys = self.required_snapshots, **self.lnKwargs
            )
            # Check if the UV band is right
            if not np.isclose(centreWaves[nBeta], 1600.):
                raise ValueError("The centre wavelength of the first rest-frame filters should be 1600 angstrom!")
            # Set other attributes
            self.nBeta = nBeta
            self.nFlux = sampler.meraxes_globals.MAGS_N_BANDS
            self.centreWaves = centreWaves
            self.logWaves = logWaves
            # Get SFR unit
            self.sfrUnit = sampler.meraxes_globals.UnitSfr_in_solar_mass_per_year
            # Convert magnitude cut into flux cut in meraxes intrinsic unit`
            self.fluxCut = 10**(-.4*(self.magCut - 8.9))/self.sfrUnit


        def _squeeze_gals(self, snapshot, gals):
            """Apply magnitude cut, and convert flux unit to Jansky, SFR unit to M_solar/yr"""
            iS = self.snapDict[snapshot]
            nFlux = self.nFlux
            iF1600 = self.nBeta
            gals = gals[gals['inBCFlux'][:, iF1600] + gals['outBCFlux'][:, iF1600] > self.fluxCut]
            inBCFlux = gals['inBCFlux'].flatten()
            outBCFlux = gals['outBCFlux'].flatten()
            #
            sfrUnit = self.sfrUnit
            inBCFlux *= sfrUnit
            outBCFlux *= sfrUnit
            gals['Sfr'] *= sfrUnit
            return inBCFlux, outBCFlux, gals


        @staticmethod
        def _save_mags(fName, ID, M1600, beta):
            df = pd.DataFrame()
            df['ID'] = ID
            df['M1600'] = M1600
            df['beta'] = beta
            df.to_hdf(fName, "w")


        def lnlikelihood(self, snapshot, gals, sampler):
            if sampler.is_controller:
                inBCFlux, outBCFlux, gals = self._squeeze_gals(snapshot, gals)
                # Compute UV properties with dust
                dustParams = self.dustParams(
                    [p.current_value for p in self.params], gals, sampler.meraxes_globals.ZZ[snapshot]
                )
                M1600, beta = compute_mags_mhysa(
                    inBCFlux, outBCFlux, self.centreWaves, self.logWaves, self.nBeta, self.nFlux, dustParams
                )
                # Compute lnlikelihood
                lnL, self.blob = self.estimator(M1600, beta, keys = snapshot)
                logging.debug('LF & CMR z=%.2f lnL=%.2e'%(self.redshift, lnL))
                #
                if self.saveMags:
                    fName = "%s/%s_snap%03d_%s_%d.hdf5"%(
                        self.outPath, self.prefix, snapshot,
                        '_'.join([str(p.current_value) for p in sampler.iter_all_params()]), sampler.i_comm
                    )
                    self._save_mags(fName, gals['ID'], M1600, beta)
                    logging.info("[i_comm %d] Save output."%sampler.i_comm)
                #
                return lnL
            else:
                return 0.


except:
    raise Exception("Cannot import mhysa!")
