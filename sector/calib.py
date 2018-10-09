import numpy as np
import pandas as pd

from .dust import compute_mags_mhysa


try:
    from mhysa import Constraint

    class LF_CMR(Constraint):
        def __init__(
            self, name, snapshots, obsData, dustParams, magCut = -15., label = None,
            saveMags = False, prefix = "mags", outPath = "outputs/"):
            #
            super().__init__(name, snapshots, label)
            self.dustParams = dustParams
            self.magCut = magCut
            self.saveMags = saveMags
            self.prefix = prefix
            self.outPath = outPath


        def add_parameter(self, *args, fixed, **kwargs):
            if 'tBC' in args or 'tBC' in kwargs.keys():
                raise KeyError("Varying tBC is not supported!")
            if fixed:
                try:
                    name, val = args
                except:
                    raise Exception("If fixed, pass parameter name and value as two positional arguments!")
                self.dustParams.set_fixed_params(**{name:val})
            else:
                super().add_parameter(*args, **kwargs)


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
                    raise Exception("Snapshots and those in the input file mismath!")
            self.snapDict = snapDict
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
            lnL = 0.
            if sampler.is_controller:
                inBCFlux, outBCFlux, gals = self._squeeze_gals(snapshot, gals)
                # Compute UV properties with dust
                dustParams = self.dustParams(
                    [p.current_value for p in self.params], gals, sampler.meraxes_globals.ZZ[snapshot]
                )
                M1600, beta = compute_mags_mhysa(
                    inBCFlux, outBCFlux, self.centreWaves, self.logWaves, self.nBeta, self.nFlux, dustParams
                )
                #
                if self.saveMags:
                    fName = "%s/%s_snap%03d_%s_%d.hdf5"%(
                        self.outPath, self.prefix, snapshot,
                        '_'.join([str(p.current_value) for p in sampler.iter_all_params()]), sampler.i_comm
                    )
                    self._save_mags(fName, gals['ID'], M1600, beta)
                    print("[i_comm %d] Save output."%sampler.i_comm)
                # compute lnlikelihood
            return lnL


except:
    raise Exception("Cannot import mhysa!")
