import numpy as np, pandas as pd
from dragons import meraxes
import magcalc as mc

# Set cosmology
h = 0.678
Om0 = 0.308

# Path to Meraxes output
fname = "/lustre/projects/p102_astro/smutch/meraxes/paper_runs/512/fiducial/output/meraxes.hdf5"

# Read galaxy properties
snapshot = 20
z = meraxes.io.grab_redshift(fname, snapshot)
gals = meraxes.io.read_gals(fname, snapshot, props = ["StellarMass", "GhostFlag"], h = h)
indices = np.where((gals["StellarMass"] > 1e-3) & (gals["GhostFlag"] == False))[0]

# Set filters:
restBands = [[1600., 100.], [2000, 100.], [9000., 200.]]
obsBands = mc.HST_filters(["B435", "V606", "i775", "I814"])

# Path to SED templates
sedPath = "/lustre/projects/p113_astro/yqiu/magcalc/input/STARBURST99-Salpeter-default"
# Compute magnitudes
mags = mc.composite_spectra(fname, snapList = snapshot, gals = indices, 
                            h = 0.678, Om0 = .308, 
                            sedPath = sedPath,
                            outType = "ph",
                            restBands = restBands,
                            obsBands = obsBands,
                            prefix = 'demo')
print mags.iloc[:5, :]

# Dust extinction to UV
MUV = mags.loc[:, ["M1600-100", "M2000-100"]]
AUV = mc.reddening([1600, 2000], MUV["M1600-100"], z)
MUV_dust = MUV + AUV


# Save SFH on the disk
mc.save_star_formation_history(fname, snapList = snapshot, idxList = indices, h = h)
# Compute magnitudes by the stored SFH
mags2 = mc.composite_spectra(fname, snapList = snapshot, gals = 'sfh_%03d.bin'%snapshot, 
                             h = 0.678, Om0 = .308, 
                             sedPath = sedPath,
                             outType = "ph",
                             restBands = restBands,
                             obsBands = obsBands,
                             prefix = 'demo')

assert mags2.equals(mags)

# Check outputs in the current directory
