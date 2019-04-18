import os

import numpy as np

__all__ = [
    'new_obs_filter',
    'HST_filters',
    'beta_filters'
]


packageDir = os.path.dirname(os.path.abspath(__file__))
filterDict = {
    "B435":os.path.join(packageDir, "HST_ACS_F435W.npy"),
    "V606":os.path.join(packageDir, "HST_ACS_F606W.npy"),
    "i775":os.path.join(packageDir, "HST_ACS_F775W.npy"),
    "I814":os.path.join(packageDir, "HST_ACS_F814W.npy"),
    "z850":os.path.join(packageDir, "HST_ACS_F850LP.npy"),
    "Y098":os.path.join(packageDir, "HST_IR_F098M.npy"),
    "Y105":os.path.join(packageDir, "HST_IR_F105W.npy"),
    "J125":os.path.join(packageDir, "HST_IR_F125W.npy"),
    "H160":os.path.join(packageDir, "HST_IR_F160W.npy"),
    "3.6":os.path.join(packageDir, "HST_IRAC_3.6.npy")
}


def new_obs_filter(name, waves, trans):
    return {'name':name, 'waves':waves, 'trans':trans}


def HST_filters(filterNames):
    """
    Quick access HST filters.

    Parameters
    ----------
    filterNames: list
        Available filters: B435, V606, i775, I814, z850, Y098, Y105,
        J125, H160, 3.6.

    Returns
    -------
    obsBands: list
        For each row, the first element is the filter name, and the
        second element is the transmission curve. The output can be
        passed to ``composite_spectra``.
    """
    obsBands = []
    for name in np.ravel(filterNames):
        waves, trans = np.load(filterDict[name])
        obsBands.append(new_obs_filter(name, waves, trans))
    return obsBands


def beta_filters():
    #=====================================================================
    # return the filters defined by Calzetti et al. 1994, which is used to
    # calculate the UV continuum slope
    #=====================================================================
    windows = np.array([[1268., 1284.],
                        [1309., 1316.],
                        [1342., 1371.],
                        [1407., 1515.],
                        [1562., 1583.],
                        [1677., 1740.],
                        [1760., 1833.],
                        [1866., 1890.],
                        [1930., 1950.],
                        [2400., 2580.]])
    return windows

