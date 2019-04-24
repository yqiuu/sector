import os


__all__ = [
    'STARBURST99_Salpeter',
    'STARBURST99_Kroupa'
]


packageDir = os.path.dirname(os.path.abspath(__file__))
STARBURST99_Salpeter = os.path.join(packageDir, 'STARBURST99_Salpeter.hdf5')
STARBURST99_Kroupa = os.path.join(packageDir, 'STARBURST99_Kroupa.hdf5')
