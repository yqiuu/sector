from distutils.core import setup, Extension
from Cython.Build import cythonize

exts = cythonize([Extension("mag_calc", sources = ["mag_calc_cext.c", "mag_calc.pyx"])])
setup(ext_modules = exts, )
