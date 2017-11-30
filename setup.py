from distutils.core import setup, Extension
from Cython.Build import cythonize

exts = cythonize([Extension("mag_calc", sources = ["mag_calc_cext.c", "mag_calc.pyx"])])

#exts = cythonize([Extension("mag_calc", sources = ["mag_calc_cext.c", "mag_calc.pyx"],
#                            extra_compile_args=['-fopenmp'], 
#                            extra_link_args=['-fopenmp'])])

setup(ext_modules = exts, )
