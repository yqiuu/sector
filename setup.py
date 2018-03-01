from distutils.core import setup, Extension
from Cython.Build import cythonize

#exts = cythonize([Extension("mag_calc", sources = ["mag_calc_cext.c", "mag_calc.pyx"])])

exts = cythonize([Extension("magcalc", 
                            sources = ["mag_calc_cext.c", "magcalc.pyx"],
                            extra_compile_args=['-fopenmp'], 
                            extra_link_args=['-fopenmp'])])

setup(name = 'magcalc',
      version = '0.2.4',
      description = 'Package to compute SEDs for a semi-analytic model',
      author = 'Yisheng Qiu',
      author_email = 'yishengq@student.unimelb.edu.au',
      license = 'MIT',
      ext_modules = exts,
      packages = ['filters'],
      package_data = {'filters':['*.npy']})
