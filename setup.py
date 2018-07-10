from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

#exts = cythonize([Extension("mag_calc", sources = ["mag_calc_cext.c", "mag_calc.pyx"])])

exts = [Extension("magcalc", 
                  sources = ["tools.c", "sector_cext.c", "dust.c", "igm.c", "magcalc.pyx"],
                  extra_compile_args=['-fopenmp', '-lhdf5', '-lhdf5_hl'], 
                  extra_link_args=['-fopenmp', '-lhdf5', '-lhdf5_hl'])]

for e in exts:
    e.cython_directives = {"embedsignature":True}

setup(name = 'magcalc',
      version = '0.3.0',
      description = 'Package to compute SEDs for semi-analytic models',
      author = 'Yisheng Qiu',
      author_email = 'yishengq@student.unimelb.edu.au',
      ext_modules = exts,
      cmdclass = {'build_ext': build_ext},
      packages = ['filters'],
      package_data = {'filters':['*.npy']})
