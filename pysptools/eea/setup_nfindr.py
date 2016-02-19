from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy


setup(
  name = 'NFINDR app',
  cmdclass = {'build_ext': build_ext},
#  ext_modules = [Extension("_nfindr", ["_nfindr.pyx"],include_dirs=[numpy.get_include()],extra_compile_args=['/openmp',])]
#  ext_modules = [Extension("_nfindr", ["_nfindr.pyx","LU.c"],include_dirs=[numpy.get_include()])]
  ext_modules = [Extension("nfindr", ["nfindr.pyx"],include_dirs=[numpy.get_include()])]
)

