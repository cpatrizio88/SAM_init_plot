from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
ext_modules = cythonize('plot_CRHsort_correlationsCY.pyx'),
              include_dirs=[numpy.get_include()]
                      
)