'''This is a cython setup file for permutation.pyx'''

from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("permute.pyx"))
