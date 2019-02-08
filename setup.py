from distutils.core import setup, Extension
from glob import glob
from Cython.Build import cythonize

ext = Extension('pyolim', ['src/pyolim.pyx'], libraries=['olim'])

setup(name='pyolim', ext_modules=cythonize([ext]))
