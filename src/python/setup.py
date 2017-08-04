from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

sourcefiles = ['spence.pyx']

extensions = [Extension("spence", sourcefiles)]

setup(
    ext_modules = cythonize(extensions)
)
