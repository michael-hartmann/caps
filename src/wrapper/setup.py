from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

extension = Extension("libcasimir", ["libcasimir.pyx"],
        include_dirs = [".."],
        library_dirs = [".."],
        libraries = ["casimir", "hodlr", "lapack"],
)

setup(
    ext_modules = cythonize(extension)
)
