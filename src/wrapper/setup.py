from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

extension = Extension("libcasimir", ["libcasimir.pyx"],
        include_dirs = [".."],
        libraries = ["casimir", "hodlr"],
        library_dirs = ["..", "/home/michael/git/hodlr_wrapper"])

setup(
    ext_modules = cythonize(extension)
)
