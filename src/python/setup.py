from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

sources = [
    "libcasimir.pyx",
    "../integration.c",
    "../libcasimir.c",
    "../matrix.c",
    "../sfunc.c",
    "../utils.c"
]

CFLAGS = [
    "-std=c99",
    "-fgnu89-inline",
    "-Wall",
    "-Wextra",
    "-O3",
    "-march=native",
    "-I..",
    "-lm",
    "-pthread",
]

LIBRARIES = [ "m", "lapack", "blas", "gfortran" ]

extensions = [Extension("libcasimir", sources, extra_compile_args=CFLAGS, libraries=LIBRARIES)]

setup(
    ext_modules = cythonize(extensions)
)
