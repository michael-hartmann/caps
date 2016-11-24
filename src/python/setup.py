from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

sources = [
    "libcasimir.pyx",
    "../hash-table.c",
    "../integration.c",
    "../libcasimir.c",
    "../matrix.c",
    "../sfunc.c",
    "../utils.c",
    "../lookup.c",
    "../cquadpack/src/dqagi.c",
    "../cquadpack/src/dqags.c",
    "../cquadpack/src/dqext.c",
    "../cquadpack/src/dqk15i.c",
    "../cquadpack/src/dqk21.c",
    "../cquadpack/src/dqsort.c"
]

CFLAGS = [
    "-std=c99",
    "-fgnu89-inline",
    "-Wall",
    "-Wextra",
    "-Wno-strict-prototypes",
    "-Wno-unused-parameter",
    "-O3",
    "-march=native",
    "-I..",
    "-I../cquadpack/include",
    "-lm",
    "-pthread",
]

LIBRARIES = [ "m", "lapack", "blas", "gfortran" ]

extensions = [Extension("libcasimir", sources, extra_compile_args=CFLAGS, libraries=LIBRARIES)]

setup(
    ext_modules = cythonize(extensions)
)
