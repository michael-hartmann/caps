from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

sources = [
    "libcasimir.pyx",
    "../gausslaguerre.c",
    "../integration_drude.c",
    "../integration_perf.c",
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

extensions = [Extension("libcasimir", sources, extra_compile_args=CFLAGS)]

setup(
    ext_modules = cythonize(extensions)
)
