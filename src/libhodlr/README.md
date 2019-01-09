# libhodlr

This is a simple wrapper for the [HODLR
library](https://github.com/sivaramambikasaran/HODLR) (commit
[3211704](https://github.com/sivaramambikasaran/HODLR/commit/321170475c8df8a21fa43bfec7b567962a0a621b))
which implements the calculation of the logarithm of the determinant of a
matrix D=Id-M.

The HODLR library also depends on Eigen. Eigen is a C++ template library for
linear algebra. In order to reduce dependencies and to improve compability, we
include [Eigen](https://eigen.tuxfamily.org) in version 3.3.7.
