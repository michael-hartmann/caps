# libhodlr

This is a simple wrapper for the HODLR library
https://github.com/sivaramambikasaran/HODLR (commit
fb5d7393616562e877325348b4be1abfe5281cdb) which implements the calculation of
the logarithm of the determinant of a matrix D=Id-M.

The HODLR library also depends on Eigen. Eigen is a C++ template library for
linear algebra. In order to reduce dependencies and due to compability, we
include [Eigen](https://eigen.tuxfamily.org) in version 3.3.7.
