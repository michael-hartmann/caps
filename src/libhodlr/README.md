# libhodlr

This is a simple wrapper for the HODLR library https://github.com/sivaramambikasaran/HODLR
which implements the calculation of the logarithm of the determinant of a matrix D=Id-M.

The HODLR library also depends on Eigen. Eigen is a C++ template library for
linear algebra. In order to reduce dependencies and due to compability, we
include Eigen in version 3.2.7.

Please note that versions of Eigen >= 3.2.8 break the library giving inaccurate
results. The reason is a change in the file FullPivLU.h described as "Make
FullPivLU::solve use rank() instead of nonzeroPivots()." in the changelog. The
diff looks like:

```
*** FullPivLU.h_327     2018-03-23 16:02:57.268313250 +0100
--- FullPivLU.h_328     2018-03-23 16:03:22.409744343 +0100
*************** struct solve_retval<FullPivLU<_MatrixTyp
*** 688,694 ****
       */
  
      const Index rows = dec().rows(), cols = dec().cols(),
!               nonzero_pivots = dec().nonzeroPivots();
      eigen_assert(rhs().rows() == rows);
      const Index smalldim = (std::min)(rows, cols);
  
--- 688,694 ----
       */
  
      const Index rows = dec().rows(), cols = dec().cols(),
!               nonzero_pivots = dec().rank();
      eigen_assert(rhs().rows() == rows);
      const Index smalldim = (std::min)(rows, cols);
```

The source file FullPivLU.h of version 3.2.7 (FullPivLU.h_327) and 3.2.8
(FullPivLU.h_328) is also in this directory.
