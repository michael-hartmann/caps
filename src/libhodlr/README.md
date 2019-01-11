# libhodlr

This is a simple wrapper for the [HODLR
library](https://github.com/sivaramambikasaran/HODLR) in version
[3.1415](https://github.com/sivaramambikasaran/HODLR/releases/tag/3.1415) which
implements the calculation of the logarithm of the determinant of a matrix
D=Id-M. The only difference is that we `max_tries` in `HODLR_Matrix.cpp` is set
to 100. In the vanilla sources `max_tries=10`. Here is a diff:
```diff
diff --git a/src/libhodlr/src/HODLR_Matrix.cpp b/src/libhodlr/src/HODLR_Matrix.cpp
index 0e8f02c..7e32ead 100644
--- a/src/libhodlr/src/HODLR_Matrix.cpp
+++ b/src/libhodlr/src/HODLR_Matrix.cpp
@@ -107,7 +107,7 @@ void HODLR_Matrix::rookPiv(int n_row_start, int n_col_start,
     dtype_base row_squared_norm, row_norm, col_squared_norm, col_norm;

     // So these would be particularly useful for poorly conditioned matrices:
-    int max_tries = 10;
+    int max_tries = 100;
     int count;

     // Repeat till the desired tolerance is obtained
```

The HODLR library depends on Eigen. Eigen is a C++ template library for linear
algebra. In order to reduce dependencies and to improve compability, we include
[Eigen](https://eigen.tuxfamily.org) in version 3.3.7.
