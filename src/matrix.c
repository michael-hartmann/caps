#include <stdio.h>

#include "libcasimir.h"
#include "sfunc.h"
#include "matrix.h"
#include "edouble.h"

MATRIX_ALLOC(matrix_sign, matrix_sign_t, sign_t);
MATRIX_FREE (matrix_sign, matrix_sign_t);

MATRIX_ALLOC(matrix_edouble, matrix_edouble_t, edouble);
MATRIX_FREE (matrix_edouble, matrix_edouble_t);
MATRIX_LOGDET_QR (matrix_edouble, matrix_edouble_t, edouble, fabsq, copysignq, sqrtq, logq);
//MATRIX_LOGDET_LU (matrix_edouble, matrix_edouble_t, edouble, fabsq, logq);
MATRIX_ABSMIN    (matrix_edouble, matrix_edouble_t, edouble, fabsq);
MATRIX_ABSMAX    (matrix_edouble, matrix_edouble_t, edouble, fabsq);
MATRIX_BALANCE   (matrix_edouble, matrix_edouble_t, edouble, fabsq);
MATRIX_LOG_BALANCE(matrix_edouble, matrix_edouble_t, edouble, logq);
MATRIX_EXP(matrix_edouble, matrix_edouble_t, expq);


#ifdef USE_LAPACK
double matrix_logdet_lapack(matrix_edouble_t *M, matrix_sign_t *signs)
{
    const int dim = M->size;
    int i,j, m = dim, n = dim, lda = dim, info;
    int *ipiv = xmalloc(dim*sizeof(int));
    double *a = xmalloc(dim*dim*sizeof(double));
    double logdet = 0;

    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            a[i+dim*j] = matrix_get(signs, i,j)*expq(matrix_get(M,i,j));

    dgetrf_(&m, &n, a, &lda, ipiv, &info);

    logdet = 0;
    for(i = 0; i < dim; i++)
        logdet += log(fabs(a[i+dim*i]));

    xfree(a);
    xfree(ipiv);

    return logdet;
}
#endif
