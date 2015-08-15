/**
 * @file   matrix.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   April, 2015
 * @brief  matrix functions
 */

#include <stdio.h>

#include "libcasimir.h"
#include "sfunc.h"
#include "matrix.h"
#include "edouble.h"

MATRIX_ALLOC(matrix_sign, matrix_sign_t, sign_t);
MATRIX_FREE (matrix_sign, matrix_sign_t);

MATRIX_ALLOC(matrix_edouble, matrix_edouble_t, edouble);
MATRIX_FREE (matrix_edouble, matrix_edouble_t);
MATRIX_LOGDET_QR (matrix_edouble, matrix_edouble_t, edouble, fabse, copysigne, sqrte, loge);
//MATRIX_LOGDET_LU (matrix_edouble, matrix_edouble_t, edouble, fabse, loge);
MATRIX_ABSMIN    (matrix_edouble, matrix_edouble_t, edouble, fabse);
MATRIX_ABSMAX    (matrix_edouble, matrix_edouble_t, edouble, fabse);
MATRIX_BALANCE   (matrix_edouble, matrix_edouble_t, edouble, fabse);
MATRIX_LOG_BALANCE(matrix_edouble, matrix_edouble_t, edouble, loge);
MATRIX_EXP(matrix_edouble, matrix_edouble_t, expe);


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
            a[i+dim*j] = matrix_get(signs, i,j)*expe(matrix_get(M,i,j));

    dgetrf_(&m, &n, a, &lda, ipiv, &info);

    for(i = 0; i < dim; i++)
        logdet += log(fabs(a[i+dim*i]));

    xfree(a);
    xfree(ipiv);

    return logdet;
}
#endif
