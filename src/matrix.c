/**
 * @file   matrix.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   September, 2015
 * @brief  matrix functions
 */

#ifdef USE_LAPACK
#include <cblas.h>
#endif
#include <stdio.h>
#include <string.h>

#include "libcasimir.h"
#include "sfunc.h"
#include "matrix.h"
#include "edouble.h"

MATRIX_ALLOC(matrix_sign, matrix_sign_t, sign_t);
MATRIX_FREE (matrix_sign, matrix_sign_t);

MATRIX_ALLOC(matrix_edouble, matrix_edouble_t, edouble);
MATRIX_FREE (matrix_edouble, matrix_edouble_t);
MATRIX_LOGDET_QR (matrix_edouble, matrix_edouble_t, edouble, fabse, copysigne, sqrte, loge);
MATRIX_LOGDET_LU (matrix_edouble, matrix_edouble_t, edouble, fabse, loge);
MATRIX_ABSMIN    (matrix_edouble, matrix_edouble_t, edouble, fabse);
MATRIX_ABSMAX    (matrix_edouble, matrix_edouble_t, edouble, fabse);
MATRIX_LOG_BALANCE(matrix_edouble, matrix_edouble_t, edouble, loge);
MATRIX_EXP(matrix_edouble, matrix_edouble_t, expe);

double matrix_edouble_logdet(matrix_edouble_t *M, matrix_sign_t *M_sign, const char *type)
{
    if(strcasecmp(type, "QR") == 0)
    {
        matrix_edouble_log_balance(M);
        matrix_edouble_exp(M, M_sign);
        return matrix_edouble_logdet_qr(M);
    }
    else if(strcasecmp(type, "LU") == 0)
    {
        matrix_edouble_log_balance(M);
        matrix_edouble_exp(M, M_sign);
        return matrix_edouble_logdet_lu(M);
    }
    #ifdef USE_LAPACK
    else if(strcasecmp(type, "block_LAPACK") == 0)
    {
        matrix_edouble_log_balance(M);
        return matrix_logdet_block_lapack(M, M_sign);
    }
    else if(strcasecmp(type, "LU_LAPACK") == 0)
    {
        matrix_edouble_log_balance(M);
        return matrix_logdet_lu_lapack(M, M_sign);
    }
    #endif
    else
    {
        TERMINATE(0, "Algorithm not supported: %s.", type);
        return 0;
    }
}

#ifdef USE_LAPACK
static void matrix_mult(double *A, double *B, double *C, double alpha, double beta, int dim)
{
    char notrans = 'N';
    dgemm_(
        &notrans,  /* don't transpose/conjugate A */
        &notrans,  /* don't transpose/conjugate B */
        &dim,      /* M: rows of A and C */
        &dim,      /* N: columns of B and C */
        &dim,      /* K: columns of A and rows of B */
        &alpha,    /* alpha: scalar */
        A,         /* A: matrix A */
        &dim,      /* LDA: leading dimension of A (columns) */
        B,         /* B: matrix B */
        &dim,      /* LDB: leading dimension of B (columns) */
        &beta,     /* beta: scalar */
        C,         /* C: matrix C */
        &dim       /* LDC: leading dimension of C (columns) */
    );
}


static int invert(double *M, int dim)
{
    int info, lwork;
    int *ipiv = NULL;
    double *work = NULL;
    double workopt;

    ipiv = xmalloc(dim*sizeof(int));

    dgetrf_(
        &dim, /* M number of rows of A */
        &dim, /* N number of columns of A */
        M,    /* matrix A to be factored */
        &dim, /* LDA: leading dimension of A */
        ipiv, /* pivot indices of dimension (min(M,N)) */
        &info
    );

    lwork = -1;
    dgetri_(&dim, M, &dim, ipiv, &workopt, &lwork, &info);
    lwork = (int)workopt;
    work = xmalloc(lwork*sizeof(double));

    dgetri_(
        &dim,   /* order of matrix A */
        M,      /* factors L and U from LU decomposition */
        &dim,   /* LDA: leading dimension of A */
        ipiv,   /* pivot indices */
        work,   /* workspace of dimension LWORK */
        &lwork, /* length of work */
        &info
    );

    xfree(ipiv);
    xfree(work);

    return info;
}

double matrix_logdet_lu_lapack(matrix_edouble_t *M, matrix_sign_t *signs)
{
    const int dim = M->size;
    int i,j, m = dim, n = dim, lda = dim, info;
    int *ipiv = xmalloc(dim*sizeof(int));
    double *a = xmalloc(dim*dim*sizeof(double));
    double logdet = 0;

    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            a[i+dim*j] = matrix_get(signs, i,j)*exp(matrix_get(M,i,j));

    dgetrf_(&m, &n, a, &lda, ipiv, &info);

    for(i = 0; i < dim; i++)
        logdet += log(fabs(a[i+dim*i]));

    xfree(a);
    xfree(ipiv);

    return logdet;
}

double matrix_logdet_block_lapack(matrix_edouble_t *M, matrix_sign_t *signs)
{
    const int dim = M->size/2;
    int i,j, m = dim, n = dim, lda = dim, info;
    int *ipiv = xmalloc(dim*sizeof(int));
    double *a = xmalloc(dim*dim*sizeof(double));
    double *b = xmalloc(dim*dim*sizeof(double));
    double *c = xmalloc(dim*dim*sizeof(double));
    double *d = xmalloc(dim*dim*sizeof(double));
    double *dinv = xmalloc(dim*dim*sizeof(double));
    double logdet = 0;

    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
        {
            a[i+dim*j] = matrix_get(signs, i,j)    *exp(matrix_get(M,i,j));
            b[i+dim*j] = matrix_get(signs, i,j+dim)*exp(matrix_get(M,i,j+dim));
            c[i+dim*j] = matrix_get(signs, i+dim,j)*exp(matrix_get(M,i+dim,j));
            d[i+dim*j] = dinv[i+dim*j] = matrix_get(signs, i+dim,j+dim)*exp(matrix_get(M,i+dim,j+dim));
        }

    /* calculate logdet of matrix d */
    dgetrf_(&m, &n, d, &lda, ipiv, &info);
    for(i = 0; i < dim; i++)
        logdet += log(fabs(d[i+dim*i]));

    /* invert d */
    invert(dinv, dim);

    /* dinv*c -> d */
    matrix_mult(dinv, c, d, 1, 0, dim);

    /* a - b*d = a - b*dinv*c */
    matrix_mult(b, d, a, -1, 1, dim);

    /* calculate logdet of matrix a */
    dgetrf_(&m, &n, a, &lda, ipiv, &info);
    for(i = 0; i < dim; i++)
        logdet += log(fabs(a[i+dim*i]));

    xfree(a);
    xfree(b);
    xfree(c);
    xfree(d);
    xfree(dinv);
    xfree(ipiv);

    return logdet;
}
#endif
