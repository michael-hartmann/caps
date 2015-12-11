/**
 * @file   matrix.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   December, 2015
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
MATRIX_EXP(matrix_edouble, matrix_edouble_t, expe);

void matrix_edouble_log_balance(matrix_edouble_t *A)
{
    size_t i,j;
    const size_t N = A->size;
    int converged = 0;
    const int p = 1;

    edouble *M = A->M;
    edouble *list_row = xmalloc(N*sizeof(edouble));
    edouble *list_col = xmalloc(N*sizeof(edouble));

    /* line 2 */ \
    while(!converged)
    {
        int f;
        edouble s;
        edouble row_norm, col_norm;

        /* line 4 */
        converged = 1;

        /* line 5 */
        for(i = 0; i < N; i++)
        {
            /* line 6 */ \
            for(j = 0; j < N; j++)
            {
                const edouble Aij = matrix_get(A,i,j);
                const edouble Aji = matrix_get(A,j,i);

                list_row[j] = p*Aij;
                list_col[j] = p*Aji;
            }

            row_norm = logadd_m(list_row, N);
            col_norm = logadd_m(list_col, N);

            if(isinf(row_norm) || isinf(col_norm))
                continue;

            /* line 7 */
            f = 0; /* log(1)=0 */
            s = logadd(p*col_norm, p*row_norm);

            /* line 8 */
            while(col_norm < (row_norm-LOG_FLOAT_RADIX))
            {
                /* line 9 */ \
                col_norm += LOG_FLOAT_RADIX;
                row_norm -= LOG_FLOAT_RADIX;
                f        += 1;
            }

            /* line 10 */ \
            while(col_norm >= (row_norm+LOG_FLOAT_RADIX))
            {
                /* line 11 */
                col_norm -= LOG_FLOAT_RADIX;
                row_norm += LOG_FLOAT_RADIX;
                f        -= 1;
            }

            /* line 12 */
            if(logadd(p*row_norm, p*col_norm) < (log(0.95)+s))
            {
                int k;
                /* line 13 */
                converged = 0;

                /* line 14 */
                for(k = 0; k < N; k++)
                {
                    M[i*N+k] -= (edouble)f*LOG_FLOAT_RADIX;
                    M[k*N+i] += (edouble)f*LOG_FLOAT_RADIX;
                }
            }
        }
    }

    xfree(list_col);
    xfree(list_row);
}

double matrix_edouble_logdet(matrix_edouble_t *M, matrix_sign_t *M_sign, const char *type)
{
    if(strcasecmp(type, "QR") == 0)
    {
        int i,j;
        edouble max, min;
        //edouble avg;

        max = min = matrix_get(M,0,0);
        for(i = 0; i < M->size; i++)
            for(j = 0; j < M->size; j++)
            {
                edouble elem = matrix_get(M,i,j);
                if(elem > max)
                    max = elem;
                if(elem < min)
                    min = elem;
            }
        printf("# vor balancing: %Lg, %Lg\n", min,max);

        matrix_edouble_log_balance(M);

        max = min = matrix_get(M,0,0);
        for(i = 0; i < M->size; i++)
            for(j = 0; j < M->size; j++)
            {
                edouble elem = matrix_get(M,i,j);
                if(elem > max)
                    max = elem;
                if(elem < min)
                    min = elem;

            }
        printf("# nach balancing: %Lg, %Lg\n", min,max);
/*
        avg = (max-min)/2;

        for(i = 0; i < M->size; i++)
            for(j = 0; j < M->size; j++)
            {
                edouble elem = matrix_get(M,i,j);
                matrix_set(M,i,j,elem+avg);
            }
            */

        matrix_edouble_exp(M, M_sign);
        return matrix_edouble_logdet_qr(M);// - M->size*avg;
    }
    else if(strcasecmp(type, "LU") == 0)
    {
        matrix_edouble_log_balance(M);
        matrix_edouble_exp(M, M_sign);
        return matrix_edouble_logdet_lu(M);
    }
    #ifdef USE_LAPACK
    else if(strcasecmp(type, "LU_LAPACK") == 0)
    {
        matrix_edouble_log_balance(M);
        return matrix_logdet_lu_lapack(M, M_sign);
    }
    #endif
    else
    {
        TERMINATE(1, "Algorithm not supported: %s.", type);
        return 0;
    }
}

#ifdef USE_LAPACK
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

#endif
