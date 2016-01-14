/**
 * @file   matrix.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   January, 2016
 * @brief  matrix functions
 */

#include <stdio.h>
#include <string.h>
#include <strings.h>

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

void matrix_edouble_log_balance(matrix_edouble_t *A, const int p)
{
    const int N = A->size;
    int converged = 0;

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
        for(int i = 0; i < N; i++)
        {
            /* line 6 */ \
            for(int j = 0; j < N; j++)
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
            if(logadd(p*row_norm, p*col_norm) < (LOG095+s))
            {
                /* line 13 */
                converged = 0;

                /* line 14 */
                for(int k = 0; k < N; k++)
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
        matrix_edouble_log_balance(M,1);
        matrix_edouble_log_balance(M,2);
        matrix_edouble_exp(M, M_sign);
        return matrix_edouble_logdet_qr(M);
    }
    else if(strcasecmp(type, "LU") == 0)
    {
        matrix_edouble_log_balance(M,1);
        matrix_edouble_log_balance(M,2);
        matrix_edouble_exp(M, M_sign);
        return matrix_edouble_logdet_lu(M);
    }
    else
    {
        TERMINATE(1, "Algorithm not supported: %s.", type);
        return 0;
    }
}
