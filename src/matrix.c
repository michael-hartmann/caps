/**
 * @file   matrix.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   January, 2016
 * @brief  matrix functions
 */

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include "libcasimir.h"
#include "sfunc.h"
#include "matrix.h"
#include "edouble.h"

MATRIX_ALLOC(matrix_sign, matrix_sign_t, sign_t);
MATRIX_FREE (matrix_sign, matrix_sign_t);
MATRIX_SAVE (matrix_sign, matrix_sign_t, sign_t);
MATRIX_LOAD (matrix_sign, matrix_sign_t, sign_t, matrix_sign_alloc, matrix_sign_free);

MATRIX_ALLOC(matrix_edouble, matrix_edouble_t, edouble);
MATRIX_FREE (matrix_edouble, matrix_edouble_t);
MATRIX_SAVE (matrix_edouble, matrix_edouble_t, edouble);
MATRIX_LOAD (matrix_edouble, matrix_edouble_t, edouble, matrix_edouble_alloc, matrix_edouble_free);

MATRIX_LOGDET_QR (matrix_edouble, matrix_edouble_t, edouble, fabse, copysigne, sqrte, loge);
MATRIX_LOGDET_LU (matrix_edouble, matrix_edouble_t, edouble, fabse, loge);
MATRIX_EXP(matrix_edouble, matrix_edouble_t, expe);

void matrix_edouble_log_balance(matrix_edouble_t *A, const int p)
{
    const int N = A->size;
    bool converged = false;

    edouble *M = A->M;
    edouble *list_row = xmalloc(N*sizeof(edouble));
    edouble *list_col = xmalloc(N*sizeof(edouble));

    /* line 2 */
    while(!converged)
    {
        /* line 4 */
        converged = true;

        /* line 5 */
        for(int i = 0; i < N; i++)
        {
            /* line 6 */
            for(int j = 0; j < N; j++)
            {
                list_row[j] = p*matrix_get(A,i,j);
                list_col[j] = p*matrix_get(A,j,i);
            }

            edouble row_norm = logadd_m(list_row, N);
            edouble col_norm = logadd_m(list_col, N);

            /*
            if(isinf(row_norm) || isinf(col_norm))
                continue;
            */

            /* line 7 */
            int f = 0; /* log(1)=0 */
            edouble s = logadd(p*col_norm, p*row_norm);

            /* line 8 */
            while(col_norm < (row_norm-LOG_FLOAT_RADIX))
            {
                /* line 9 */
                col_norm += LOG_FLOAT_RADIX;
                row_norm -= LOG_FLOAT_RADIX;
                f        += 1;
            }

            /* line 10 */
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
                converged = false;

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
        matrix_edouble_exp(M, M_sign);
        return matrix_edouble_logdet_qr(M);
    }
    else if(strcasecmp(type, "LU") == 0)
    {
        matrix_edouble_log_balance(M,1);
        matrix_edouble_exp(M, M_sign);
        return matrix_edouble_logdet_lu(M);
    }
    else
    {
        TERMINATE(1, "Algorithm not supported: %s.", type);
        return 0;
    }
}
