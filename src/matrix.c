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
#include "floattypes.h"

#include "utils.h"

extern double matrix_floatdd_logdet(matrix_float80 *M, matrix_sign_t *M_sign);

#ifdef FLOAT128
MATRIX_ALLOC (matrix_float128, matrix_float128, float128);
MATRIX_FREE  (matrix_float128, matrix_float128);
MATRIX_SAVE  (matrix_float128, matrix_float128, float128);
MATRIX_LOAD  (matrix_float128, matrix_float128, float128, matrix_float128_alloc, matrix_float128_free);
MATRIX_MINMAX(matrix_float128, matrix_float128, float128);
#endif

MATRIX_ALLOC(matrix_sign, matrix_sign_t, sign_t);
MATRIX_FREE (matrix_sign, matrix_sign_t);
MATRIX_SAVE (matrix_sign, matrix_sign_t, sign_t);
MATRIX_LOAD (matrix_sign, matrix_sign_t, sign_t, matrix_sign_alloc, matrix_sign_free);

MATRIX_ALLOC (matrix_float80, matrix_float80, float80);
MATRIX_FREE  (matrix_float80, matrix_float80);
MATRIX_SAVE  (matrix_float80, matrix_float80, float80);
MATRIX_LOAD  (matrix_float80, matrix_float80, float80, matrix_float80_alloc, matrix_float80_free);
MATRIX_MINMAX(matrix_float80, matrix_float80, float80);

MATRIX_LOGDET_LU (matrix_float80, matrix_float80, float80, fabs80, log80);
MATRIX_EXP(matrix_float80, matrix_float80, exp80);

#ifdef FLOAT128
double matrix_float128_logdet_qr(matrix_float128 *M)
{
    const int dim = M->size;
    float128 *m = M->M;

    for(int j = 0; j < dim-1; j++)
        for(int i = j+1; i < dim; i++)
        {
            float128 c,s, Mij = m[i*dim+j];

            if(Mij != 0)
            {
                const float128 a = m[j*dim+j];
                const float128 b = Mij;

                if(b == 0)
                {
                    c = copysign128(1,a);
                    s = 0;
                }
                else if(a == 0)
                {
                    c = 0;
                    s = -copysign128(1, b);
                }
                else if(fabs128(b) > fabs128(a))
                {
                    const float128 t = a/b;
                    const float128 u = copysign128(sqrt128(1+t*t),b);
                    s = -1/u;
                    c = -s*t;
                }
                else
                {
                    const float128 t = b/a;
                    const float128 u = copysign128(sqrt128(1+t*t),a);
                    c = 1/u;
                    s = -c*t;
                }

                for(int n = j; n < dim; n++)
                {
                    const float128 Min = m[i*dim+n];
                    const float128 Mjn = m[j*dim+n];

                    m[i*dim+n] = c*Min + s*Mjn;
                    m[j*dim+n] = c*Mjn - s*Min;
                }

                /* m[i*dim+j] = 0; */
            }
        }

    float128 det = 0;
    for(int i = 0; i < dim; i++)
        det += log128(fabs128(m[i*dim+i]));

    return det;
}
#endif


double matrix_float80_log_logdet_qr(matrix_float80 *M, matrix_sign_t *M_sign)
{
    const int dim = M->size;
    float80 *m = M->M;
    sign_t *m_sign = M_sign->M;

    for(int j = 0; j < dim-1; j++)
        for(int i = j+1; i < dim; i++)
        {
            float80 c;
            const float80 t = m[i*dim+j]-m[j*dim+j];
            const sign_t sign_c =  m_sign[j*dim+j];
            const sign_t sign_s = -m_sign[i*dim+j];

            if(t > 0)
                c = -t - log1p80(exp80(-2*t))/2;
            else
                c = -log1p80(exp80(2*t))/2;

            const float80 s = c+t;  /* sign(s) = -sign(c)*sign(t) = -sign(M[i,j]) */

            for(int n = j; n < dim; n++)
            {
                const float80 Min = m[i*dim+n];
                const float80 Mjn = m[j*dim+n];
                const sign_t sign_Min = m_sign[i*dim+n];
                const sign_t sign_Mjn = m_sign[j*dim+n];

                m[i*dim+n] = logadd_s(c+Min, sign_c*sign_Min, s+Mjn,  sign_s*sign_Mjn, &m_sign[i*dim+n]);
                m[j*dim+n] = logadd_s(c+Mjn, sign_c*sign_Mjn, s+Min, -sign_s*sign_Min, &m_sign[j*dim+n]);
            }
        }

    float80 det = 0;
    for(int i = 0; i < dim; i++)
        det += m[i*dim+i];

    return det;
}
double matrix_float80_logdet_qr(matrix_float80 *M)
{
    const int dim = M->size;
    float80 *m = M->M;

    for(int j = 0; j < dim-1; j++)
        for(int i = j+1; i < dim; i++)
        {
            float80 c,s, Mij = m[i*dim+j];

            if(Mij != 0)
            {
                const float80 a = m[j*dim+j];
                const float80 b = Mij;

                if(b == 0)
                {
                    c = copysign80(1,a);
                    s = 0;
                }
                else if(a == 0)
                {
                    c = 0;
                    s = -copysign80(1, b);
                }
                else if(fabs80(b) > fabs80(a))
                {
                    const float80 t = a/b;
                    const float80 u = copysign80(sqrt80(1+t*t),b);
                    s = -1/u;
                    c = -s*t;
                }
                else
                {
                    const float80 t = b/a;
                    const float80 u = copysign80(sqrt80(1+t*t),a);
                    c = 1/u;
                    s = -c*t;
                }

                for(int n = j; n < dim; n++)
                {
                    const float80 Min = m[i*dim+n];
                    const float80 Mjn = m[j*dim+n];

                    m[i*dim+n] = c*Min + s*Mjn;
                    m[j*dim+n] = c*Mjn - s*Min;
                }

                /* m[i*dim+j] = 0; */
            }
        }

    float80 det = 0;
    for(int i = 0; i < dim; i++)
        det += log80(fabs80(m[i*dim+i]));

    return det;
}

void matrix_float80_log_balance(matrix_float80 *A)
{
    const int N = A->size;
    bool converged = false;

    float80 *M = A->M;
    float80 *list_row = xmalloc(N*sizeof(float80));
    float80 *list_col = xmalloc(N*sizeof(float80));

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
                list_row[j] = matrix_get(A,i,j);
                list_col[j] = matrix_get(A,j,i);
            }

            float80 row_norm = logadd_m(list_row, N);
            float80 col_norm = logadd_m(list_col, N);

            /* line 7 */
            int f = 0; /* log(1)=0 */
            float80 s = logadd(col_norm, row_norm);

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
            if(logadd(row_norm, col_norm) < (LOG095+s))
            {
                /* line 13 */
                converged = false;

                /* line 14 */
                for(int k = 0; k < N; k++)
                {
                    M[i*N+k] -= (float80)f*LOG_FLOAT_RADIX;
                    M[k*N+i] += (float80)f*LOG_FLOAT_RADIX;
                }
            }
        }
    }

    xfree(list_col);
    xfree(list_row);
}

double matrix_float80_logdet(matrix_float80 *M, matrix_sign_t *M_sign, const char *type)
{
    if(strcasecmp(type, "LU_FLOAT80") == 0)
    {
        matrix_float80_log_balance(M);
        matrix_float80_exp(M, M_sign);
        return matrix_float80_logdet_lu(M);
    }
    else if(strcasecmp(type, "QR_FLOAT80") == 0)
    {
        matrix_float80_log_balance(M);
        matrix_float80_exp(M, M_sign);
        return matrix_float80_logdet_qr(M);
    }
    else if(strcasecmp(type, "QR_FLOATDD") == 0)
    {
        matrix_float80_log_balance(M);
        return matrix_floatdd_logdet(M, M_sign);
    }
    else if(strcasecmp(type, "QR_LOG80") == 0)
    {
        matrix_float80_log_balance(M);
        return matrix_float80_log_logdet_qr(M, M_sign);
    }
    #ifdef FLOAT128
    else if(strcasecmp(type, "QR_FLOAT128") == 0)
    {
        const int dim = M->size;
        matrix_float128 *M128 = matrix_float128_alloc(dim);

        matrix_float80_log_balance(M);

        for(int i = 0; i < dim*dim; i++)
            M128->M[i] = M_sign->M[i]*exp128(M->M[i]);

        const double logdet = matrix_float128_logdet_qr(M128);

        matrix_float128_free(M128);

        return logdet;
    }
    #endif
    else
    {
        TERMINATE(1, "Algorithm not supported: %s.", type);
        return 0;
    }
}
