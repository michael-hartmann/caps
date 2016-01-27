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

extern double matrix_floatdd_logdet(matrix_float80 *M, matrix_sign_t *M_sign, const char *type);

/*
MATRIX_ALLOC(matrix_sfloat, matrix_sfloat_t, sfloat_t);
MATRIX_FREE (matrix_sfloat, matrix_sfloat_t);
MATRIX_SAVE (matrix_sfloat, matrix_sfloat_t, sfloat_t);
MATRIX_LOAD (matrix_sfloat, matrix_sfloat_t, sfloat_t, matrix_sfloat_alloc, matrix_sfloat_free);
*/

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

                for(int n = 0; n < dim; n++)
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

                for(int n = 0; n < dim; n++)
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

#if 0
double matrix_sfloat_logdet_qr(matrix_sfloat_t *M)
{
    const int dim = M->size;
    sfloat_t zero,one;
    sfloat_from_double(&zero, 0);
    sfloat_from_double(&one,  1);

    sfloat_t *m = M->M;

    for(int j = 0; j < dim-1; j++)
        for(int i = j+1; i < dim; i++)
        {
            sfloat_t c,s;
            sfloat_t *Mij = &m[i*dim+j];


            if(!sfloat_iszero(Mij))
            {
                sfloat_t *a = &m[j*dim+j];
                sfloat_t *b = Mij;

                if(sfloat_iszero(b))
                {
                    sfloat_copysign(&one, a, &c);
                    s = zero;
                }
                else if(sfloat_iszero(a))
                {
                    c = zero;
                    sfloat_copysign(&one, b, &s);
                    sfloat_neg(&s);
                }
                else if(sfloat_abscompare(b, a) > 0)
                {
                    sfloat_t t,u,v,w;

                    /* t = a/b */
                    sfloat_div(a, b, &t);

                    /* v = sqrt(1+t*t)
                     * u = copysign80(v,b);
                     */
                    sfloat_mul(&t, &t, &v);
                    sfloat_add(&one, &v, &w);
                    sfloat_sqrt(&w);
                    sfloat_copysign(&w, b, &u);

                    /* s = -1/u; */
                    sfloat_div(&one, &u, &s);
                    sfloat_neg(&s);

                    /* c = -s*t; */
                    sfloat_mul(&s, &t, &c);
                    sfloat_neg(&c);
                }
                else
                {
                    sfloat_t t,u,v,w;

                    /* t = b/a*/
                    sfloat_div(b,a,&t);

                    /* v = sqrt(1+t*t)
                     * u = copysign80(v,a)
                     */
                    sfloat_mul(&t,&t, &v);
                    sfloat_add(&one, &v, &w);
                    sfloat_sqrt(&w);
                    sfloat_copysign(&w,a,&u);

                    /* c = 1/u */
                    sfloat_div(&one, &u, &c);

                    /* s = -c*t */
                    sfloat_mul(&c,&t,&s);
                    sfloat_neg(&s);
                }

                for(int n = 0; n < dim; n++)
                {
                    sfloat_t t1,t2;
                    sfloat_t Min = m[i*dim+n];
                    sfloat_t Mjn = m[j*dim+n];

                    /* m[i*dim+n] = c*Min + s*Mjn; */
                    sfloat_mul(&c,&Min,&t1);
                    sfloat_mul(&s,&Mjn,&t2);
                    sfloat_add(&t1, &t2, &m[i*dim+n]);

                    /* m[j*dim+n] = c*Mjn - s*Min; */
                    sfloat_mul(&c,&Mjn,&t1);
                    sfloat_mul(&s,&Min,&t2);
                    sfloat_sub(&t1,&t2,&m[j*dim+n]);
                }
            }
        }

    double det = 0;
    for(int i = 0; i < dim; i++)
    {
        sfloat_t *x = &m[i*dim+i];
        det += log(fabs(x->mantisse))+x->exponent*M_LN2;
    }

    return det;
}
#endif

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
    if(strcasecmp(type, "QR") == 0)
    {
        matrix_float80_log_balance(M);
        matrix_float80_exp(M, M_sign);
        return matrix_float80_logdet_qr(M);
    }
    else if(strcasecmp(type, "LU") == 0)
    {
        matrix_float80_log_balance(M);
        matrix_float80_exp(M, M_sign);
        return matrix_float80_logdet_lu(M);
    }
    else if(strcasecmp(type, "DEBUG") == 0)
    {
        double tstart, tend;
        double logdet_dd, logdet_float80;
        #ifdef FLOAT128
        double logdet_float128;
        #endif
        float80 minimum, maximum;

        matrix_float80_minmax(M, &minimum, &maximum);
        printf("# before balancing: min=%Lg, max=%Lg\n", (long double)minimum, (long double)maximum);

        tstart = now();
        matrix_float80_log_balance(M);
        tend = now();
        matrix_float80_minmax(M, &minimum, &maximum);
        printf("# after balancing:  min=%Lg, max=%Lg, (t=%g)\n", (long double)minimum, (long double)maximum, tend-tstart);

        tstart = now();
        logdet_dd = matrix_floatdd_logdet(M, M_sign, type);
        tend = now();
        printf("# matrix_floatdd_logdet:  %+.15g (%g)\n", logdet_dd, tend-tstart);

        #ifdef FLOAT128
        const int dim = M->size;
        tstart = now();
        matrix_float128 *M128 = matrix_float128_alloc(dim);
        for(int i = 0; i < dim*dim; i++)
            M128->M[i] = M_sign->M[i] * exp128(M->M[i]);

        logdet_float128 = matrix_float128_logdet_qr(M128);
        matrix_float128_free(M128);
        tend = now();
        printf("# matrix_float128_logdet: %+.15g (%g)\n", logdet_float128, tend-tstart);
        #endif

        tstart = now();
        matrix_float80_exp(M, M_sign);
        logdet_float80 = matrix_float80_logdet_qr(M);
        tend = now();
        printf("# matrix_float80_logdet:  %+.15g (%g)\n", logdet_float80, tend-tstart);

        return logdet_float80;
    }
    else
    {
        TERMINATE(1, "Algorithm not supported: %s.", type);
        return 0;
    }
}
