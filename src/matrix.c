/**
 * @file   matrix.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   February, 2016
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

#ifdef FLOATDD
extern double matrix_floatdd_logdet(matrix_float80 *M, matrix_sign_t *M_sign);
#endif

/** functions for float128 */
#ifdef FLOAT128
MATRIX_ALLOC (matrix_float128, matrix_float128, float128);
MATRIX_FREE  (matrix_float128, matrix_float128);
MATRIX_SAVE  (matrix_float128, matrix_float128, float128);
MATRIX_LOAD  (matrix_float128, matrix_float128, float128, matrix_float128_alloc, matrix_float128_free);
MATRIX_MINMAX(matrix_float128, matrix_float128, float128);
#endif

/** functions for sign_t */
MATRIX_ALLOC(matrix_sign, matrix_sign_t, sign_t);
MATRIX_FREE (matrix_sign, matrix_sign_t);
MATRIX_SAVE (matrix_sign, matrix_sign_t, sign_t);
MATRIX_LOAD (matrix_sign, matrix_sign_t, sign_t, matrix_sign_alloc, matrix_sign_free);

/** functions for float80 */
MATRIX_ALLOC (matrix_float80, matrix_float80, float80);
MATRIX_FREE  (matrix_float80, matrix_float80);
MATRIX_SAVE  (matrix_float80, matrix_float80, float80);
MATRIX_LOAD  (matrix_float80, matrix_float80, float80, matrix_float80_alloc, matrix_float80_free);
MATRIX_MINMAX(matrix_float80, matrix_float80, float80);
MATRIX_EXP   (matrix_float80, matrix_float80, exp80);
MATRIX_LOGDET_LU(matrix_float80, matrix_float80, float80, fabs80, log80);

void matrix_float80_swap(matrix_float80 *M, const int i, const int j)
{
    const int dim = M->size;

    /* swap columns */
    for(int k = 0; k < dim; k++)
    {
        const float80 Mkj = matrix_get(M, k,j);
        matrix_set(M, k,j, matrix_get(M, k,i)); /* Mkj = Mki */
        matrix_set(M, k,i, Mkj);                /* Mki = Mkj */
    }

    /* swap rows */
    for(int k = 0; k < dim; k++)
    {
        const float80 Mjk = matrix_get(M, j,k);
        matrix_set(M, j,k, matrix_get(M, i,k)); /* Mjk = Mik */
        matrix_set(M, i,k, Mjk);                /* Mik = Mjk */
    }
}


void matrix_float80_pivot(matrix_float80 *M)
{
    const int dim = M->size;

    for(int k = 0; k < dim; k++)
    {
        /* find maximum */
        int index = k;
        float80 elem = matrix_get(M,k,k);

        for(int z = k+1; z < dim; z++)
        {
            const float80 Mzz = matrix_get(M,z,z);
            if(fabs80(Mzz) < fabs80(elem))
            {
                elem = Mzz;
                index = z;
            }
        }

        /* swap k <-> index */
        if(k != index)
            matrix_float80_swap(M,k,index);
    }
}

#ifdef FLOAT128
/* calculate QR decomposition of M */
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


/* calculate QR decomposition of M */
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
                const float80 b = Mij; /* b != 0 */

                if(a == 0)
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

/* balance a matrix that elements are give by log */
void matrix_float80_log_balance(matrix_float80 *A)
{
    matrix_float80_log_balance_stop(A, LOG095);
}

/* balance a matrix that elements are give by log with stop criterion */
void matrix_float80_log_balance_stop(matrix_float80 *A, const double stop)
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
            if(logadd(row_norm, col_norm) < (stop+s))
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

/* calculate log(det(1-M)) */
double matrix_logdet1mM(matrix_float80 *M, matrix_sign_t *M_sign, const char *type, const bool pivot)
{
    const int dim = M->size;
    #ifdef TRACE
    float80 minimum, maximum;

    matrix_float80_minmax(M, &minimum, &maximum);
    printf("# before balancing: min=%Lg, max=%Lg\n", minimum, maximum);
    #endif

    /* balance matrix */
    matrix_float80_log_balance(M);

    #ifdef TRACE
    matrix_float80_minmax(M, &minimum, &maximum);
    printf("# after balancing: min=%Lg, max=%Lg\n", minimum, maximum);
    #endif

    if(strcasecmp(type, "LU_FLOAT80") == 0)
    {
        /* exponentiate */
        matrix_float80_exp(M, M_sign);

        /* add unity matrix */
        for(int i = 0; i < dim; i++)
            M->M[i*dim+i] += 1;

        /* calculate log(det(M)) */
        return matrix_float80_logdet_lu(M);
    }
    #ifdef FLOAT128
    else if(strcasecmp(type, "QR_FLOAT128") == 0)
    {
        matrix_float128 *M128 = matrix_float128_alloc(dim);

        for(int i = 0; i < dim*dim; i++)
            M128->M[i] = M_sign->M[i]*exp128(M->M[i]);

        for(int i = 0; i < dim; i++)
            M128->M[i*dim+i] += 1;

        const double logdet = matrix_float128_logdet_qr(M128);

        matrix_float128_free(M128);

        return logdet;
    }
    #endif
    else
    {
        if(strcasecmp(type, "QR_FLOAT80") != 0)
            WARN(1, "Algorithm \"%s\" not supported. Defaulting to QR_FLOAT80.", type);

        matrix_float80_exp(M, M_sign);

        /* add identity matrix */
        for(int i = 0; i < dim; i++)
            M->M[i*dim+i] += 1;

        /* pivot */
        if(pivot)
        {
            #ifdef TRACE
            const double t0 = now();
            printf("pivoting...\n");
            #endif

            matrix_float80_pivot(M);

            #ifdef TRACE
            printf("pivoting: t=%g\n", now()-t0);
            #endif
        }

        return matrix_float80_logdet_qr(M);
    }
}
