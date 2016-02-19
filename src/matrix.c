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


/**
 * @brief swap rows and columns
 *
 * This function is used for pivoting a matrix.
 *
 * First, columns i and j are swapped, then rows i and j are swapped.
 *
 * @param [in,out] M matrix
 * @param [in] i row
 * @param [in] j column
 */
void matrix_float80_swap(matrix_float80 *M, const int i, const int j)
{
    const int dim = M->dim;

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


/**
 * @brief Pivot matrix
 *
 * This function pivots a matrix using swapping of rows and columns. After
 * pivoting \f$|M_{00}| < |M_{11}| < |M_{22}| \dots f$.
 *
 * For some reason this pivoting contradicts what the literature suggests.
 * However, this way of pivoting seems to work fine.  One should probably read
 * carefully [1].
 *
 * [1] Algebraic and numerical techniques for the computation of matrix
 * determinants, Pan, Yu, Stewart, Computers & Mathematics with Applications,
 * 1997, http://dx.doi.org/10.1016/S0898-1221(97)00097-7
 *
 * @param [in,out] M matrix
 */
void matrix_float80_pivot(matrix_float80 *M)
{
    const int dim = M->dim;

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


/**
 * @brief calculate QR decomposition of matrix
 *
 * This function calculates the QR decomposition of a matrix M and \f$\log|\det(M)|\f$.
 *
 * @param [in,out] M matrix
 * @retval logdet \f$\log|\det M|\f$
 */
double matrix_float80_logdet_qr(matrix_float80 *M)
{
    const int dim = M->dim;
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
    const int dim = A->dim;
    bool converged = false;

    float80 *M = A->M;
    float80 *list_row = xmalloc(dim*sizeof(float80));
    float80 *list_col = xmalloc(dim*sizeof(float80));

    /* line 2 */
    while(!converged)
    {
        /* line 4 */
        converged = true;

        /* line 5 */
        for(int i = 0; i < dim; i++)
        {
            /* line 6 */
            for(int j = 0; j < dim; j++)
            {
                list_row[j] = matrix_get(A,i,j);
                list_col[j] = matrix_get(A,j,i);
            }

            float80 row_norm = logadd_m(list_row, dim);
            float80 col_norm = logadd_m(list_col, dim);

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
                for(int k = 0; k < dim; k++)
                {
                    M[i*dim+k] -= (float80)f*LOG_FLOAT_RADIX;
                    M[k*dim+i] += (float80)f*LOG_FLOAT_RADIX;
                }
            }
        }
    }

    xfree(list_col);
    xfree(list_row);
}

/* balance a matrix that elements are give by log with stop criterion */
void matrix_float80_log_balance_fast(matrix_float80 *A);

void matrix_float80_log_balance_fast(matrix_float80 *A)
{
    const int dim = A->dim;
    float80 list_row[dim];
    float80 list_col[dim];

    float80 *M = A->M;

    /* line 2 */
    while(true)
    {
        float80 max = matrix_get(A, 0,0); /* value */
        float80 min = matrix_get(A, 0,0); /* value */
        int index_max[2] = {0,0};
        int index_min[2] = {0,0};

        for(int i = 0; i < dim; i++)
            for(int j = 0; j < dim; j++)
            {
                const float80 elem = matrix_get(A, i,j);

                if(elem > max)
                {
                    max = elem;
                    index_max[0] = i;
                    index_max[1] = j;
                }
                else if(elem < min)
                {
                    min = elem;
                    index_min[0] = i;
                    index_min[1] = j;
                }
            }

        /* stop criterion */
        if(max < 5000 && min > -5000)
            break;

        int indices[4] = { index_max[0], index_max[1], index_min[0], index_min[1] };

        /* line 5 */
        for(int n = 0; n < 4; n++)
        {
            const int i = indices[n];

            /* line 6 */
            for(int j = 0; j < dim; j++)
            {
                list_row[j] = matrix_get(A,i,j);
                list_col[j] = matrix_get(A,j,i);
            }

            float80 row_norm = logadd_m(list_row, dim);
            float80 col_norm = logadd_m(list_col, dim);

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
            if(logadd(row_norm, col_norm) < (log(0.97)+s))
            {
                /* line 14 */
                for(int k = 0; k < dim; k++)
                {
                    M[i*dim+k] -= (float80)f*LOG_FLOAT_RADIX;
                    M[k*dim+i] += (float80)f*LOG_FLOAT_RADIX;
                }
            }
        }
    }
}

/* calculate log(det(1-M)) */
double matrix_logdet1mM(casimir_t *casimir, matrix_float80 *M, matrix_sign_t *M_sign)
{
    const char *detalg = casimir->detalg;
    const bool pivot   = casimir->pivot;
    #define TRACE

    const int dim = M->dim;
    #ifdef TRACE
    float80 minimum, maximum;

    printf("now=%.1f\n", now());
    matrix_float80_minmax(M, &minimum, &maximum);
    printf("# before balancing: min=%Lg, max=%Lg\n", minimum, maximum);
    #endif

    /* balance matrix */
    matrix_float80_log_balance_fast(M);

    #ifdef TRACE
    matrix_float80_minmax(M, &minimum, &maximum);
    printf("# after balancing: min=%Lg, max=%Lg\n", minimum, maximum);
    #endif

    if(strcasecmp(detalg, "LU_FLOAT80") == 0)
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
    else if(strcasecmp(detalg, "QR_FLOAT128") == 0)
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
        if(strcasecmp(detalg, "QR_FLOAT80") != 0)
            WARN(1, "Algorithm \"%s\" not supported. Defaulting to QR_FLOAT80.", detalg);

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
