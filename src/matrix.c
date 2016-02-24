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

/* balance a matrix that elements are give by log with stop criterion */
void matrix_float80_log_balance(matrix_float80 *A, float80 *minimum, float80 *maximum)
{
    const int dim = A->dim;

    float80 *M = A->M;
    float80 list_row[dim];
    float80 list_col[dim];

    /* line 2 */
    while(true)
    {
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

            /* line 14 */
            for(int k = 0; k < dim; k++)
            {
                M[i*dim+k] -= (float80)f*LOG_FLOAT_RADIX;
                M[k*dim+i] += (float80)f*LOG_FLOAT_RADIX;
            }

            float80 min,max;
            matrix_float80_minmax(A,&min,&max);
            if(min > -4000 && max < 4000)
            {
                if(minimum != NULL)
                    *minimum = min;
                if(maximum != NULL)
                    *maximum = max;

                return;
            }
        }
    }
}

/* calculate log(det(1-M)) */
double matrix_logdet1mM(casimir_t *casimir, matrix_float80 *M, matrix_sign_t *M_sign)
{
    const size_t dim = M->dim;
    const char *detalg = casimir->detalg;
    double t = now();
    float80 minimum, maximum;
    int h,m,s;

    sec2human(t-casimir->birthtime, &h, &m, &s);
    casimir_printf(casimir, 2, "# calculating matrix elements: %02d:%02d:%02d\n", h,m,s);


    matrix_float80_minmax(M,&minimum,&maximum);
    casimir_printf(casimir, 2, "# before balancing: min=%Lg, max=%Lg\n", h,m,s, minimum, maximum);

    /* balance matrix */
    matrix_float80_log_balance(M, &minimum, &maximum);
    sec2human(now()-t, &h, &m, &s);
    casimir_printf(casimir, 2, "# balancing: %02d:%02d:%02d (min=%Lg, max=%Lg)\n", h,m,s, minimum, maximum);

    if(strcasecmp(detalg, "LU_FLOAT80") == 0)
    {
        /* exponentiate */
        matrix_float80_exp(M, M_sign);

        /* add unity matrix */
        for(size_t i = 0; i < dim; i++)
            M->M[i*dim+i] += 1;

        /* calculate log(det(M)) */
        t = now();
        const double logdet = matrix_float80_logdet_lu(M);
        TERMINATE(logdet > 0, "logdet > 0: %g", logdet);

        sec2human(now()-t, &h, &m, &s);
        casimir_printf(casimir, 2, "# QR-decomposition: %02d:%02d:%02d\n", h,m,s);
        casimir_printf(casimir, 2, "#\n");

        return logdet;
    }
    #ifdef FLOAT128
    else if(strcasecmp(detalg, "QR_FLOAT128") == 0)
    {
        size_t dim2 = pow_2(dim);
        matrix_float128 *M128 = matrix_float128_alloc(dim);

        for(size_t i = 0; i < dim2; i++)
            M128->M[i] = M_sign->M[i]*exp128(M->M[i]);

        for(size_t i = 0; i < dim; i++)
            M128->M[i*dim+i] += 1;

        const double logdet = matrix_float128_logdet_qr(M128);
        TERMINATE(logdet > 0, "logdet > 0: %g", logdet);

        matrix_float128_free(M128);

        sec2human(now()-t, &h, &m, &s);
        casimir_printf(casimir, 2, "# QR-decomposition: %02d:%02d:%02d\n", h,m,s);
        casimir_printf(casimir, 2, "#\n");

        return logdet;
    }
    #endif
    else
    {
        float80 *A = M->M;
        if(strcasecmp(detalg, "QR_FLOAT80") != 0)
            WARN(1, "Algorithm \"%s\" not supported. Defaulting to QR_FLOAT80.", detalg);

        matrix_float80_exp(M, M_sign);

        float80 trace  = 0;
        float80 trace2 = 0;
        for(size_t i = 0; i < dim; i++)
        {
            for(size_t k = 0; k < dim; k++)
                trace2 += A[i*dim+k]*A[k*dim+i];

            trace += A[i*dim+i];
        }
        casimir_printf(casimir, 2, "# Mercator (1): %Lg\n", +trace);
        casimir_printf(casimir, 2, "# Mercator (2): %Lg\n", +trace-trace2/2);

        /* add identity matrix */
        for(size_t i = 0; i < dim; i++)
            M->M[i*dim+i] += 1;

        /* pivot */
        if(casimir->pivot)
            matrix_float80_pivot(M);

        t = now();
        const double logdet = matrix_float80_logdet_qr(M);
        WARN(logdet > trace-trace2/2, "value of logdet > truncated Mercator series: logdet=%g, Mercator (2): %Lg", logdet, trace-trace2/2);
        TERMINATE(logdet > 0, "logdet > 0: %g", logdet);

        sec2human(now()-t, &h, &m, &s);
        casimir_printf(casimir, 2, "# QR-decomposition: %02d:%02d:%02d\n", h,m,s);
        casimir_printf(casimir, 2, "#\n");

        return logdet;
    }
}
