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
 * pivoting \f$|M_{00}| < |M_{11}| < |M_{22}| \dots \f$.
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
            if(fabs80(Mzz) > fabs80(elem))
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
/**
 * @brief calculate QR decomposition of matrix
 *
 * This function calculates the QR decomposition of a matrix M and \f$\log|\det(M)|\f$.
 *
 * @param [in,out] M matrix
 * @retval logdet \f$\log|\det M|\f$
 */
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


/** @brief Precondition matrix
 *
 * This function preconditions the matrix A. The scaling of the matrix elements
 * ~b^(l1-l2) is estimated using linear regression of the anti diagonal of
 * M_EE. Then a similarity transform is performed.
 *
 * @param [in,out] A matrix
 */
void matrix_precondition(matrix_float80 *A)
{
    double b;
    const int dim = A->dim;
    const int dimby2 = dim/2;

    /* do linear regression => y(x) = a + b*x */
    {
        double y[dimby2];
        const double xm = (dimby2-1)/2.; /* sum formula: (0+1+2+...)/(dimby-1) **/
        double ym = 0;

        for(int i = 0; i < dimby2; i++)
        {
            y[i] = matrix_get(A, dimby2-1-i, i);
            ym += y[i];
        }

        ym /= dimby2;

        double num = 0, denom = 0;
        for(int i = 0; i < dimby2; i++)
        {
            num   += (i-xm)*(y[i]-ym);
            denom += pow_2(i-xm);
        }

        b = num/denom;
        /* a = ym-b*xm; */
    }

    for(int i = 0; i < dimby2; i++)
        for(int j = 0; j < dimby2; j++)
        {
            const float80 scale = 2*(i-j)*log80(b); /* ??? */
            float80 elem;

            /* EE */
            elem = matrix_get(A, i,j);
            matrix_set(A, i,j, elem+scale);

            /* EM */
            elem = matrix_get(A, i+dimby2,j);
            matrix_set(A, i+dimby2,j, elem+scale);

            /* ME */
            elem = matrix_get(A, i,j+dimby2);
            matrix_set(A, i,j+dimby2, elem+scale);

            /* MM */
            elem = matrix_get(A, i+dimby2,j+dimby2);
            matrix_set(A, i+dimby2,j+dimby2, elem+scale);
        }
}

void matrix_float80_balance(matrix_float80 *A)
{
    const float80 beta = 2;
    const int dim = A->dim;
    bool converged = false;
    float80 *M = A->M;

    while(!converged)
    {
        converged = true;

        for(int i = 0; i < dim; i++)
        {
            float80 c = 0, r = 0;

            for(int j = 0; j < dim; j++)
            {
                c += fabs80(matrix_get(A,j,i));
                r += fabs80(matrix_get(A,i,j));
            }

            const float80 s = c+r;
            float80 f = 1;

            while(c < r/beta)
            {
                c *= beta;
                r /= beta;
                f *= beta;
            }

            while(c >= beta*r)
            {
                c /= beta;
                r *= beta;
                f /= beta;
            }

            if((c+r) < BALANCE_STOP*s)
            {
                converged = false;

                /* line 14 */
                for(int k = 0; k < dim; k++)
                {
                    M[i*dim+k] /= f;
                    M[k*dim+i] *= f;
                }
            }
        }
    }
}

/** @brief Balance matrix A
 *
 * Balance a matrix that elements are give by logarithms.
 *
 * If minimum/maximum is not NULL, the minimum/maximum of A after balancing.
 *
 * @param [in,out] A matrix (elements given as logarithms)
 * @param [out]    minimum after balancing smallest matrix element (logarithm)
 * @param [out]    maximum after balancing largest matrix element (logarithm)
 */
void matrix_float80_log_balance(matrix_float80 *A)
{
    const double stop = log(0.95);
    const int dim = A->dim;
    bool converged = false;

    float80 *M = A->M;
    float80 *list_row = xmalloc(dim*sizeof(float80));
    float80 *list_col = xmalloc(dim*sizeof(float80));

    /* line 2 */
    while(!converged)
    {
        float80 minimum,maximum;
        matrix_float80_minmax(A,&minimum,&maximum);
        if(minimum > -6000 && maximum < 6000)
            break;

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


/**
 * @brief Calculate \f$\log\det(\mathrm{Id}+M)\f$ for matrix M
 *
 * This function calculates \f$\log\det(\mathrm{Id}+M)\f$ for a matrix \f$M\f$
 * which elements are in logarithmic representation and which signs are stored
 * in M_sign. A matrix element of \f$M\f$ is given by \f$\tilde M_{ij} =
 * (\mathrm{M\_sign})_{ij} \exp(M_{ij})\f$
 *
 * \f$\log\det(\mathrm{Id}+M)\f$ can also be calculated using the Mercator series:
 * \f$\log\det(\mathrm{Id}+M) \approx \mathrm{trace}(M) - \mathrm{trace}(M/2) + \dots\f$
 * We compute the first two terms of this series and check if the truncated
 * Mercator series \f$\mathrm{mercator2} < \log\det(\mathrm{Id}+M)\f$. If this
 * is not true a warning will be printed.
 *
 * If \f$\log\det(\mathrm{Id}+M) > 0\f$ the program is terminated and an error
 * is printed to stderr.
 *
 * As \f$M\f$ is usually a bad conditioned matrix, the calculation is performed
 * in the following way (in this process \f$M\f$ will be overwritten):
 *
 * - 1) Preconditioning
 *
 *   Preconditioning will drastically reduce the orders of magnitudes of the
 *   largest and smallest matrix elements. See function \ref
 *   matrix_precondition.  (This step will be performed if casimir->precondtion
 *   is true.)
 *
 * - 2) Balancing of logarithmic elements
 *
 *   This will further reduce the orders of
 *   magnitudes of largest and smallest matrix elements. Balancing is stopped as
 *   soon as it is possible to exponentiate all matrix elements without loss of
 *   significance. (This step will be performed if casimir->balance is true.)
 *
 * - 3) Exponentiating
 *
 *   All matrix elements are exponentiated and the signs of M_signs are
 *   multiplied: \f$M_{ij} = (\mathrm{M\_sign})_{ij} \exp{(M_{ij})} \f$
 *
 * - 4) Adding identity matrix
 *
 *   Add identify matrix: \f$M_{ij} = \mathrm{Id} + M_{ij}\f$
 *
 * - 5) Balancing
 *
 *   We further balance the matrix to make the QR decomposition more stable.
 *   But now we can operate on a matrix in "ordinary representation" and the
 *   algorithm is much faster than in step 2). (This step will be performed if
 *   casimir->balance is true.)
 *
 * - 6) Pivoting
 *
 *   See \ref matrix_float80_pivot. (The step will be performed if
 *   casimir->pivot is true.)
 *
 * - 7) QR decomposition
 *
 *   We use a series of Givens rotations to perform a QR decomposition: \f$M =
 *   QR\f$. Note that we only calculate \f$R\f$ and \f$M\f$ will be overwritten: \f$M=R\f$
 *
 * - 8) Calculating \f$\log\det M\f$
 *
 *   We calcalculate \f$\log\det(QR) = \log\det(R) = \log\det(M) = \sum_i
 *   \log|M_{ii}|\f$ (We assume that the determinant is positive and therefore
 *   \f$\log\det M\f$ is a real number.)
 *.
 *
 * This function will also print debugging information to stderr if
 * casimir->debug is set to true.
 *
 * @param [in]     casimir casimir object
 * @param [in,out] M round trip matrix M (matrix elements given as logarithms); M will be overwritten.
 * @param [in]     M_sign signs of matrix elements M
 * @retval logdet  \f$\log\det(\mathrm{Id}+M)\f$
 */
double matrix_logdetIdpM(casimir_t *casimir, matrix_float80 *M, matrix_sign_t *M_sign)
{
    const size_t dim = M->dim;
    const bool debug   = casimir->debug;
    const char *detalg = casimir->detalg;
    double logdetD;

    if(casimir->precondition)
    {
        const double start = now();
        matrix_precondition(M);

        if(debug)
        {
            float80 minimum, maximum;
            matrix_float80_minmax(M,&minimum,&maximum);
            casimir_debug(casimir, "# precondition: %gs (min=%Lg, max=%Lg)\n", now()-start, minimum, maximum);
        }
    }

    /* balance matrix */
    if(casimir->balance)
    {
        const double start = now();
        matrix_float80_log_balance(M);

        if(debug)
        {
            float80 minimum, maximum;
            matrix_float80_minmax(M,&minimum,&maximum);
            casimir_debug(casimir, "# log balancing: %gs (min=%Lg, max=%Lg)\n", now()-start, minimum, maximum);
        }
    }

    if(strcasecmp(detalg, "LU_FLOAT80") == 0)
    {
        /* exponentiate */
        matrix_float80_exp(M, M_sign);

        /* add unity matrix */
        for(size_t i = 0; i < dim; i++)
            M->M[i*dim+i] += 1;

        /* calculate log(det(M)) */
        logdetD = matrix_float80_logdet_lu(M);
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

        logdetD = matrix_float128_logdet_qr(M128);

        matrix_float128_free(M128);

        return logdet;
    }
    #endif
    else
    {
        float80 *A = M->M;
        if(strcasecmp(detalg, "QR_FLOAT80") != 0)
            WARN(1, "Algorithm \"%s\" not supported. Defaulting to QR_FLOAT80.", detalg);

        matrix_float80_exp(M, M_sign);

        float80 traceM  = 0; /* trace(M)  */
        float80 traceM2 = 0; /* trace(M²) */
        for(size_t i = 0; i < dim; i++)
        {
            for(size_t k = 0; k < dim; k++)
                traceM2 += A[i*dim+k]*A[k*dim+i];

            traceM += A[i*dim+i];
        }

        /* The mercator series for matrices is
         * log(Id+M) = M - M²/2 + ...
         *
         * Here:
         * log(det(Id+M)) = trace(log(Id+M)) = trace(M) - trace(M²)/2 + ...
         */
        const float80 mercator1 = traceM;
        const float80 mercator2 = traceM-traceM2/2;

        casimir_debug(casimir, "# Mercator: log(det(Id+M)) = trace(M)               = %.10Lg\n", mercator1);
        casimir_debug(casimir, "# Mercator: log(det(Id+M)) = trace(M) - trace(M²)/2 = %.10Lg\n", mercator2);

        /* balance */
        if(casimir->balance)
        {
            const double start = now();
            matrix_float80_balance(M);
            casimir_debug(casimir, "# balancing: %gs\n", now()-start);
        }

        /* add identity matrix */
        for(size_t i = 0; i < dim; i++)
            M->M[i*dim+i] += 1;

        /* pivot */
        if(casimir->pivot)
        {
            const double start = now();
            matrix_float80_pivot(M);
            casimir_debug(casimir, "# pivot: %gs\n", now()-start);
        }

        /* QR decomposition */
        {
            const double start = now();
            logdetD = matrix_float80_logdet_qr(M);
            casimir_debug(casimir, "# QR decomposition: %gs\n", now()-start);
        }

        WARN(logdetD > mercator2 && fabs80(logdetD-mercator2) > 1e-8, "value of logdet > truncated Mercator series: logdet=%.14g, Mercator (2): %.14Lg", logdetD, mercator2);
    }

    return logdetD;
}
