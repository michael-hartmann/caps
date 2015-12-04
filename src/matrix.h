#ifndef __MATRIX_H
#define __MATRIX_H

#include <stdio.h>
#include <math.h>

#include "edouble.h"
#include "libcasimir.h"
#include "utils.h"

#define FLOAT_RADIX       2.0
#define FLOAT_RADIX_SQ    (FLOAT_RADIX * FLOAT_RADIX)
#define LOG_FLOAT_RADIX   M_LN2
#define LOG_FLOAT_RADIX_SQ 2*M_LN2
#define LOG_095 -0.05129329438755058

/* LAPACK LU decomposition */
int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
/* LAPACK invert */
int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
/* matrix multiplication */
int dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c__, int *ldc);

#define MATRIX_TYPEDEF(NAME, MATRIX_TYPE) \
    typedef struct { \
        size_t size; \
        MATRIX_TYPE *M; \
    } NAME


MATRIX_TYPEDEF(matrix_sign_t, sign_t);
MATRIX_TYPEDEF(matrix_t, double);
MATRIX_TYPEDEF(matrix_edouble_t, edouble);

#define MATRIX_ALLOC(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) \
    MATRIX_TYPE *FUNCTION_PREFIX ## _alloc(size_t size)  \
    { \
        MATRIX_TYPE *matrix = xmalloc(sizeof(MATRIX_TYPE)); \
        if(matrix == NULL) \
            return NULL; \
 \
        matrix->size = size; \
        matrix->M = xmalloc(size*size*sizeof(TYPE)); \
        if(matrix->M == NULL) \
        { \
            FUNCTION_PREFIX ## _free(matrix); \
            return NULL; \
        } \
 \
        return matrix; \
    }

#define MATRIX_ALLOC_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) MATRIX_TYPE *FUNCTION_PREFIX ## _alloc(size_t size)

#define MATRIX_FREE(FUNCTION_PREFIX, MATRIX_TYPE) \
    void FUNCTION_PREFIX ## _free(MATRIX_TYPE *m) \
    { \
        if(m->M != NULL) \
        { \
            xfree(m->M); \
            m->M = NULL; \
        } \
        xfree(m); \
    }

#define MATRIX_FREE_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) void FUNCTION_PREFIX ## _free(MATRIX_TYPE *m)

#define MATRIX_LOGDET_LU(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION, LOG_FUNCTION) \
    TYPE FUNCTION_PREFIX ## _logdet_lu(MATRIX_TYPE *M) \
    { \
        const int dim = M->size; \
        int i,j,k; \
        TYPE sum, det = 0; \
        TYPE *a = M->M; \
\
        for(j = 0; j < dim; j++) \
        { \
            for(i = 0; i < j+1; i++) \
            { \
                sum = 0; \
                for(k = 0; k < i; k++) \
                    sum += a[i*dim+k]*a[k*dim+j]; \
                a[i*dim+j] -= sum; \
            } \
            for(i = j+1; i < dim; i++) \
            { \
                sum = 0; \
                for(k = 0; k < j; k++) \
                    sum += a[i*dim+k]*a[k*dim+j]; \
                a[i*dim+j] = (a[i*dim+j]-sum)/a[j*dim+j]; \
            } \
        } \
\
        for(i = 0; i < dim; i++) \
            det += LOG_FUNCTION(ABS_FUNCTION(a[i*dim+i])); \
        return det; \
    }

#define MATRIX_LOGDET_QR(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION, COPYSIGN_FUNCTION, SQRT_FUNCTION, LOG_FUNCTION) \
    TYPE FUNCTION_PREFIX ## _logdet_qr(MATRIX_TYPE *M) \
    { \
        size_t i, j, n, dim = M->size; \
        TYPE det = 0; \
        TYPE *m = M->M; \
\
        for(j = 0; j < dim-1; j++) \
            for(i = j+1; i < dim; i++) \
            {\
                TYPE c,s, Mij = m[i*dim+j]; \
\
                if(Mij != 0) \
                { \
                    const TYPE a = m[j*dim+j]; \
                    const TYPE b = Mij; \
\
                    if(b == 0) \
                    { \
                        c = COPYSIGN_FUNCTION(1,a); \
                        s = 0; \
                    } \
                    else if(a == 0) \
                    { \
                        c = 0; \
                        s = -COPYSIGN_FUNCTION(1, b); \
                    } \
                    else if(ABS_FUNCTION(b) > ABS_FUNCTION(a)) \
                    { \
                        const TYPE t = a/b; \
                        const TYPE u = COPYSIGN_FUNCTION(SQRT_FUNCTION(1+t*t),b); \
                        s = -1/u; \
                        c = -s*t; \
                    } \
                    else \
                    { \
                        const TYPE t = b/a; \
                        const TYPE u = COPYSIGN_FUNCTION(SQRT_FUNCTION(1+t*t),a); \
                        c = 1/u; \
                        s = -c*t; \
                    } \
 \
                    for(n = 0; n < dim; n++) \
                    { \
                        const TYPE Min = m[i*dim+n]; \
                        const TYPE Mjn = m[j*dim+n]; \
 \
                        m[i*dim+n] = c*Min + s*Mjn; \
                        m[j*dim+n] = c*Mjn - s*Min; \
                    } \
 \
                    /* m[i*dim+j] = 0; */ \
                } \
            } \
 \
        for(i = 0; i < dim; i++) \
            det += LOG_FUNCTION(ABS_FUNCTION(m[i*dim+i])); \
        return det; \
    }

#define MATRIX_LOGDET_QR_HEADER(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) TYPE FUNCTION_PREFIX ## _logdet_qr(MATRIX_TYPE *M)

#define MATRIX_LOGDET_LU_HEADER(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) TYPE FUNCTION_PREFIX ## _logdet_lu(MATRIX_TYPE *M)

#define MATRIX_ABSMAX(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION) \
    TYPE FUNCTION_PREFIX ## _absmax(MATRIX_TYPE *M) \
    { \
        size_t i,j, dim = M->size; \
        TYPE max = ABS_FUNCTION(matrix_get(M, 0,0)); \
 \
        for(i = 0; i < dim; i++) \
            for(j = 0; j < dim; j++) \
                max = MAX(max, ABS_FUNCTION(matrix_get(M, i,j))); \
 \
        return max; \
    }

#define MATRIX_ABSMAX_HEADER(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) TYPE FUNCTION_PREFIX ## _absmax(MATRIX_TYPE *M) \

#define MATRIX_ABSMIN(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION) \
    TYPE FUNCTION_PREFIX ## _absmin(MATRIX_TYPE *M) \
    { \
        size_t i,j, dim = M->size; \
        TYPE min = ABS_FUNCTION(matrix_get(M, 0,0)); \
 \
        for(i = 0; i < dim; i++) \
            for(j = 0; j < dim; j++) \
                min = MIN(min, ABS_FUNCTION(matrix_get(M, i,j))); \
 \
        return min; \
    }

#define MATRIX_ABSMIN_HEADER(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) TYPE FUNCTION_PREFIX ## _absmin(MATRIX_TYPE *M) \


#define MATRIX_EXP(FUNCTION_PREFIX, MATRIX_TYPE, EXP_FUNCTION) \
    void FUNCTION_PREFIX ## _exp(MATRIX_TYPE *M, matrix_sign_t *M_sign) \
    { \
        size_t i,j, dim = M->size; \
 \
        for(i = 0; i < dim; i++) \
            for(j = 0; j < dim; j++) \
                matrix_set(M, i,j, matrix_get(M_sign,i,j)*EXP_FUNCTION(matrix_get(M,i,j))); \
    }

#define MATRIX_EXP_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) void FUNCTION_PREFIX ## _exp(MATRIX_TYPE *M, matrix_sign_t *M_sign) \


/* Balance a general matrix by scaling the rows and columns, so the
 * new row and column norms are the same order of magnitude.
 *
 * B =  D^-1 A D
 *
 * where D is a diagonal matrix
 * 
 * This is necessary for the unsymmetric eigenvalue problem since the
 * calculation can become numerically unstable for unbalanced
 * matrices.  
 *
 * See Golub & Van Loan, "Matrix Computations" (3rd ed), Section 7.5.7
 * and Wilkinson & Reinsch, "Handbook for Automatic Computation", II/11 p320.
 */
#define MATRIX_BALANCE(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION) \
    void FUNCTION_PREFIX ## _balance(MATRIX_TYPE *A) \
    { \
        size_t i,j; \
        const size_t N = A->size; \
        int not_converged = 1; \
 \
        while(not_converged) \
        { \
            TYPE g, f, s; \
            TYPE row_norm, col_norm; \
 \
            not_converged = 0; \
 \
            for (i = 0; i < N; ++i) \
            { \
                row_norm = 0; \
                col_norm = 0; \
 \
                for (j = 0; j < N; ++j) \
                    if (j != i) \
                    { \
                      col_norm += ABS_FUNCTION(matrix_get(A, j, i)); \
                      row_norm += ABS_FUNCTION(matrix_get(A, i, j)); \
                    } \
 \
                if ((col_norm == 0.0) || (row_norm == 0.0)) \
                  continue; \
 \
                g = row_norm / FLOAT_RADIX; \
                f = 1.0; \
                s = col_norm + row_norm; \
 \
                /* \
                 * find the integer power of the machine radix which \
                 * comes closest to balancing the matrix \
                 */ \
                while (col_norm < g) \
                { \
                    f *= FLOAT_RADIX; \
                    col_norm *= FLOAT_RADIX_SQ; \
                } \
 \
                g = row_norm * FLOAT_RADIX; \
 \
                while (col_norm > g) \
                { \
                    f /= FLOAT_RADIX; \
                    col_norm /= FLOAT_RADIX_SQ; \
                } \
 \
                if ((row_norm + col_norm) < 0.95 * s * f) \
                { \
                    int k; \
                    not_converged = 1; \
 \
                    g = 1.0 / f; \
 \
                    /* \
                     * apply similarity transformation D, where \
                     * D_{ij} = f_i * delta_{ij} \
                     */ \
 \
                    /* multiply by D^{-1} on the left */ \
                    for(k = 0; k < N; k++) \
                        matrix_set(A, i,k, g*matrix_get(A,i,k)); \
 \
 \
                    /* multiply by D on the right */ \
                    for(k = 0; k < N; k++) \
                        matrix_set(A, k,i, f*matrix_get(A,k,i)); \
                } \
            } \
        } \
    }

#define MATRIX_BALANCE_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) void FUNCTION_PREFIX ## _balance(MATRIX_TYPE *A)


#define MATRIX_BALANCE_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) void FUNCTION_PREFIX ## _balance(MATRIX_TYPE *A)

/* We implement algorithm 3 (Balancing, proposed) from [1]. This algorithm
 * balances a square matrix A, so that on output A is nearly balanced in the
 * 1-norm.
 *
 * Notes:
 *  - The matrix is stored as the logarithm of the actual matrix elements.
 *  - We use the 1-norm.
 *  - We don't calculate D.
 *
 * [1] Rodney James, Julien Langou, Bradley R. Lowery, On matrix balancing and
 * eigenvector computation
 */
#define MATRIX_LOG_BALANCE(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, LOG_FUNCTION) \
    void FUNCTION_PREFIX ## _log_balance(MATRIX_TYPE *A) \
    { \
        size_t i,j; \
        const size_t N = A->size; \
        int converged = 0; \
\
        TYPE *M = A->M; \
        double *list_row    = xmalloc(N*sizeof(double)); \
        double *list_column = xmalloc(N*sizeof(double)); \
 \
        /* line 2 */ \
        while(!converged) \
        { \
            double f, s, col_max, row_max; \
            double row_norm, col_norm; \
\
            /* line 4 */ \
            converged = 1; \
\
            /* line 5 */ \
            for(i = 0; i < N; i++) \
            { \
                /* line 6 */ \
                col_max = row_max = list_column[0] = list_row[0] = M[0]; \
                for(j = 1; j < N; j++) \
                { \
                    const double Aji = M[j*N+i]; \
                    const double Aij = M[i*N+j]; \
                    list_column[j] = Aji; \
                    list_row[j]    = Aij; \
                    col_max = MAX(Aji,col_max); \
                    row_max = MAX(Aij,row_max); \
                } \
\
                /* faster than logadd_m */ \
                row_norm = col_norm = 0; \
                for(j = 0; j < N; j++) \
                { \
                    row_norm += exp(list_row[j]-row_max); \
                    col_norm += exp(list_column[j]-col_max); \
                } \
                row_norm = row_max + log(fabs(row_norm)); \
                col_norm = col_max + log(fabs(col_norm)); \
\
                if ((col_norm == LOG_FUNCTION(0)) || (row_norm == LOG_FUNCTION(0))) \
                    continue; \
\
                /* line 7 */ \
                f = 0; /* log(1)=0 */ \
                if(col_norm > row_norm) \
                    s = col_norm + log1p(exp(row_norm-col_norm)); \
                else \
                    s = row_norm + log1p(exp(col_norm-row_norm)); \
\
                /* line 8 */ \
                while(col_norm < (row_norm-LOG_FLOAT_RADIX)) \
                { \
                    /* line 9 */ \
                    col_norm += LOG_FLOAT_RADIX; \
                    row_norm -= LOG_FLOAT_RADIX; \
                    f        += LOG_FLOAT_RADIX; \
                } \
\
                /* line 10 */ \
                while(col_norm >= (row_norm+LOG_FLOAT_RADIX)) \
                { \
                    /* line 11 */ \
                    col_norm -= LOG_FLOAT_RADIX; \
                    row_norm += LOG_FLOAT_RADIX; \
                    f        -= LOG_FLOAT_RADIX; \
                } \
\
                /* line 12 */ \
                if(logadd(row_norm, col_norm) < (LOG_095+s)) \
                { \
                    int k; \
                    /* line 13 */ \
                    converged = 0; \
 \
                    /* line 14 */ \
                    for(k = 0; k < N; k++) \
                    { \
                        M[i*N+k] -= f; \
                        M[k*N+i] += f; \
                    } \
                } \
            } \
        } \
\
        xfree(list_column); \
        xfree(list_row); \
    }

#define MATRIX_LOG_BALANCE_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) void FUNCTION_PREFIX ## _log_balance(MATRIX_TYPE *A)


#define matrix_get(m, i, j)   ((m)->M[(i)*m->size+(j)])
#define matrix_set(m, i, j,v) ((m)->M[(i)*m->size+(j)]=(v))

MATRIX_ALLOC_HEADER(matrix_edouble, matrix_edouble_t);
MATRIX_FREE_HEADER (matrix_edouble, matrix_edouble_t);
MATRIX_ABSMIN_HEADER    (matrix_edouble, matrix_edouble_t, edouble);
MATRIX_ABSMAX_HEADER    (matrix_edouble, matrix_edouble_t, edouble);
MATRIX_BALANCE_HEADER   (matrix_edouble, matrix_edouble_t);
MATRIX_LOG_BALANCE_HEADER(matrix_edouble, matrix_edouble_t);
MATRIX_EXP_HEADER(matrix_edouble, matrix_edouble_t);

MATRIX_LOGDET_QR_HEADER(matrix_edouble, matrix_edouble_t, edouble);
MATRIX_LOGDET_LU_HEADER(matrix_edouble, matrix_edouble_t, edouble);

MATRIX_ALLOC_HEADER(matrix_sign, matrix_sign_t);
MATRIX_FREE_HEADER (matrix_sign, matrix_sign_t);

double matrix_edouble_logdet(matrix_edouble_t *M, matrix_sign_t *M_sign, const char *type);
double matrix_logdet_lu_lapack(matrix_edouble_t *M, matrix_sign_t *signs);
double matrix_logdet_block_lapack(matrix_edouble_t *M, matrix_sign_t *signs);

#endif
