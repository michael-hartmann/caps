#ifndef __MATRIX_H
#define __MATRIX_H

#include <stdio.h>
#include <math.h>
#ifdef FLOAT128
#include <quadmath.h>
#endif

#include "floattypes.h"
#include "libcasimir.h"
#include "utils.h"


#define LOG_FLOAT_RADIX   0.6931471805599453094172321214581765680755001343602552L
#define LOG095 -0.05129329438755058

#define BALANCE_STOP 0.95

/** define matrix type */
#define MATRIX_TYPEDEF(NAME, TYPE) \
    typedef struct { \
        int dim; \
        TYPE *M; \
    } NAME


/** macro to access matrix elements */
#define matrix_get(m, i, j)   ((m)->M[(i)*m->dim+(j)])
/** macro to set matrix elements */
#define matrix_set(m, i, j,v) ((m)->M[(i)*m->dim+(j)]=(v))

/** define various matrix types */
MATRIX_TYPEDEF(matrix_float80, float80);
MATRIX_TYPEDEF(matrix_sign_t, sign_t);
#ifdef FLOAT128
MATRIX_TYPEDEF(matrix_float128, float128);
#endif

/** matrix allocation */
#define MATRIX_ALLOC(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) \
    MATRIX_TYPE *FUNCTION_PREFIX ## _alloc(int dim)  \
    { \
        if(dim <= 0) \
            return NULL; \
\
        MATRIX_TYPE *matrix = xmalloc(sizeof(MATRIX_TYPE)); \
        if(matrix == NULL) \
            return NULL; \
 \
        const int dim2 = (size_t)dim*(size_t)dim; \
        matrix->dim = dim; \
        matrix->M = xmalloc(dim2*sizeof(TYPE)); \
        if(matrix->M == NULL) \
        { \
            FUNCTION_PREFIX ## _free(matrix); \
            return NULL; \
        } \
 \
        return matrix; \
    }

#define MATRIX_ALLOC_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) MATRIX_TYPE *FUNCTION_PREFIX ## _alloc(int dim)

#define MATRIX_FREE(FUNCTION_PREFIX, MATRIX_TYPE) \
    void FUNCTION_PREFIX ## _free(MATRIX_TYPE *m) \
    { \
        if(m != NULL) \
        { \
            xfree(m->M); \
            m->M = NULL; \
            xfree(m); \
        } \
    }

#define MATRIX_FREE_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) void FUNCTION_PREFIX ## _free(MATRIX_TYPE *m)


#define MATRIX_MINMAX(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) \
    void FUNCTION_PREFIX ## _minmax(MATRIX_TYPE *M, TYPE *min, TYPE *max) \
    { \
        if(min == NULL && max == NULL) \
            return; \
\
        const int dim = M->dim; \
        const size_t dim2 = (size_t)dim*(size_t)dim; \
        TYPE minimum = M->M[0]; \
        TYPE maximum = M->M[0]; \
\
        for(size_t i = 1; i < dim2; i++) \
        { \
            const TYPE elem = M->M[i]; \
            if     (elem < minimum) minimum = elem; \
            else if(elem > maximum) maximum = elem; \
        } \
\
        if(min != NULL) *min = minimum; \
        if(max != NULL) *max = maximum; \
    }

#define MATRIX_MINMAX_HEADER(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) void FUNCTION_PREFIX ## _minmax(MATRIX_TYPE *M, TYPE *min, TYPE *max)

#define MATRIX_SAVE(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) \
    int FUNCTION_PREFIX ## _save(MATRIX_TYPE *M, const char *path) \
    { \
        FILE *f; \
        const int dim = M->dim; \
        const size_t dim2 = (size_t)dim*(size_t)dim; \
        const TYPE *ptr = M->M; \
\
        if((f = fopen(path, "w")) == NULL) \
            goto fail; \
\
        if(fwrite(&dim, sizeof(int), 1, f) != 1) \
            goto fail; \
\
        if(fwrite(ptr, sizeof(TYPE), dim2, f) != dim2) \
            goto fail; \
\
        if(fclose(f) == 0) \
            return 0; \
\
        fail: \
        if(f != NULL) \
            fclose(f); \
\
        return 1; \
    }

#define MATRIX_SAVE_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) int FUNCTION_PREFIX ## _save(MATRIX_TYPE *M, const char *path)

#define MATRIX_LOAD(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, MATRIX_ALLOC, MATRIX_FREE) \
    MATRIX_TYPE *FUNCTION_PREFIX ## _load(const char *path) \
    { \
        MATRIX_TYPE *M = NULL; \
        FILE *f; \
        int dim; \
\
        if((f = fopen(path, "r")) == NULL) \
            goto fail; \
\
        if(fread(&dim, sizeof(dim), 1, f) != 1) \
            goto fail; \
\
        if(dim <= 0) \
            goto fail; \
\
        const size_t dim2 = (size_t)dim*(size_t)dim; \
\
        M = MATRIX_ALLOC(dim); \
        if(M == NULL) \
            goto fail; \
\
        if(fread(M->M, sizeof(TYPE), dim2, f) != dim2) \
            goto fail; \
\
        if(fclose(f) == 0) \
            return M; \
\
        fail: \
        if(M != NULL) \
            MATRIX_FREE(M); \
        if(f != NULL) \
            fclose(f); \
\
        return NULL; \
    }

#define MATRIX_LOAD_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) MATRIX_TYPE *FUNCTION_PREFIX ## _load(const char *path)

#define MATRIX_LOGDET_LU(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION, LOG_FUNCTION) \
    TYPE FUNCTION_PREFIX ## _logdet_lu(MATRIX_TYPE *M) \
    { \
        const int dim = M->dim; \
        TYPE sum, det = 0; \
        TYPE *a = M->M; \
\
        for(int j = 0; j < dim; j++) \
        { \
            for(int i = 0; i < j+1; i++) \
            { \
                sum = 0; \
                for(int k = 0; k < i; k++) \
                    sum += a[i*dim+k]*a[k*dim+j]; \
                a[i*dim+j] -= sum; \
            } \
            for(int i = j+1; i < dim; i++) \
            { \
                sum = 0; \
                for(int k = 0; k < j; k++) \
                    sum += a[i*dim+k]*a[k*dim+j]; \
                a[i*dim+j] = (a[i*dim+j]-sum)/a[j*dim+j]; \
            } \
        } \
\
        for(int i = 0; i < dim; i++) \
            det += LOG_FUNCTION(ABS_FUNCTION(a[i*dim+i])); \
        return det; \
    }

#define MATRIX_LOGDET_QR(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION, COPYSIGN_FUNCTION, SQRT_FUNCTION, LOG_FUNCTION) \
    TYPE FUNCTION_PREFIX ## _logdet_qr(MATRIX_TYPE *M) \
    { \
        const int dim = M->size; \
        TYPE det = 0; \
        TYPE *m = M->M; \
\
        for(int j = 0; j < dim-1; j++) \
            for(int i = j+1; i < dim; i++) \
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
                    for(int n = 0; n < dim; n++) \
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
        for(int i = 0; i < dim; i++) \
            det += LOG_FUNCTION(ABS_FUNCTION(m[i*dim+i])); \
        return det; \
    }

//#define MATRIX_LOGDET_QR_HEADER(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) TYPE FUNCTION_PREFIX ## _logdet_qr(MATRIX_TYPE *M)

#define MATRIX_LOGDET_LU_HEADER(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) TYPE FUNCTION_PREFIX ## _logdet_lu(MATRIX_TYPE *M)

#define MATRIX_EXP(FUNCTION_PREFIX, MATRIX_TYPE, EXP_FUNCTION) \
    void FUNCTION_PREFIX ## _exp(MATRIX_TYPE *M, matrix_sign_t *M_sign) \
    { \
        const int dim = M->dim; \
        const size_t dim2 = (size_t)dim*(size_t)dim; \
\
        for(size_t i = 0; i < dim2; i++) \
            M->M[i] = M_sign->M[i]*EXP_FUNCTION(M->M[i]); \
    }

#define MATRIX_EXP_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) void FUNCTION_PREFIX ## _exp(MATRIX_TYPE *M, matrix_sign_t *M_sign) \

/** matrix functions for float80 */
MATRIX_ALLOC_HEADER (matrix_float80, matrix_float80);
MATRIX_FREE_HEADER  (matrix_float80, matrix_float80);
MATRIX_LOAD_HEADER  (matrix_float80, matrix_float80);
MATRIX_SAVE_HEADER  (matrix_float80, matrix_float80);
MATRIX_EXP_HEADER   (matrix_float80, matrix_float80);
MATRIX_MINMAX_HEADER(matrix_float80, matrix_float80, float80);
MATRIX_LOGDET_LU_HEADER(matrix_float80, matrix_float80, float80);

/** matrix functions for sign_t */
MATRIX_ALLOC_HEADER (matrix_sign, matrix_sign_t);
MATRIX_FREE_HEADER  (matrix_sign, matrix_sign_t);
MATRIX_LOAD_HEADER  (matrix_sign, matrix_sign_t);
MATRIX_SAVE_HEADER  (matrix_sign, matrix_sign_t);

/** matrix functions for float128 */
#ifdef FLOAT128
MATRIX_ALLOC_HEADER (matrix_float128, matrix_float128);
MATRIX_FREE_HEADER  (matrix_float128, matrix_float128);
MATRIX_LOAD_HEADER  (matrix_float128, matrix_float128);
MATRIX_SAVE_HEADER  (matrix_float128, matrix_float128);
MATRIX_MINMAX_HEADER(matrix_float128, matrix_float128, float128);
#endif

/* prototypes */
double matrix_logdetIdpM(casimir_t *casimir, matrix_float80 *M, matrix_sign_t *M_sign);

double matrix_float80_logdet_qr(matrix_float80 *M);
#ifdef FLOAT128
double matrix_float128_logdet_qr(matrix_float128 *M);
#endif

void matrix_precondition(matrix_float80 *A);

void matrix_float80_log_balance(matrix_float80 *A);
void matrix_float80_balance(matrix_float80 *A);

void matrix_float80_swap(matrix_float80 *M, const int i, const int j);
void matrix_float80_pivot(matrix_float80 *M);

#endif
