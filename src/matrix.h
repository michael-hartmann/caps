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
//#define LOG_FLOAT_RADIX 1L
#define LOG095 -0.05129329438755058

#define MATRIX_TYPEDEF(NAME, MATRIX_TYPE) \
    typedef struct { \
        int size; \
        MATRIX_TYPE *M; \
    } NAME


#define matrix_get(m, i, j)   ((m)->M[(i)*m->size+(j)])
#define matrix_set(m, i, j,v) ((m)->M[(i)*m->size+(j)]=(v))

MATRIX_TYPEDEF(matrix_sign_t, sign_t);
MATRIX_TYPEDEF(matrix_float80, float80);

#ifdef FLOAT128
MATRIX_TYPEDEF(matrix_float128, float128);
#endif

void matrix_float80_log_balance(matrix_float80 *A);

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


#define MATRIX_MINMAX(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) \
    void FUNCTION_PREFIX ## _minmax(MATRIX_TYPE *M, TYPE *min, TYPE *max) \
    { \
        const int dim = M->size; \
        TYPE minimum = matrix_get(M, 0,0); \
        TYPE maximum = matrix_get(M, 0,0); \
\
        for(int i = 0; i < dim; i++) \
            for(int j = 0; j < dim; j++) \
            { \
                const TYPE elem = matrix_get(M, i,j); \
                if(elem < minimum) minimum = elem; \
                if(elem > maximum) maximum = elem; \
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
        const int size = M->size; \
        const TYPE *ptr = M->M; \
\
        if((f = fopen(path, "w")) == NULL) \
            goto fail; \
\
        if(fwrite(&size, sizeof(int), 1, f) != 1) \
            goto fail; \
\
        if(fwrite(ptr,   sizeof(TYPE), pow_2(size), f) != (size_t)pow_2(size)) \
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
        int size; \
\
        if((f = fopen(path, "r")) == NULL) \
            goto fail; \
\
        if(fread(&size, sizeof(size), 1, f) != 1) \
            goto fail; \
\
        if(size <= 0) \
            goto fail; \
\
        M = MATRIX_ALLOC(size); \
        if(M == NULL) \
            goto fail; \
\
        if(fread(M->M, sizeof(TYPE), pow_2(size), f) != (size_t)pow_2(size)) \
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
        int dim = M->size; \
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
        const int dim = M->size; \
 \
        for(int i = 0; i < dim; i++) \
            for(int j = 0; j < dim; j++) \
                matrix_set(M, i,j, matrix_get(M_sign,i,j)*EXP_FUNCTION(matrix_get(M,i,j))); \
    }

#define MATRIX_EXP_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) void FUNCTION_PREFIX ## _exp(MATRIX_TYPE *M, matrix_sign_t *M_sign) \

MATRIX_ALLOC_HEADER (matrix_float80, matrix_float80);
MATRIX_FREE_HEADER  (matrix_float80, matrix_float80);
MATRIX_LOAD_HEADER  (matrix_float80, matrix_float80);
MATRIX_SAVE_HEADER  (matrix_float80, matrix_float80);
MATRIX_EXP_HEADER   (matrix_float80, matrix_float80);
MATRIX_MINMAX_HEADER(matrix_float80, matrix_float80, float80);

MATRIX_LOGDET_LU_HEADER(matrix_float80, matrix_float80, float80);

MATRIX_ALLOC_HEADER (matrix_sign, matrix_sign_t);
MATRIX_FREE_HEADER  (matrix_sign, matrix_sign_t);
MATRIX_LOAD_HEADER  (matrix_sign, matrix_sign_t);
MATRIX_SAVE_HEADER  (matrix_sign, matrix_sign_t);

/*
MATRIX_ALLOC_HEADER(matrix_sfloat, matrix_sfloat_t);
MATRIX_FREE_HEADER (matrix_sfloat, matrix_sfloat_t);
MATRIX_LOAD_HEADER (matrix_sfloat, matrix_sfloat_t);
MATRIX_SAVE_HEADER (matrix_sfloat, matrix_sfloat_t);
*/

#ifdef FLOAT128
MATRIX_ALLOC_HEADER (matrix_float128, matrix_float128);
MATRIX_FREE_HEADER  (matrix_float128, matrix_float128);
MATRIX_LOAD_HEADER  (matrix_float128, matrix_float128);
MATRIX_SAVE_HEADER  (matrix_float128, matrix_float128);
MATRIX_MINMAX_HEADER(matrix_float128, matrix_float128, float128);
#endif

double matrix_float80_logdet(matrix_float80 *M, matrix_sign_t *M_sign, const char *type);
double matrix_float80_logdet_qr(matrix_float80 *M);
#ifdef FLOAT128
double matrix_float128_logdet_qr(matrix_float128 *M);
#endif

//double matrix_floatdd_logdet(matrix_float80 *M, matrix_sign_t *M_sign, const char *type);
double matrix_float80_log_logdet_qr(matrix_float80 *M, matrix_sign_t *M_sign);

#endif
