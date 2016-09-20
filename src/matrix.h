#ifndef __MATRIX_H
#define __MATRIX_H

#include <stdio.h>
#include <math.h>

#include "floattypes.h"
#include "utils.h"


/** define matrix type */
typedef struct {
    size_t dim,dim2;
    float64 *M;
} matrix_t;


/** macro to access matrix elements */
#define matrix_get(m, i, j)   ((m)->M[(i)*m->dim+(j)])

/** macro to set matrix elements */
#define matrix_set(m, i, j,v) ((m)->M[(i)*m->dim+(j)]=(v))


/* prototypes */
matrix_t *matrix_alloc(const size_t dim);
void matrix_free(matrix_t *A);

float64 matrix_logdet_lu(matrix_t *A);
float64 matrix_logdet_qr(matrix_t *M);
float64 matrix_logdet(matrix_t *A, float64 z, const char *detalg);

#endif
