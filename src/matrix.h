#ifndef __MATRIX_H
#define __MATRIX_H

#include <stdio.h>
#include <math.h>

#include "utils.h"


/** define matrix type */
typedef struct {
    size_t dim,dim2;
    double *M;
} matrix_t;


/** macro to access matrix elements */
#define matrix_get(m, i, j)   ((m)->M[(i)*m->dim+(j)])

/** macro to set matrix elements */
#define matrix_set(m, i, j,v) ((m)->M[(i)*m->dim+(j)]=(v))


/* prototypes */
matrix_t *matrix_alloc(const size_t dim);
void matrix_free(matrix_t *A);
void matrix_setall(matrix_t *A, double z);

int matrix_save_to_stream(matrix_t *A, FILE *stream);
int matrix_save_to_file(matrix_t *A, const char *filename);

double matrix_logdet_lu(matrix_t *A);
double matrix_logdet_qr(matrix_t *M);
double matrix_logdet(matrix_t *A, double z, const char *detalg);
double matrix_logdet_lu_lapack(matrix_t *A);
double matrix_logdetIdmM_eig_lapack(matrix_t *A, double z);

#endif
