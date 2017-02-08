#ifndef __MATRIX_H
#define __MATRIX_H

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "utils.h"

typedef enum { DETALG_LU, DETALG_QR, DETALG_EIG } detalg_t;

/** define matrix type */
typedef struct {
    size_t dim,dim2;
    size_t lda;
    double *M;
    bool free_memory;
} matrix_t;


/** macro to access matrix elements */
#define matrix_get(m, i, j)   ((m)->M[(i)*(m)->lda+(j)])

/** macro to set matrix elements */
#define matrix_set(m, i, j,v) ((m)->M[(i)*(m)->lda+(j)]=(v))


/* prototypes */
matrix_t *matrix_alloc(const size_t dim);
matrix_t *matrix_view(double *ptr, size_t dim, size_t lda);

void matrix_free(matrix_t *A);
void matrix_setall(matrix_t *A, double z);

int matrix_save_to_stream(matrix_t *A, FILE *stream);
int matrix_save_to_file(matrix_t *A, const char *filename);

double matrix_trace(matrix_t *A);
double matrix_trace2(matrix_t *A);

double matrix_norm_frobenius(matrix_t *A);

double matrix_logdet_triangular(matrix_t *A);
double matrix_logdet(matrix_t *A, double z, detalg_t detalg);
double matrix_logdet_lu(matrix_t *A);
double matrix_logdet_qr(matrix_t *A);
double matrix_logdetIdmM_eig(matrix_t *A, double z);

#endif
