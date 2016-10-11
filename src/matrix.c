/**
 * @file   matrix.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   September, 2016
 * @brief  matrix functions
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include "matrix.h"
#include "sfunc.h"
#include "utils.h"

#include "clapack.h"

matrix_t *matrix_alloc(const size_t dim)
{
    const size_t dim2 = dim*dim;

    if(dim == 0)
        return NULL;

    matrix_t *A = xmalloc(sizeof(matrix_t));
    if(A == NULL)
        return NULL;

    A->dim  = dim;
    A->dim2 = dim*dim;
    A->lda  = dim;

    A->M = xmalloc(dim2*sizeof(double));
    if(A->M == NULL)
    {
        matrix_free(A);
        return NULL;
    }

    return A;
}

void matrix_free(matrix_t *A)
{
    if(A != NULL)
    {
        xfree(A->M);
        A->M = NULL;
        xfree(A);
    }
}

int matrix_save_to_stream(matrix_t *A, FILE *stream)
{
    /* dump matrix */
    char d_str[512] = { 0 };
    uint16_t len = 0;
    const size_t dim = A->dim;

    /* write magic string, major number and minor number */
    fwrite("\x93NUMPY\x01\x00", sizeof(char), 8, stream);

    /* write length of header and header */
    snprintf(d_str, sizeof(d_str)/sizeof(d_str[0]), "{'descr': '<f8', 'fortran_order': True, 'shape': (%zu, %zu), }", dim, dim);

    len = strlen(d_str);

    fwrite(&len,  sizeof(len),  1,   stream);
    fwrite(d_str, sizeof(char), len, stream);

    /* write matrix */
    fwrite(A->M, sizeof(double), A->dim2, stream);

    return 0;
}

int matrix_save_to_file(matrix_t *A, const char *filename)
{
    FILE *f = fopen(filename, "w");
    if(f == NULL)
        return 1;

    int ret = matrix_save_to_stream(A, f);

    fclose(f);

    return ret;
}

void matrix_setall(matrix_t *A, double z)
{
    const size_t dim2 = A->dim2;
    double *a = A->M;

    for(size_t i = 0; i < dim2; i++)
        a[i] = z;
}

double matrix_logdet_triangular(matrix_t *A)
{
    size_t dim = A->dim;
    double logdet[dim];

    for(size_t i = 0; i < dim; i++)
        logdet[i] = log(fabs(matrix_get(A, i, i)));

    return kahan_sum(logdet, dim);
}

/**
 * @brief Calculate \f$\log\det(\mathrm{Id}+z*M)\f$ for matrix M
 *
 * Detalg may be:
 *  - LU_LAPACK (default)
 *  - QR_LAPACK
 *  - LU
 *  - QR_GIVENS
 *  - EIG_LAPACK
 *
 * @param [in,out] M round trip matrix M; M will be overwritten.
 * @param [in]     z factor z in log(det(Id+z*M))
 * @param [in]     detalg algorithm to be used
 * @retval logdet  \f$\log\det(\mathrm{Id}+z*M)\f$
 */
double matrix_logdet(matrix_t *A, double z, const char *detalg)
{
    if(strcmp(detalg, "EIG_LAPACK") == 0)
        return matrix_logdetIdmM_eig_lapack(A, z);

    /* M = Id+z*M */
    {
        const size_t dim  = A->dim;
        const size_t dim2 = A->dim2;
        double *a = A->M;

        /* multiply by z */
        /* dscal_(&dim2, &z, a, &one); // BLAS 1 */
        for(size_t i = 0; i < dim2; i++)
            a[i] *= z;

        /* add identity matrix */
        for(size_t i = 0; i < dim; i++)
            a[i*dim+i] += 1;
    }

    if(strcmp(detalg, "LU_LAPACK") == 0)
        return matrix_logdet_lu_lapack(A);
    if(strcmp(detalg, "QR_LAPACK") == 0)
        return matrix_logdet_qr_lapack(A);

    WARN(1, "Unknown algorithm %s, defaulting to LU_LAPACK", detalg);
    return matrix_logdet_lu_lapack(A);
}

double matrix_logdet_lu_lapack(matrix_t *A)
{
    int info = 0;
    int dim = (int)A->dim;
    int ipiv[dim];
    double *a = A->M;

    dgetrf_(
        &dim, /* M number of rows of A */
        &dim, /* N number of columns of A */
        a,    /* matrix A to be factored */
        &dim, /* LDA: leading dimension of A */
        ipiv, /* pivot indices of dimension min(rows,cols) */
        &info
    );

    WARN(info != 0, "dgetrf returned %d", info);

    return matrix_logdet_triangular(A);
}

double matrix_logdet_qr_lapack(matrix_t *A)
{
    int dim = (int)A->dim;
    double tau[dim];
    int info = 0;
    int lwork = -1;
    double opt_lwork = 0;
    double *work = NULL;

    /* determine optimal value for workspace */
    dgeqrf_(
        &dim,   /* number of rows of the matrix A */
        &dim,   /* number of columns of the matrix A */
        A->M,   /* matrix A */
        &dim,   /* leading dimension of the array A */
        tau,   /* scalar factors of the elementary reflectors */
        &opt_lwork,/* workspace */
        &lwork, /* dimension of the array WORK */
        &info
    );

    /* allocate memory for workspace */
    lwork = opt_lwork;
    work = xmalloc(lwork*sizeof(double));

    /* perform QR decomposition */
    dgeqrf_(
        &dim,   /* number of rows of the matrix A */
        &dim,   /* number of columns of the matrix A */
        A->M,   /* matrix A */
        &dim,   /* leading dimension of the array A */
        tau,    /* scalar factors of the elementary reflectors */
        work,   /* workspace */
        &lwork, /* dimension of the array WORK */
        &info
    );

    xfree(work);

    WARN(info != 0, "dgeqrf returned %d", info);

    return matrix_logdet_triangular(A);
}

double matrix_logdetIdmM_eig_lapack(matrix_t *A, double z)
{
    int dim = A->dim;
    char jobvl = 'N'; /* don't compute left eigenvectors */
    char jobvr = 'N'; /* don't compute right eigenvectors */
    double wr[dim]; /* real parts of eigenvalues */
    double wi[dim]; /* imaginary parts of eigenvalues */
    int lwork = 100*dim;
    double logdet[dim];
    double *work = xmalloc(lwork*sizeof(double));
    int info = 0;

    dgeev_(
        &jobvl, /* JOBVL */
        &jobvr, /* JOBVR */
        &dim,   /* N */
        A->M,   /* A */
        &dim,   /* LDA */
        wr,     /* WR */
        wi,     /* WI */
        NULL,   /* VL */
        &dim,   /* LDVL */
        NULL,   /* VR */
        &dim,   /* LDVR */
        work,   /* WORK */
        &lwork, /* LWORK */
        &info   /* INFO */
    );

    xfree(work);

    WARN(info != 0, "dgetrf returned %d", info);

    for(int i = 0; i < dim; i++)
    {
        double lambda = wr[i];

        if(lambda < -z)
            logdet[i] = log1p(+z*lambda);
        else
            logdet[i] = log(fabs(1-z*lambda));
    }

    return kahan_sum(logdet, dim);
}
