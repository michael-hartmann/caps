/**
 * @file   matrix.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   November, 2016
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


/**
 * @brief Create new matrix object
 *
 * Create a new square matrix with dimension dim x dim. The matrix will not be
 * initialized.
 *
 * @param [in] dim dimension of square matrix
 * @retval A matrix
 */
matrix_t *matrix_alloc(const size_t dim)
{
    const size_t dim2 = dim*dim;

    if(dim == 0)
        return NULL;

    matrix_t *A = xmalloc(sizeof(matrix_t));
    if(A == NULL)
        return NULL;

    A->dim  = dim;
    A->dim2 = dim2;
    A->lda  = dim;

    A->M = xmalloc(dim2*sizeof(double));
    if(A->M == NULL)
    {
        matrix_free(A);
        return NULL;
    }

    A->free_memory = true;

    return A;
}

/**
 * @brief Create matrix view
 *
 * Create a matrix view from an existing matrix.
 *
 * Note that you still have to call matrix_free once you don't need the view
 * anymore. The actual data given by a will not be freed.
 *
 * @param [in] a double array with matrix data
 * @param [in] dim dimension of square matrix
 * @param [in] lda leading dimension of a
 * @retval A matrix view
 */
matrix_t *matrix_view(double *a, size_t dim, size_t lda)
{
    matrix_t *A = xmalloc(sizeof(matrix_t));
    if(A == NULL)
        return NULL;

    A->dim  = dim;
    A->dim2 = dim*dim;
    A->lda  = lda;
    A->M    = a;
    A->free_memory = false;

    return A;
}


/**
 * @brief Free matrix
 *
 * This function frees memory allocated for the matrix A.
 *
 * Note that you also have to call matrix_free on matrix views. see \ref
 * matrix_view.
 *
 * @param [in,out] A matrix
 */
void matrix_free(matrix_t *A)
{
    if(A != NULL)
    {
        if(A->free_memory)
        {
            xfree(A->M);
            A->M = NULL;
        }
        xfree(A);
    }
}

/**
 * @brief Save matrix to stream
 *
 * This function saves the matrix A to the stream given by stream. The output
 * is in the numpy .npy format.
 *
 * This function does not support matrix views at the moment.
 *
 * @param [in] A matrix
 * @param [in] stream stream
 * @retval 0
 */
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


/**
 * @brief Save matrix to file
 *
 * Save matrix A to file filename. See \ref matrix_save_to_stream for more
 * information.
 *
 * @param [in] A matrix
 * @param [in] filename filename of output file
 * @retval 0
 */
int matrix_save_to_file(matrix_t *A, const char *filename)
{
    FILE *f = fopen(filename, "w");
    if(f == NULL)
        return 1;

    int ret = matrix_save_to_stream(A, f);

    fclose(f);

    return ret;
}


/**
 * @brief Set all matrix elements to value z
 *
 * @param [in,out] A matrix
 * @param [in] z value
 */
void matrix_setall(matrix_t *A, double z)
{
    const size_t dim = A->dim;

    if(A->lda == dim)
    {
        double *a = A->M;
        const size_t dim2 = dim*dim;

        for(size_t i = 0; i < dim2; i++)
            a[i] = z;
    }
    else
    {
        for(size_t m = 0; m < dim; m++)
            for(size_t n = 0; n < dim; n++)
                matrix_set(A, m,n, z);
    }
}

/**
 * @brief Calculate trace of matrix
 * 
 * This function uses kahan sumation to decrease rounding errors.
 *
 * @param [in] A matrix
 * @retval trace trace of A
 */
double matrix_trace(matrix_t *A)
{
    const size_t dim = A->dim;
    double array[dim];

    for(size_t i = 0; i < dim; i++)
        array[i] = matrix_get(A, i,i);

    return kahan_sum(array, dim);
}


/**
 * @brief Calculate log(|det(A)|) for A triangular
 *
 * This function calculates the logarithm of the determinant of the matrix A
 * assuming A is upper or lower triangular.
 *
 * @param [in] A triangular matrix
 * @retval logdet log(|det(A)|)
 */
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
 *  - LU
 *  - QR
 *  - EIG
 *
 * @param [in,out] M round trip matrix M; M will be overwritten.
 * @param [in]     z factor z in log(det(Id+z*M))
 * @param [in]     detalg algorithm to be used
 * @retval logdet  \f$\log\det(\mathrm{Id}+z*M)\f$
 */
double matrix_logdet(matrix_t *A, double z, detalg_t detalg)
{
    const double trace = z*matrix_trace(A);

    /* log(det(Id+A)) ≈ tr(A) - 1/2 tr(A²) + ... */
    if(fabs(trace) < 1e-8)
        return trace;

    if(detalg == DETALG_EIG)
        return matrix_logdetIdmM_eig(A, z);

    /* M = Id+z*M */
    {
        const size_t dim  = A->dim;
        const size_t dim2 = A->dim2;
        double *a = A->M;

        /* multiply by z and add identity matrix */
        if(A->lda == dim)
        {
            /* dscal_(&dim2, &z, a, &one); // BLAS 1 */
            for(size_t i = 0; i < dim2; i++)
                a[i] *= z;

            for(size_t i = 0; i < dim; i++)
                a[i*dim+i] += 1;
        }
        else
        {
            for(size_t m = 0; m < dim; m++)
            {
                for(size_t n = 0; n < dim; n++)
                    matrix_set(A, m,n, z*matrix_get(A,m,n));

                matrix_set(A,m,m, 1+matrix_get(A,m,m));
            }
        }
    }

    if(detalg == DETALG_LU)
        return matrix_logdet_lu(A);
    if(detalg == DETALG_QR)
        return matrix_logdet_qr(A);

    WARN(1, "Unknown algorithm %d, defaulting to LU", detalg);
    return matrix_logdet_lu(A);
}

/**
 * @brief Calculate log(|det(A)|) using LU decomposition
 *
 * Calculate LU decomposition of A and use \ref matrix_logdet_triangular to
 * calculate log(|det(A)|).
 *
 * @param [in,out] A matrix
 * @retval logdet log(|det(A)|)
 */
double matrix_logdet_lu(matrix_t *A)
{
    int info = 0;
    int dim = (int)A->dim;
    int lda = (int)A->lda;
    int ipiv[dim];
    double *a = A->M;

    dgetrf_(
        &dim, /* M number of rows of A */
        &dim, /* N number of columns of A */
        a,    /* matrix A to be factored */
        &lda, /* LDA: leading dimension of A */
        ipiv, /* pivot indices of dimension min(rows,cols) */
        &info
    );

    WARN(info != 0, "dgetrf returned %d", info);

    return matrix_logdet_triangular(A);
}


/**
 * @brief Calculate log(|det(A)|) using QR decomposition
 *
 * Calculate QR decomposition of A and use \ref matrix_logdet_triangular to
 * calculate log(|det(A)|).
 *
 * @param [in,out] A matrix
 * @retval logdet log(|det(A)|)
 */
double matrix_logdet_qr(matrix_t *A)
{
    int dim = (int)A->dim;
    double tau[dim];
    int lda = A->lda;
    int info = 0;
    int lwork = -1;
    double opt_lwork = 0;
    double *work = NULL;

    /* determine optimal value for workspace */
    dgeqrf_(
        &dim,   /* number of rows of the matrix A */
        &dim,   /* number of columns of the matrix A */
        A->M,   /* matrix A */
        &lda,   /* leading dimension of the array A */
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
        &lda,   /* leading dimension of the array A */
        tau,    /* scalar factors of the elementary reflectors */
        work,   /* workspace */
        &lwork, /* dimension of the array WORK */
        &info
    );

    xfree(work);

    WARN(info != 0, "dgeqrf returned %d", info);

    return matrix_logdet_triangular(A);
}

/**
 * @brief Calculate log(|det(Id+z*A)|) using eigenvalues
 *
 * Calculate eigenvalues of A and calculate log(|det(Id+z*A)|).
 *
 * @param [in,out] A matrix
 * @param [in] z factor
 * @retval logdet log(|det(A)|)
 */
double matrix_logdetIdmM_eig(matrix_t *A, double z)
{
    int dim = A->dim;
    char jobvl = 'N'; /* don't compute left eigenvectors */
    char jobvr = 'N'; /* don't compute right eigenvectors */
    double wr[dim]; /* real parts of eigenvalues */
    double wi[dim]; /* imaginary parts of eigenvalues */
    int lda = (int)A->lda;
    int lwork = 100*dim;
    double logdet[dim];
    double *work = xmalloc(lwork*sizeof(double));
    int info = 0;

    dgeev_(
        &jobvl, /* JOBVL */
        &jobvr, /* JOBVR */
        &dim,   /* N */
        A->M,   /* A */
        &lda,   /* LDA */
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
