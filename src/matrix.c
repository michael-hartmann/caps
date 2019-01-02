/**
 * @file   matrix.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   January, 2019
 * @brief  Matrix functions
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include <hodlr.h>

#include "matrix.h"
#include "misc.h"
#include "utils.h"

#include "clapack.h"


/** @brief Compute \f$\log \det(1-A)\f$
 *
 * This function computes \f$\log \det(1-A)\f$ using either the HODLR approach or
 * LU decomposition. The matrix \f$A\f$ is given as a callback function. This
 * callback accepts two integers, the row and the column of the matrix entry
 * (starting from 0), and a pointer to args. The callback returns the
 * corresponding matrix element.
 *
 * If the matrix elements of \f$A\f$ are small, i.e., if the modulus of the
 * trace is smaller than 1e-8, the trace will be used as an approximation to
 * prevent a loss of significance. If the modulus of the trace is larger than
 * the modulus of the value computed using HODLR, the trace approximation is
 * returned.
 *
 * If the determinant is not computed using the HODLR approach, all matrix
 * elements have to be computed. In this case the matrix \f$A\f$ is written to
 * the filesystem if the environment variable CASIMIR_DUMP is set. If the
 * variable is set, the matrix will be stored in the path given by CASIMIR_DUMP
 * as a two-dimensional numpy array (npy). This option might be useful for
 * debugging. Also note that if detalg is CHOLESKY, only the upper half of the
 * matrix will be initialized.
 *
 * @param [in] dim       dimension of matrix
 * @param [in] kernel    callback function that returns matrix elements of \f$A\f$
 * @param [in] args      pointer given to callback function kernel
 * @param [in] symmetric bool indicating whether matrix is symmetric
 * @param [in] detalg    algorithm (DETALG_HODLR, DETALG_LU, DETALG_QR, DETALG_CHOLESKY)
 * @retval logdet \f$\log \det(1-A)\f$
 */
double kernel_logdet(int dim, double (*kernel)(int,int,void *), void *args, int symmetric, detalg_t detalg)
{
    double logdet = NAN;
    double *diagonal = xmalloc(((size_t)(dim))*sizeof(double));

    /* calculate diagonal elements */
    for(int n = 0; n < dim; n++)
        diagonal[n] = kernel(n,n,args);

    const double trace = kahan_sum(diagonal, dim);

    /* use trace approximation to avoid cancellation
     *
     * log det(Id-M) = tr log(Id-M) = -tr(M + M²/2 + M³/3 + ...) = -tr(M) + R
     *
     * tr(M)   = Σ_i λ_i                where 0<λ_i<1 are the eigenvalues of M
     * tr(M^r) = Σ_i λ_i^r < tr(M)^r
     *
     * Therefore we find:
     * |R| = |tr(M²/2 + M³/3 + ...)| < |tr(M²+M³+...)|/2 < |tr(M)²+tr(M)³|/2 = |tr(M)²/(2-2tr(M))| =~ |tr(M)²|/2
     */
    if(fabs(trace) < 1e-8)
    {
        xfree(diagonal);
        return -trace;
    }

    if(detalg != DETALG_HODLR)
    {
        /* allocate space for matrix M */
        matrix_t *M = matrix_alloc(dim);

        /* set matrix elements to 0 */
        matrix_setall(M, 0);

        /* copy diagonal */
        for(size_t k = 0; k < (size_t)dim; k++)
            matrix_set(M, k,k, diagonal[k]);

        /* n-th minor diagonal */
        for(size_t md = 1; md < (size_t)dim; md++)
            for(size_t k = 0; k < (size_t)dim-md; k++)
            {
                /* for cholesky decomposition we only need the upper part of
                 * the matrix */
                if(detalg != DETALG_CHOLESKY)
                    matrix_set(M, k,md+k, kernel(k,md+k,args));
                matrix_set(M, md+k,k, kernel(md+k,k,args));
            }

        /* dump */
        const char *filename = getenv("CASIMIR_DUMP");
        if(filename != NULL)
            matrix_save_to_file(M, filename);

        /* compute logdet */
        logdet = matrix_logdet_dense(M, -1, detalg);

        matrix_free(M);
        xfree(diagonal);

        return logdet;
    }
    else
    {
        /* HODLR */

        /* nLeaf is the size (number of rows of the matrix) of the smallest
         * block at the leaf level. The number of levels in the tree is given
         * by n_levels=log_2(N/nLeaf) where N denotes the dimension of the
         * matrix.
         */
        const unsigned int nLeaf = 50;

        /* Choose relative error to compute the determinant as ~1e-13.
         * The value of the determinant is estimated using the trace. As
         *      |log det(Id-M)| < trace(M)
         * the estimate trace*1e-13 gives actually a lower error than 1e-13.
         */
        const double tolerance = fmax(1e-13, trace*1e-13);

        /* calculate log(det(D)) using HODLR approach */
        logdet = hodlr_logdet_diagonal(dim, kernel, args, diagonal, nLeaf, tolerance, symmetric);

        xfree(diagonal);

        /* if |trace| > |log(det(D))|, then the trace result is more accurate */
        if(fabs(trace) > fabs(logdet))
            return -trace;

        return logdet;
    }
}

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

    return A;
}

/**
 * @brief Free matrix
 *
 * This function frees the memory allocated for the matrix A.
 *
 * @param [in,out] A matrix
 */
void matrix_free(matrix_t *A)
{
    if(A != NULL)
    {
        xfree(A->M);
        xfree(A);
    }
}

/**
 * @brief Save matrix to stream
 *
 * This function saves the matrix \f$A\f$ to the stream given by stream. The
 * output is in the numpy .npy format.
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
 * Save matrix \f$A\f$ to file filename. See \ref matrix_save_to_stream for
 * more information.
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
 * @brief Load matrix from stream
 *
 * This function loads a matrix from a given stream. The input must be in .npy
 * format. The input matrix must be a square matrix.
 *
 * The function will rudimentary parse the description string and abort if an
 * error occures. Do not use this function on untrusted data.
 *
 * @param [in] stream stream
 * @retval A matrix if successful
 * @retval NULL if an error occured
 */
matrix_t *matrix_load_from_stream(FILE *stream)
{
    matrix_t *A = NULL;
    size_t dim1, dim2;
    char *p = NULL;
    char header[8] = { 0 };
    char d_str[512] = { 0 };
    uint16_t len = 0;
    size_t ret;

    /* check if header is correct */
    ret = fread(header, sizeof(char), 8, stream);
    if(ret != 8)
        return NULL;
    if(strcmp(header, "\x93NUMPY\x01\x00") != 0)
        return NULL;

    /* read len */
    ret = fread(&len, sizeof(len), 1, stream);
    if(ret != 1)
        return NULL;
    /* read description */
    ret = fread(d_str, sizeof(char), len, stream);
    if(ret != len)
        return NULL;

    if(len < 2)
        return NULL;

    if(d_str[0] != '{' && d_str[len-1] != '}')
        return NULL;

    if(strstr(d_str, "'descr': '<f8'") == NULL)
        return NULL;

    if(strstr(d_str, "'fortran_order': True") == NULL)
        return NULL;

    p = strstr(d_str, "'shape': (");
    if(p == NULL)
        return NULL;

    p += 10;
    dim1 = atoi(p);
    p = strstr(p, ",");
    if(p == NULL)
        return NULL;

    dim2 = atoi(p+1);

    if(dim1 != dim2)
        return NULL;

    A = matrix_alloc(dim1);
    ret = fread(A->M, sizeof(double), dim1*dim1, stream);

    if(ret != dim1*dim1)
        return NULL;

    return A;
}

/**
 * @brief Load matrix from file
 *
 * Load matrix matrix from file filename. See \ref matrix_load_from_stream for
 * more information.
 *
 * @param [in] filename filename of output file
 * @retval A matrix if successful
 * @retval NULL if an error occured
 */
matrix_t *matrix_load_from_file(const char *filename)
{
    FILE *f = fopen(filename, "r");
    if(f == NULL)
        return NULL;

    matrix_t *A = matrix_load_from_stream(f);

    fclose(f);

    return A;
}

/**
 * @brief Set all matrix elements to value \f$z\f$
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
 * This function uses Kahan sumation (see \ref kahan_sum) to reduce rounding
 * errors.
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
 * @brief Calculate trace of \f$A^2\f$
 *
 * This function uses Kahan sumation (see \ref kahan_sum) to reduce rounding
 * errors.
 *
 * The function needs \f$\mathcal{O}(N^2)\f$ operation for an \f$N\times N\f$
 * matrix.
 *
 * @param [in] A matrix
 * @retval trace \f$\mathrm{tr}\left(A^2\right)\f$
 */
double matrix_trace2(matrix_t *A)
{
    const size_t dim = A->dim;
    double *M = A->M;
    double array[dim];

    int N = dim;
    int incx = 1;
    int incy = N;

    for(size_t i = 0; i < dim; i++)
        array[i] = ddot_(&N, &M[dim*i], &incx, &M[i], &incy);

    return kahan_sum(array, dim);
}

/**
 * @brief Calculate Frobenius norm of \f$A\f$
 *
 * @param [in] A matrix
 * @retval |A| Frobenius norm of \f$A\f$
 */
double matrix_norm_frobenius(matrix_t *A)
{
    double norm = 0;
    double *M = A->M;

    for(size_t i = 0; i < A->dim2; i++)
        norm += pow_2(M[i]);

    return sqrt(norm);
}



/**
 * @brief Calculate \f$\log\det A\f$ for triangular matrix \f$A\f$
 *
 * This function calculates the logarithm of the determinant of the matrix
 * \f$A\f$ assuming \f$A\f$ is upper or lower triangular:
 * \f[
 *   \log\det A = \log \prod_j A_{jj} = \sum_j \log A_{jj}
 * \f]
 *
 * @param [in] A triangular matrix
 * @retval logdet \f$\log\det A\f$
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
 * @brief Calculate \f$\log\det(1+zA)\f$ for matrix \f$A\f$
 *
 * Compute \f$\log\det(1+zA)\f$ using LAPACK. The algorithm is chosen by detalg
 * and may be DETALG_QR, DETALG_LU or DETALG_CHOLESKY.
 *
 * If the Frobenius norm of \f$zA\f$ is smaller than 1, the function tries to
 * approximate \f$\log\det A\f$ using a Mercator series (if possible) to reduce
 * the complexity for an \f$N\times N\f$ matrix \f$A\f$ from
 * \f$\mathcal{O}(N^3)\f$ to \f$\mathcal{O}(N^2)\f$.
 *
 * @param [in,out] A matrix; will be overwritten.
 * @param [in]     z factor \f$z\f$
 * @param [in]     detalg algorithm to use (cholesky, lu or qr)
 * @retval logdet  \f$\log\det(1+zA)\f$
 */
double matrix_logdet_dense(matrix_t *A, double z, detalg_t detalg)
{
    /* ||zA|| = |z| ||A|| */
    const double norm = fabs(z)*matrix_norm_frobenius(A);

    if(norm < 1)
    {
        /* log(det(Id+zA)) ≈ z*tr(A) - z²/2 tr(A²) + ... */
        const double trA  = z*matrix_trace(A);
        const double trA2 = pow_2(z)*matrix_trace2(A);
        const double mercator = trA-trA2/2;
        const double error = fabs(pow_2(norm)/2+norm+log1p(-norm));
        const double rel_error = fabs(error/mercator);

        if(rel_error < 1e-8)
            return mercator;
    }

    /* A = Id+z*A */
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

    if(detalg == DETALG_CHOLESKY)
        /* cholesky decomposition */
        return matrix_logdet_cholesky(A, 'U');
    else if(detalg == DETALG_QR)
        /* QR decomposition */
        return matrix_logdet_qr(A);
    else
        /* LU decomposition */
        return matrix_logdet_lu(A);
}

/**
 * @brief Calculate \f$\log\det A\f$ using LU decomposition
 *
 * Calculate LU decomposition of \f$A\f$ and use \ref matrix_logdet_triangular to
 * calculate \f$\log\det A\f$.
 *
 * @param [in,out] A matrix
 * @retval logdet \f$\log\det A\f$
 */
double matrix_logdet_lu(matrix_t *A)
{
    int info = 0;
    int dim = (int)A->dim;
    if(dim <= 0)
        return NAN;
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

    TERMINATE(info != 0, "dgetrf returned %d", info);

    return matrix_logdet_triangular(A);
}


/**
 * @brief Calculate \f$\log \det A\f$ using Cholesky decomposition
 *
 * Calculate Cholesky decomposition of \f$A\f$ and use \ref matrix_logdet_triangular to
 * calculate \f$\log \det A\f$.
 *
 * Only the lower part of the matrix (uplo=L) or the upper part of the matrix
 * (uplo=U) are used.
 *
 * @param [in,out] A matrix
 * @param [in] uplo L or U
 * @retval logdet \f$\log \det A\f$
 */
double matrix_logdet_cholesky(matrix_t *A, char uplo)
{
    int info = 0;
    int dim = (int)A->dim;
    if(dim <= 0)
        return NAN;
    int lda = (int)A->lda;
    double *a = A->M;

    if(uplo != 'U' && uplo != 'u')
        uplo = 'L';

    dpotrf_(
        &uplo, /* UPLO U: Upper triangle is stored; L: Lower triangle is stored */
        &dim,  /* N number of columns of A */
        a,     /* matrix A to be factored */
        &lda,  /* LDA leading dimension of A */
        &info
    );

    TERMINATE(info != 0, "dpotrf returned %d", info);

    return 2*matrix_logdet_triangular(A);
}


/**
 * @brief Calculate \f$\log\det A\f$ using QR decomposition
 *
 * Calculate QR decomposition of \f$A\f$ and use \ref matrix_logdet_triangular
 * to calculate \f$\log\det A\f$.
 *
 * @param [in,out] A matrix
 * @retval logdet \f$\log\det A\f$
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

    TERMINATE(info != 0, "dgeqrf returned %d", info);

    return matrix_logdet_triangular(A);
}
