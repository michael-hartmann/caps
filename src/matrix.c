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
#include "floattypes.h"

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

    A->M = xmalloc(dim2*sizeof(float64));
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

    fclose(stream);

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

void matrix_setall(matrix_t *A, float64 z)
{
    const size_t dim2 = A->dim2;
    float64 *a = A->M;

    for(size_t i = 0; i < dim2; i++)
        a[i] = z;
}

float64 matrix_logdet_lu(matrix_t *A)
{
    const size_t dim = A->dim;
    float64 sum, det = 0;
    float64 *a = A->M;

    for(size_t j = 0; j < dim; j++)
    {
        for(size_t i = 0; i < j+1; i++)
        {
            sum = 0;
            for(size_t k = 0; k < i; k++)
                sum += a[i*dim+k]*a[k*dim+j];
            a[i*dim+j] -= sum;
        }
        for(size_t i = j+1; i < dim; i++)
        {
            sum = 0;
            for(size_t k = 0; k < j; k++)
                sum += a[i*dim+k]*a[k*dim+j];
            a[i*dim+j] = (a[i*dim+j]-sum)/a[j*dim+j];
        }
    }

    for(size_t i = 0; i < dim; i++)
        det += log64(fabs64(a[i*dim+i]));

    return det;
}


/**
 * @brief calculate QR decomposition of matrix
 *
 * This function calculates the QR decomposition of a matrix M and \f$\log|\det(M)|\f$.
 *
 * @param [in,out] M matrix
 * @retval logdet \f$\log|\det M|\f$
 */
float64 matrix_logdet_qr(matrix_t *M)
{
    const size_t dim = M->dim;
    float64 *m = M->M;

    for(size_t j = 0; j < dim-1; j++)
        for(size_t i = j+1; i < dim; i++)
        {
            float64 c,s, Mij = m[i*dim+j];

            if(Mij != 0)
            {
                const float64 a = m[j*dim+j];
                const float64 b = Mij; /* b != 0 */

                if(a == 0)
                {
                    c = 0;
                    s = -copysign64(1, b);
                }
                else if(fabs64(b) > fabs64(a))
                {
                    const float64 t = a/b;
                    const float64 u = copysign64(sqrt64(1+t*t),b);
                    s = -1/u;
                    c = -s*t;
                }
                else
                {
                    const float64 t = b/a;
                    const float64 u = copysign64(sqrt64(1+t*t),a);
                    c = 1/u;
                    s = -c*t;
                }

                for(size_t n = j; n < dim; n++)
                {
                    const float64 Min = m[i*dim+n];
                    const float64 Mjn = m[j*dim+n];

                    m[i*dim+n] = c*Min + s*Mjn;
                    m[j*dim+n] = c*Mjn - s*Min;
                }

                /* m[i*dim+j] = 0; */
            }
        }

    float64 logdet = 0;
    for(size_t i = 0; i < dim; i++)
        logdet += log64(fabs64(m[i*dim+i]));

    return logdet;
}


/**
 * @brief Calculate \f$\log\det(\mathrm{Id}+z*M)\f$ for matrix M
 *
 * @param [in]     casimir casimir object
 * @param [in,out] M round trip matrix M (matrix elements given as logarithms); M will be overwritten.
 * @param [in]     M_sign signs of matrix elements M
 * @param [in]     z factor z in log(det(Id+z*M))
 * @retval logdet  \f$\log\det(\mathrm{Id}+z*M)\f$
 */
float64 matrix_logdet(matrix_t *A, float64 z, const char *detalg)
{
    float64 logdetD = 0;

    /* M = Id+z*M */
    /* XXX use BLAS routine XXX */
    {
        const size_t dim  = A->dim;
        const size_t dim2 = A->dim2;
        float64 *a = A->M;

        for(size_t i = 0; i < dim2; i++)
            a[i] *= z;

        for(size_t i = 0; i < dim; i++)
            a[i*dim+i] += 1;
    }

    logdetD = matrix_logdet_lu(A);
    //logdetD = matrix_logdet_lu_lapack(A);
    //logdetD = matrix_logdet_qr(A);

    return logdetD;
}

float64 matrix_logdet_lu_lapack(matrix_t *A)
{
    int info = 0;
    int dim = (int)A->dim;
    int rows = dim;
    int cols = dim;
    int ipiv[dim];
    double *a = A->M;

    dgetrf_(
        &rows, /* M number of rows of A */
        &cols, /* N number of columns of A */
        a,     /* matrix A to be factored */
        &dim,  /* LDA: leading dimension of A */
        ipiv,  /* pivot indices of dimension min(rows,cols) */
        &info
    );

    float64 logdet = 0;
    for(int i = 0; i < dim; i++)
        logdet += log64(fabs64(matrix_get(A, i, i)));

    return logdet;
}

float64 matrix_logdetIdmM_eig_lapack(matrix_t *A, double z)
{
    int dim = A->dim;
    char jobvl = 'N'; /* don't compute left eigenvectors */
    char jobvr = 'N'; /* don't compute right eigenvectors */
    double wr[dim]; /* real parts of eigenvalues */
    double wi[dim]; /* imaginary parts of eigenvalues */
    int lwork = 100*dim;
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

    double logdet = 0;
    for(int i = 0; i < dim; i++)
    {
        double lambda = wr[i];

        if(lambda < -z)
            logdet += log1p64(+z*lambda);
        else
            logdet += log64(fabs64(1-z*lambda));
    }

    return logdet;
}
