#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "clapack.h"

#define index(m,n) ((n)*dim+(m))

#define matrix_set(M,m,n,v) M[(kl+ku+(m)-(n))+(n)*dim] = v;

void print_matrix(double *M, int dim)
{
    for(int m = 0; m < dim; m++)
    {
        double prod = 1;
        for(int n = 0; n < dim; n++)
        {
            prod *= M[n*dim+m];
            printf("%g ", M[n*dim+m]);
        }
        printf(" | %g\n", prod);
    }
}


int main(int argc, char *argv[])
{
    int info;
    int dim = 7;
    int ipiv[dim];
    double *M = malloc(dim*dim*sizeof(double));
    int kl = 2;
    int ku = 2;

    /* set matrix elements to 0 */
    for(int m = 0; m < dim*dim; m++)
        M[m] = 0;

    /* subdiagonals */
    for(int m = 0; m < dim; m++)
    {
        //M[index(kl+ku+m-m,m)] = 2;
        matrix_set(M,m,m,2);

        if((m-1) >= 0)
            //M[index(kl+ku+m-(m-1),m-1)] = 1;
            matrix_set(M,m,m-1,1);

        if((m+1) < dim)
            //M[index(kl+ku+m-(m+1),m+1)] = 1;
            matrix_set(M,m,m+1,1);
    }

    print_matrix(M,dim);


    int lda = dim;
    dgbtrf_(
        &dim, /* M number of rows of A */
        &dim, /* N number of columns of A */
        &kl,
        &ku,
        M,    /* matrix A to be factored */
        &lda, /* LDA: leading dimension of A */
        ipiv, /* pivot indices of dimension min(rows,cols) */
        &info
    );

    printf("info = %d\n", info);

    print_matrix(M,dim);

    return 0;
}
