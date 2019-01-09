#include <stdlib.h>
#include <math.h>

#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"

#include "hodlr.h"

class Kernel : public HODLR_Matrix
{
private:
    void *args;
    double *diag;
    double (*get_kernel_element)(int, int, void *);

public:
    Kernel(unsigned N, double (kernel)(int,int,void *), double *diag, void *args) : HODLR_Matrix(N)
    {
        this->args = args;
        this->diag = diag;
        this->get_kernel_element = kernel;
    };

    double getMatrixEntry(int m, int n)
    {
        if(m==n)
            return 1-diag[m];

        return -get_kernel_element(m,n,args);
    };
};


double hodlr_logdet_diagonal(int dim, double (*callback)(int,int,void *), void *args, double *diagonal, unsigned int nLeaf, double tolerance, int sym_spd)
{
    /* nLeaf is the size (number of rows of the matrix) of the smallest block
     * at the leaf level. The number of levels in the tree is given by
     * n_levels=log_2(dim/nLeaf).
     */
    const int n_levels = fmax(1,log2(((double)dim)/nLeaf));

    Kernel K = Kernel(dim, callback, diagonal, args);

    HODLR_Tree *T = new HODLR_Tree(n_levels, tolerance, &K);

    int sym = sym_spd > 0 ? 1 : 0; /* matrix is symmetric */
    int spd = sym_spd > 1 ? 1 : 0; /* matrix is symmetric positive definite */
    T->assembleTree(sym, spd);
    T->factorize();
    double logdet = T->logDeterminant();

    delete T;

    return logdet;
}

/** Wrapper for hodlr_logdet_diagonal
 *
 * See \ref hodlr_logdet_diagonal.
 */
double hodlr_logdet(int dim, double (*callback)(int,int,void *), void *args, unsigned int nLeaf, double tolerance, int sym_spd)
{
    /* allocate memory */
    double *diag = (double *)malloc(((size_t)dim)*sizeof(double));
    if(diag == NULL)
        return NAN;

    /* save diagonal elements into diag */
    for(int m = 0; m < dim; m++)
        diag[m] = callback(m,m,args);

    double logdet = hodlr_logdet_diagonal(dim, callback, args, diag, nLeaf, tolerance, sym_spd);

    /* free memory */
    free(diag);

    return logdet;
}
