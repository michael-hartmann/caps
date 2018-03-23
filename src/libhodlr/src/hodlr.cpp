#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"

#include "hodlr.h"

class Kernel : public HODLR_Matrix
{
private:
    void *args;
    double (*get_kernel_element)(int, int, void *);

public:
    Kernel(unsigned N, double (kernel_)(int,int,void *), void *args_)
    {
        args = args_;
        get_kernel_element = kernel_;
    };

    double get_Matrix_Entry(const unsigned m, const unsigned n)
    {
        return -get_kernel_element((int)m,(int)n,args);
    };
};


double hodlr_logdet_diagonal(int dim, double (*callback)(int,int,void *), void *args, double *diagonal, unsigned int nLeaf, double tolerance, int is_symmetric)
{
    double logdet = NAN;
    char s = 'n'; /* non symmetric by default */

    if(is_symmetric)
        s = 's';

    Kernel kernel(dim, callback, args);

    HODLR_Tree<Kernel>* A = new HODLR_Tree<Kernel>(&kernel, dim, nLeaf);

    VectorXd d = VectorXd::Ones(dim);
    for(int n = 0; n < dim; n++)
        d(n) = 1-diagonal[n];

    A->assemble_Matrix(d, tolerance, s);
    A->compute_Factor();
    A->compute_Determinant(logdet);

    delete A;

    return logdet;
}

double hodlr_logdet(int dim, double (*callback)(int,int,void *), void *args, unsigned int nLeaf, double tolerance, int is_symmetric)
{
    double diagonal[dim];

    for(int n = 0; n < dim; n++)
        diagonal[n] = callback(n,n,args);

    return hodlr_logdet_diagonal(dim, callback, args, diagonal, nLeaf, tolerance, is_symmetric);
}
