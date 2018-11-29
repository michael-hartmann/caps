#include <stdlib.h>
#include <math.h>

#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"

#include "hodlr.h"

class Kernel : public HODLR_Matrix
{
private:
    void *args;
    double (*get_kernel_element)(int, int, void *);

public:
    Kernel(unsigned N, double (kernel_)(int,int,void *), void *args_) : HODLR_Matrix(N)
    {
        args = args_;
        get_kernel_element = kernel_;
    };

    double getMatrixEntry(int m, int n)
    {
        double delta = (m==n) ? 1 : 0;
        return delta-get_kernel_element((int)m,(int)n,args);
    };
};


double hodlr_logdet(int dim, double (*callback)(int,int,void *), void *args, unsigned int nLeaf, double tolerance, int is_symmetric)
{
    // XXX
    int n_levels = 1+log2(((double)dim)/nLeaf);

    Kernel *K = new Kernel(dim, callback, args);

    HODLR_Tree *T = new HODLR_Tree(n_levels, tolerance, K);

    T->assembleTree(true);
    T->factorize();
    double logdet = T->logDeterminant();

    //delete T;
    delete K;

    return logdet;
}
