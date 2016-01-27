#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <qd/qd_real.h>
#include "floattypes.h"
#include "matrix.h"
#include "utils.h"

using std::cout;
using std::endl;

typedef struct {
    int size;
    dd_real *M;
} matrix_floatdd;

matrix_floatdd *matrix_floatdd_alloc(int size)
{
    matrix_floatdd *A = (matrix_floatdd *)malloc(sizeof(matrix_floatdd));
    if(A == NULL)
        return NULL;

    A->size = size;
    A->M    = (dd_real *)malloc(size*size*sizeof(dd_real));
    if(A->M == NULL)
    {
        free(A);
        return NULL;
    }

    return A;
}

void matrix_floatdd_free(matrix_floatdd *A)
{
    if(A != NULL)
    {
        if(A->M != NULL)
            free(A->M);

        free(A);
    }
}

double matrix_floatdd_logdet_qr(matrix_floatdd *M)
{
    const int dim = M->size;
    dd_real *m = M->M;

    for(int j = 0; j < dim-1; j++)
        for(int i = j+1; i < dim; i++)
        {
            dd_real c,s, Mij = m[i*dim+j];

            if(!Mij.is_zero())
            {
                const dd_real a = m[j*dim+j];
                const dd_real b = Mij;

                if(b.is_zero())
                {
                    //c = copysign128(1,a);
                    c = a.is_positive() ? dd_real(1) : dd_real(-1);
                    s = dd_real(0);
                }
                else if(a.is_zero())
                {
                    c = dd_real(0);
                    //s = -copysign128(1, b);
                    s = b.is_positive() ? dd_real(-1.) : dd_real(1.);
                }
                else if(fabs(b) > fabs(a))
                {
                    const dd_real t = a/b;
                    const dd_real x = sqrt(dd_real(1)+t*t);
                    //const dd_real u = copysign128(sqrt(1+t*t),b);
                    const dd_real u = b.is_positive() ? x : -x;
                    s = -dd_real(1)/u;
                    c = -s*t;
                }
                else
                {
                    const dd_real t = b/a;
                    const dd_real x = sqrt(dd_real(1)+t*t);
                    //const dd_real u = copysign128(sqrt128(1+t*t),a);
                    const dd_real u = a.is_positive() ? x : -x;

                    c = dd_real(1)/u;
                    s = -c*t;
                }

                for(int n = 0; n < dim; n++)
                {
                    const dd_real Min = m[i*dim+n];
                    const dd_real Mjn = m[j*dim+n];

                    m[i*dim+n] = c*Min + s*Mjn;
                    m[j*dim+n] = c*Mjn - s*Min;
                }

                /* m[i*dim+j] = 0; */
            }
        }

    dd_real det = 0;
    for(int i = 0; i < dim; i++)
        det += log(fabs(m[i*dim+i]));

    return to_double(det);
}

extern "C" double matrix_floatdd_logdet(matrix_float80 *M, matrix_sign_t *M_sign, const char *type)
{
    unsigned int oldcw;
    fpu_fix_start(&oldcw);

    double logdet;
    const int dim = M->size;

    matrix_floatdd *Add = matrix_floatdd_alloc(dim);

    for(int i = 0; i < dim*dim; i++)
    {
        const float80 elem = M->M[i];
        const double hi = elem;
        const double lo = elem-hi;
        Add->M[i] = M_sign->M[i]*exp(dd_real(hi,lo));
    }

    logdet = matrix_floatdd_logdet_qr(Add);
    matrix_floatdd_free(Add);

    fpu_fix_end(&oldcw);

    return logdet;
}
