#include <stdlib.h>
#include <stdio.h>

#include "casimir_cylinder.h"

#include "constants.h"
#include "bessel.h"
#include "matrix.h"


casimir_cp_t *casimir_cp_init(double R, double d)
{
    casimir_cp_t *self = malloc(sizeof(casimir_cp_t));

    self->R = R;
    self->d = d;
    self->H = d+R;

    self->lmax = MAX(5*ceil(R/d),10);

    return self;
}

int casimir_cp_get_lmax(casimir_cp_t *self, int lmax)
{
    return self->lmax;
}

int casimir_cp_set_lmax(casimir_cp_t *self, int lmax)
{
    if(lmax < 1)
        return 0;

    self->lmax = lmax;
    return 1;
}

double casimir_cp_dirichlet(casimir_cp_t *self, double q)
{
    const int lmax = self->lmax;
    const int dim = 2*lmax+1;
    const double H = self->H, R = self->R;

    matrix_t *D = matrix_alloc(dim);

    for(int i = 0; i < dim; i++)
    {
        int l1 = i-lmax;

        for(int j = i; j < dim; j++)
        {
            int l2 = j-lmax;

            /* matrix element ℳ_{l1,l2} */
            double elem = exp(0.5*(bessel_logIn(l1,R*q)+bessel_logIn(l2,R*q)-bessel_logKn(l1,R*q)-bessel_logKn(l2,R*q))+bessel_logKn(l1+l2,2*H*q));
            matrix_set(D,i,j, -elem);
            matrix_set(D,j,i, -elem);
        }
    }

    for(int i = 0; i < dim; i++)
        matrix_set(D,i,i, 1+matrix_get(D,i,i));

    const double logdet = matrix_logdet_lu(D);

    matrix_free(D);

    return logdet;
}

double casimir_cp_neumann(casimir_cp_t *self, double q)
{
    const int lmax = self->lmax;
    const int dim = 2*lmax+1;
    const double H = self->H, R = self->R;

    matrix_t *D = matrix_alloc(dim);

    for(int i = 0; i < dim; i++)
    {
        int l1 = i-lmax;

        for(int j = 0; j < dim; j++)
        {
            int l2 = j-lmax;

            /* matrix element ℳ_{l1,l2} */
            /* XXX symmetrize XXX */
            //double elem = exp( (iv(l1-1,z)+iv(l1+1,z)) / (kv(l2-1,z)+kv(l2+1,z))* kv(l1+l2,2*H*q) )
            double elem = (bessel_In(l1-1,R*q)+bessel_In(l1+1,R*q))/(bessel_Kn(l2-1,R*q)+bessel_Kn(l2+1,R*q)) * bessel_Kn(l1+l2,2*H*q);
            matrix_set(D,i,j, -elem);
        }
    }

    for(int i = 0; i < dim; i++)
        matrix_set(D,i,i, 1+matrix_get(D,i,i));

    const double logdet = matrix_logdet_lu(D);

    matrix_free(D);

    return logdet;
}

void casimir_cp_free(casimir_cp_t *self)
{
    free(self);
}

int main(int argc, char *argv[])
{
    const double q = 1;
    const double R = 100e-6;
    const double d = 200e-6;
    //double H = R+a

    printf("R = %g\n", R);
    printf("d = %g\n", d);
    printf("q = %g\n", q);

    printf("\n");

    casimir_cp_t *c = casimir_cp_init(R,d);
    double dirichlet = casimir_cp_dirichlet(c, q);
    double neumann   = casimir_cp_neumann(c, q);
    casimir_cp_free(c);

    printf("dirichlet=%.15g\n", dirichlet);
    printf("neumann  =%.15g\n", neumann);

    return 0;
}
