#include <stdlib.h>
#include <stdio.h>

#include "casimir_cylinder.h"

#include "constants.h"
#include "bessel.h"
#include "matrix.h"

#include "cquadpack/include/cquadpack.h"


casimir_cp_t *casimir_cp_init(double R, double d)
{
    casimir_cp_t *self = malloc(sizeof(casimir_cp_t));

    self->R = R;
    self->d = d;
    self->H = R+d;

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

            double log_dIl1 = bessel_logIn(l1-1,R*q)+log1p(exp(bessel_logIn(l1+1,R*q)-bessel_logIn(l1-1,R*q)))-log(2);
            double log_dIl2 = bessel_logIn(l2-1,R*q)+log1p(exp(bessel_logIn(l2+1,R*q)-bessel_logIn(l2-1,R*q)))-log(2);

            double log_dKl1 = bessel_logKn(l1+1,R*q)+log1p(exp(bessel_logKn(l1-1,R*q)-bessel_logKn(l1+1,R*q)))-log(2);
            double log_dKl2 = bessel_logKn(l2+1,R*q)+log1p(exp(bessel_logKn(l2-1,R*q)-bessel_logKn(l2+1,R*q)))-log(2);

            /* matrix element ℳ_{l1,l2} */
            double elem = 0.5*(log_dIl1+log_dIl2-log_dKl1-log_dKl2) + bessel_logKn(l1+l2,2*H*q);

            matrix_set(D,i,j, -exp(elem));
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


static double __integrand_dirichlet(double x, void *args)
{
    casimir_cp_t *self = (casimir_cp_t *)args;
    double q = x/(2*self->d);
    return q*casimir_cp_dirichlet(self, q);
}

static double __integrand_neumann(double x, void *args)
{
    casimir_cp_t *self = (casimir_cp_t *)args;
    double q = x/(2*self->d);
    return q*casimir_cp_neumann(self, q);
}


int main(int argc, char *argv[])
{
    const double T = 0;
    const double R = 100e-6;
    const double d = 25e-6;
    const int lmax = 50;

    /* in units of hbar*c*L, i.e., E_PFA^DN / (hbar*c*L) */
    double E_PFA_DN = -M_PI*M_PI*M_PI/1920*sqrt(R/(2*d))/(d*d);
    double E_PFA = 2*E_PFA_DN;

    casimir_cp_t *c = casimir_cp_init(R,d);
    casimir_cp_set_lmax(c,lmax);

    double abserr;
    int neval, ier;

    /* in units of hbar*c*L */
    double E_D = dqagi(__integrand_dirichlet, 0, 1, 1e-8, 1e-8, &abserr, &neval, &ier, c)/(4*M_PI*2*d);

    /* in units of hbar*c*L */
    double E_N = dqagi(__integrand_neumann, 0, 1, 1e-8, 1e-8, &abserr, &neval, &ier, c)/(4*M_PI*2*d);

    /* in units of hbar*c*L */
    double E_EM = E_D+E_N;

    printf("# d/R, d, R, lmax, T, E_PFA/(L*hbar*c), E_D/E_PFA, E_N/E_PFA, E_EM/E_PFA\n");
    printf("%.15g, %.15g, %.15g, %d, %.15g, %.15g, %.15g, %.15g, %.15g\n", d/R, d, R, lmax, T, E_PFA, E_D/E_PFA, E_N/E_PFA, E_EM/E_PFA);

    casimir_cp_free(c);

    return 0;
}
