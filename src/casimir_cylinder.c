#include <stdlib.h>
#include <stdio.h>

#include "casimir_cylinder.h"

#include "argparse.h"
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

int casimir_cp_get_lmax(casimir_cp_t *self)
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

    double *cacheI = malloc((lmax+2)*sizeof(double));
    double *cacheK1 = malloc((lmax+2)*sizeof(double));
    double *cacheK2 = malloc(2*(lmax+1)*sizeof(double));
    for(int i = 0; i < lmax+2; i++)
    {
        cacheI[i]  = bessel_logIn(i, R*q);
        cacheK1[i] = bessel_logKn(i, R*q);
    }

    for(int i = 0; i < 2*(lmax+1); i++)
        cacheK2[i] = bessel_logKn(i,2*H*q);

    matrix_t *D = matrix_alloc(dim);

    for(int i = 0; i < dim; i++)
    {
        int l1 = i-lmax;

        double term_l1 = cacheI[abs(l1)]-cacheK1[abs(l1)];

        for(int j = i; j < dim; j++)
        {
            int l2 = j-lmax;

            /* matrix element ℳ_{l1,l2} */
            double term_l2 = cacheI[abs(l2)]-cacheK1[abs(l2)];

            double elem = exp(0.5*(term_l1+term_l2)+cacheK2[abs(l1+l2)]);
            matrix_set(D,i,j, -elem);
            matrix_set(D,j,i, -elem);
        }
    }

    free(cacheI);
    free(cacheK1);
    free(cacheK2);

    for(int i = 0; i < dim; i++)
        matrix_set(D,i,i, 1+matrix_get(D,i,i));

    const double logdet = matrix_logdet_cholesky(D, 'U');

    matrix_free(D);

    return logdet;
}

double casimir_cp_neumann(casimir_cp_t *self, double q)
{
    const int lmax = self->lmax;
    const int dim = 2*lmax+1;
    const double H = self->H, R = self->R;

    double *cacheI = malloc((lmax+2)*sizeof(double));
    double *cacheK1 = malloc((lmax+2)*sizeof(double));
    double *cacheK2 = malloc(2*(lmax+1)*sizeof(double));
    for(int i = 0; i < lmax+2; i++)
    {
        cacheI[i]  = bessel_logIn(i, R*q);
        cacheK1[i] = bessel_logKn(i, R*q);
    }

    for(int i = 0; i < 2*(lmax+1); i++)
        cacheK2[i] = bessel_logKn(i,2*H*q);

    matrix_t *D = matrix_alloc(dim);

    for(int i = 0; i < dim; i++)
    {
        int l1 = i-lmax;

        double log_dIl1 = cacheI [abs(l1-1)]+log1p(exp(cacheI [abs(l1+1)]-cacheI [abs(l1-1)]))-log(2);
        double log_dKl1 = cacheK1[abs(l1-1)]+log1p(exp(cacheK1[abs(l1+1)]-cacheK1[abs(l1-1)]))-log(2);

        for(int j = i; j < dim; j++)
        {
            int l2 = j-lmax;

            double log_dIl2 = cacheI [abs(l2-1)]+log1p(exp(cacheI [abs(l2+1)]-cacheI [abs(l2-1)]))-log(2);
            double log_dKl2 = cacheK1[abs(l2-1)]+log1p(exp(cacheK1[abs(l2+1)]-cacheK1[abs(l2-1)]))-log(2);

            /* matrix element ℳ_{l1,l2} */
            double elem = exp(0.5*(log_dIl1+log_dIl2-log_dKl1-log_dKl2) + cacheK2[abs(l1+l2)]);

            matrix_set(D,i,j, -elem);
            matrix_set(D,j,i, -elem);
        }
    }

    free(cacheI);
    free(cacheK1);
    free(cacheK2);

    for(int i = 0; i < dim; i++)
        matrix_set(D,i,i, 1+matrix_get(D,i,i));

    const double logdet = matrix_logdet_cholesky(D, 'U');

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


int main(int argc, const char *argv[])
{
    double epsrel = 1e-8, T = 0, R = NAN, d = NAN;
    int lmax = 0;

    const char *const usage[] = {
        "casimir_cylinder [options] [[--] args]",
        "casimir_cylinder [options]",
        NULL,
    };

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_GROUP("Mandatory options"),
        OPT_DOUBLE('R', "radius",      &R, "radius of cylinder in m", NULL, 0, 0),
        OPT_DOUBLE('d', "separation",  &d, "seperation between cylinder and plate in m", NULL, 0, 0),
        OPT_DOUBLE('T', "temperature", &T, "temperature in K", NULL, 0, 0),
        OPT_GROUP("Further options"),
        OPT_INTEGER('l', "lmax", &lmax, "dimension of vector space", NULL, 0, 0),
        OPT_DOUBLE('e', "epsrel", &epsrel, "relative error for integration", NULL, 0, 0),
        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usage, 0);
    argparse_describe(&argparse, "\nCompute the Casimir interaction in the cylinder-plane geometry.", NULL);
    argc = argparse_parse(&argparse, argc, argv);

    /* check arguments */
    if(isnan(R))
    {
        fprintf(stderr, "missing mandatory option: -R radius of sphere\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(R <= 0)
    {
        fprintf(stderr, "radius of sphere must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(isnan(d))
    {
        fprintf(stderr, "missing mandatory option: -d separation between cylinder and plate\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(d <= 0)
    {
        fprintf(stderr, "separation between cylinder and plate must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(T < 0)
    {
        fprintf(stderr, "temperature must be non-negative\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(epsrel <= 0)
    {
        fprintf(stderr, "epsrel must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }

    /* PFA for Dirichlet/Neumann in units of hbar*c*L, i.e., E_PFA^DN / (hbar*c*L) */
    double E_PFA_DN = -M_PI*M_PI*M_PI/1920*sqrt(R/(2*d))/(d*d);
    /* PFA for EM in units of hbar*c*L, i.e., E_PFA^DN / (hbar*c*L) */
    double E_PFA = 2*E_PFA_DN;

    casimir_cp_t *c = casimir_cp_init(R,d);
    
    if(lmax > 0)
        casimir_cp_set_lmax(c,lmax);

    double abserr;
    int neval, ier;

    /* energy Dirichlet in units of hbar*c*L */
    double E_D = dqagi(__integrand_dirichlet, 0, 1, 1e-8, epsrel, &abserr, &neval, &ier, c)/(4*M_PI*2*d);

    /* energy Neumann in units of hbar*c*L */
    double E_N = dqagi(__integrand_neumann, 0, 1, 1e-8, epsrel, &abserr, &neval, &ier, c)/(4*M_PI*2*d);

    /* energy EM in units of hbar*c*L */
    double E_EM = E_D+E_N;

    printf("# d/R, d, R, T, lmax, E_PFA/(L*hbar*c), E_D/E_PFA, E_N/E_PFA, E_EM/E_PFA\n");
    printf("%.15g, %.15g, %.15g, %.15g, %d, %.15g, %.15g, %.15g, %.15g\n", d/R, d, R, T, casimir_cp_get_lmax(c), E_PFA, E_D/E_PFA, E_N/E_PFA, E_EM/E_PFA);

    casimir_cp_free(c);

    return 0;
}
