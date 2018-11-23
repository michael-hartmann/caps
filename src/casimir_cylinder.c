#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>

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
    self->dim = 2*self->lmax+1;

    self->verbose = 0;

    return self;
}

int casimir_cp_set_detalg(casimir_cp_t *self, detalg_t detalg)
{
    self->detalg = detalg;
    return 1;
}

detalg_t casimir_cp_get_detalg(casimir_cp_t *self)
{
    return self->detalg;
}

int casimir_cp_get_verbose(casimir_cp_t *self)
{
    return self->verbose;
}

int casimir_cp_set_verbose(casimir_cp_t *self, int verbose)
{
    return self->verbose = verbose;
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

static double __kernel(int i, int j, void *args_)
{
    kernel_args_t *args = (kernel_args_t *)args_;

    const int lmax = args->lmax;
    const int mu1 = i-lmax, mu2 = j-lmax;

    return exp(0.5*(args->cache_ratio[abs(mu1)]+args->cache_ratio[abs(mu2)]) + args->cacheK[abs(mu1+mu2)]);
}

kernel_args_t *kernel_init(casimir_cp_t *self, double q, char DN)
{
    const int lmax = self->lmax;
    const double R = self->R, H = self->H; /* H = R+d */

    kernel_args_t *args = malloc(sizeof(kernel_args_t));

    args->lmax = self->lmax;
    args->DN = DN;

    /* K_{µ1+µ2}(2*(R+d)*q) */
    args->cacheK = malloc(2*(lmax+1)*sizeof(double));
    for(int j = 0; j < 2*(lmax+1); j++)
        args->cacheK[j] = bessel_logKn(j,2*H*q);

    double *logIn = malloc((lmax+3)*sizeof(double));
    double *logKn = malloc((lmax+3)*sizeof(double));

    for(int j = 0; j < lmax+3; j++)
    {
        logIn[j] = bessel_logIn(j, R*q);
        logKn[j] = bessel_logKn(j, R*q);
    }

    args->cache_ratio = malloc((lmax+2)*sizeof(double));
    if(DN == 'D')
    {
        /* Dirichlet */
        for(int j = 0; j < lmax+2; j++)
            args->cache_ratio[j] = logIn[j]-logKn[j];
    }
    else
    {
        args->cache_ratio[0] = logIn[1]-logKn[1];

        /* Neumann */
        for(int j = 1; j < lmax+2; j++)
        {
            /* denom = -2K'_j(x); K'_j(x) = -1/2*[ K_{j+1}(x) + K_{j-1}(x) ] */
            double denom = logKn[j+1]+log1p(exp(logKn[j-1]-logKn[j+1]));

            /* num = 2I'_j(x); I'_j(x) = = 1/2*[ I_{j+1}(x) + I_{j-1}(x) ] = dI */
            double num = logIn[j-1]+log1p(exp(logIn[j+1]-logIn[j-1]));

            args->cache_ratio[j] = num-denom;
        }
    }

    free(logIn);
    free(logKn);

    return args;
}

void kernel_free(kernel_args_t *args)
{
    free(args->cache_ratio);
    free(args->cacheK);
    free(args);
}

double casimir_cp_logdetD(casimir_cp_t *self, double q, char DN)
{
    kernel_args_t *args = kernel_init(self, q, DN);

    double logdet = kernel_logdet(self->dim, __kernel, args, true, self->detalg);
    TERMINATE(isnan(logdet), "bc=%c, q=%.15g, logdet=nan", DN, q);

    kernel_free(args);

    if(self->verbose)
        printf("# %c, q=%.15g, logdet=%.15g\n", DN, q, logdet);

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
    return q*casimir_cp_logdetD(self, q, 'D');
}

static double __integrand_neumann(double x, void *args)
{
    casimir_cp_t *self = (casimir_cp_t *)args;
    double q = x/(2*self->d);
    return q*casimir_cp_logdetD(self, q, 'N');
}


int main(int argc, const char *argv[])
{
    double epsrel = 1e-8, eta = 6, T = 0, R = NAN, d = NAN;
    const char *detalg = NULL;
    int verbose = 0, lmax = 0;

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
        OPT_DOUBLE('n', "eta", &eta, "set eta", NULL, 0, 0),
        OPT_BOOLEAN('v', "verbose", &verbose, "be verbose", NULL, 0, 0),
        OPT_STRING('D', "detalg", &detalg, "algorithm to compute determinants (HODLR, LU, QR or CHOLESKY", NULL, 0, 0),
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
    if(eta <= 0)
    {
        fprintf(stderr, "eta must be positive\n\n");
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
    else
        casimir_cp_set_lmax(c,MAX(20,ceil(eta/d*R)));

    if(verbose)
        casimir_cp_set_verbose(c, verbose);

    printf("# R/d = %.15g\n", R/d);
    printf("# d = %.15g\n", d);
    printf("# R = %.15g\n", R);
    printf("# T = %.15g\n", T);
    printf("# lmax = %d\n", lmax);
    printf("# epsrel = %g\n", epsrel);

    if(detalg)
    {
        if(strcasecmp(detalg, "LU") == 0)
        {
            printf("# detalg = LU\n");
            casimir_cp_set_detalg(c, DETALG_LU);
        }
        else if(strcasecmp(detalg, "QR") == 0)
        {
            printf("# detalg = QR\n");
            casimir_cp_set_detalg(c, DETALG_QR);
        }
        else if(strcasecmp(detalg, "Cholesky") == 0)
        {
            printf("# detalg = Cholesky\n");
            casimir_cp_set_detalg(c, DETALG_QR);
        }
    }
    else
    {
        printf("# detalg = HODLR\n");
        casimir_cp_set_detalg(c, DETALG_HODLR);
    }

    printf("#\n");

    double integral, abserr, E_D, E_N;
    int neval, ier;

    /* energy Dirichlet in units of hbar*c*L */
    integral = dqagi(__integrand_dirichlet, 0, 1, 0, epsrel, &abserr, &neval, &ier, c);
    printf("# D: ier=%d, neval=%d, I=%.15g, absrel=%g\n", ier, neval, integral, fabs(abserr/integral));
    E_D = integral/(4*M_PI*2*d);

    /* energy Neumann in units of hbar*c*L */
    integral = dqagi(__integrand_neumann, 0, 1, 0, epsrel, &abserr, &neval, &ier, c);
    printf("# N: ier=%d, neval=%d, I=%.15g, absrel=%g\n", ier, neval, integral, fabs(abserr/integral));
    E_N = integral/(4*M_PI*2*d);

    /* energy EM in units of hbar*c*L */
    double E_EM = E_D+E_N;

    printf("#\n");

    printf("# d/R, d, R, T, lmax, E_PFA/(L*hbar*c), E_D/E_PFA, E_N/E_PFA, E_EM/E_PFA\n");
    printf("%.15g, %.15g, %.15g, %.15g, %d, %.15g, %.15g, %.15g, %.15g\n", d/R, d, R, T, casimir_cp_get_lmax(c), E_PFA, E_D/E_PFA, E_N/E_PFA, E_EM/E_PFA);

    casimir_cp_free(c);

    return 0;
}
