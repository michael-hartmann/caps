#include <stdio.h>
#include <math.h>

#include "argparse.h"
#include "cquadpack.h"

#include "libcaps.h"
#include "misc.h"

/* Compute the Casimir free energy at T=0 in the sphere-sphere geometry for
 * perfect reflectors. It is assumed that both spheres have the same radius
 * R=R1=R2.
 *
 * For more information see the user manual.
 */

typedef struct {
    double cutoff;
    caps_t *self;
} args_t;

static double logdetD2(caps_t *caps, int m, double xi_)
{
    size_t lmin, lmax;
    caps_estimate_lminmax(caps, m, &lmin, &lmax);
    int ldim = lmax-lmin;

    matrix_t *M = matrix_alloc(2*ldim);
    matrix_setall(M, NAN);

    /* compute the round-trip matrix M of the sphere-plane geometry */
    caps_M_t *self = caps_M_init(caps, m, xi_);
    for(int i = 0; i < ldim; i++)
    {
        const int l1 = i+lmin;
        for(int j = 0; j < ldim; j++)
        {
            const int l2 = j+lmin;
            
            matrix_set(M, i, j,           caps_M_elem(self, l1, l2, 'E', 'E'));
            matrix_set(M, i, j+ldim,      caps_M_elem(self, l1, l2, 'E', 'M'));
            matrix_set(M, i+ldim, j,      caps_M_elem(self, l1, l2, 'M', 'E'));
            matrix_set(M, i+ldim, j+ldim, caps_M_elem(self, l1, l2, 'M', 'M'));
        }
    }
    caps_M_free(self);

    /* M2 = M² = M*M: round-trip operator in the sphere-sphere geometry */
    matrix_t *M2 = matrix_mult(M, M, 1);
    matrix_free(M);

    /* check if we can use the trace approximation: logdet(1-A) =~ trace(A) */
    double trace2 = matrix_trace(M2);
    if(trace2 < 1e-8)
        return -trace2;

    /* logdet(Id-M²) */
    double logdet2 = matrix_logdet_dense(M2, -1, DETALG_LU);

    matrix_free(M2);

    return logdet2;
}

/* d = 2*L */
static double integrand(double xidbyc, void *args_)
{
    args_t *args = (args_t *)args_;
    caps_t *self = args->self;
    double xi_ = xidbyc*(self->L+self->R)/(2*self->L);
    double terms[4096] = { 0 };

    double cutoff = args->cutoff;

    for(unsigned int m = 0; m < sizeof(terms)/sizeof(terms[0]); m++)
    {
        terms[m] = logdetD2(self, m, xi_);
        //printf("m=%d, xi_=%g, %g\n", m, xi_, terms[m]);

        if(terms[m] == 0 || terms[m]/terms[0] < cutoff)
            break;
    }

    double logdet = 2*kahan_sum(terms, sizeof(terms)/sizeof(terms[0]))-terms[0];

    printf("# xi*d/c=%.12g, logdet=%.12g\n", xidbyc, logdet);

    return logdet;
}

int main(int argc, char *argv[])
{
    double R = NAN, d = NAN;
    double cutoff = 1e-10;
    double epsrel = 1e-6;
    double iepsrel = 1e-9;
    double eta = 5;
    int ldim = 0;

    const char *const usage[] = {
        "cass [options] [[--] args]",
        "cass [options]",
        NULL,
    };

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_GROUP("Mandatory options"),
        OPT_DOUBLE('R', "radius", &R, "radius of the two spheres in m", NULL, 0, 0),
        OPT_DOUBLE('d', "separation", &d, "smallest seperation between spheres in m", NULL, 0, 0),
        OPT_GROUP("Further options"),
        OPT_INTEGER('l', "ldim", &ldim, "dimension of vector space", NULL, 0, 0),
        OPT_DOUBLE('e', "epsrel", &epsrel, "relative error for integration", NULL, 0, 0),
        OPT_DOUBLE('i', "iepsrel", &iepsrel, "relative error for internal integration", NULL, 0, 0),
        OPT_DOUBLE('c', "cutoff", &cutoff, "cutoff for summation over m", NULL, 0, 0),
        OPT_DOUBLE('n', "eta", &eta, "set eta", NULL, 0, 0),
        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usage, 0);
    argparse_describe(&argparse, "\nCompute the Casimir interaction in the sphere-sphere geometry at T=0.", NULL);
    argc = argparse_parse(&argparse, argc, argv);

    /* check arguments */
    if(isnan(R))
    {
        fprintf(stderr, "missing mandatory option: -R radius of spheres\n\n");
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
        fprintf(stderr, "missing mandatory option: -d separation between spheres\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(d <= 0)
    {
        fprintf(stderr, "separation between cylinder and plate must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(ldim < 0)
    {
        fprintf(stderr, "ldim must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(epsrel <= 0)
    {
        fprintf(stderr, "epsrel must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(iepsrel <= 0)
    {
        fprintf(stderr, "iepsrel must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(cutoff <= 0)
    {
        fprintf(stderr, "cutoff must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }
    if(eta <= 0)
    {
        fprintf(stderr, "eta must be positive\n\n");
        argparse_usage(&argparse);
        return 1;
    }

    const double L = d*0.5;

    if(ldim == 0)
        ldim = fmax(20,R/d*eta);

    /* print information to stdout */
    caps_build(stdout, "# ");
    printf("#\n");
    printf("# sphere-sphere geometry\n");
    printf("# model: perfect reflectors\n");
    printf("# R1 = R2 = %.14g\n", R);
    printf("# d = %.14g\n", d);
    printf("# T = 0\n");
    printf("# ldim = %d\n", ldim);
    printf("# epsrel = %g\n", epsrel);
    printf("# iepsrel = %g\n", iepsrel);
    printf("# cutoff = %g\n", cutoff);
    printf("#\n");

    caps_t *self = caps_init(R, L);

    caps_set_epsrel(self, iepsrel);
    caps_set_ldim(self, ldim);
    caps_set_detalg(self, DETALG_LU);

    args_t args = {
        .cutoff = cutoff,
        .self   = self
    };

    double abserr;
    int neval, ier;
    double integral = dqagi(integrand, 0, +1, 0, epsrel, &abserr, &neval, &ier, &args);

    printf("#\n");
    printf("# ier = %d, neval = %d\n", ier, neval);
    printf("#\n");

    double E_ = integral/(2*CAPS_PI);

    printf("# R1, R2, L, E*d/(hbar*c), relative error (due to integration)\n");
    printf("%.12g, %.12g, %.12g, %.12g, %g\n", R, R, d, E_, fabs(abserr/integral));

    return 0;
}
