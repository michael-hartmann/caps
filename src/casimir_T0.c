#include <ctype.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "casimir_T0.h"

#include "integration_perf.h"
#include "gausslaguerre.h"
#include "libcasimir.h"
#include "sfunc.h"
#include "edouble.h"

#define PRECISION 1e-12
#define ORDER 50
#define LFAC 6.

void usage(FILE *stream) {
    fprintf(stderr,
"Usage: casimir_T0 [OPTIONS]\n"
"This program will calculate the free Casimir energy F(T=0,L/R) for the\n"
"plane-sphere geometry for given L/R and temperature T. The output is in scaled\n"
"units.\n"
"\n"
"Mandatory options:\n"
"    -x, --LbyR L/R\n"
"        Separation L between sphere and plane divided by radius of sphere,\n"
"        where L/R > 0.\n"
"\n"
"Further options:\n"
"    -l, --lscale\n"
"        Specify parameter lscale. The vector space has to be truncated for\n"
"        some value lmax. This program will use lmax=(R/L*lscale) (default: %g)\n"
"\n"
"    -L LMAX\n"
"        Set lmax to the value LMAX. When -L is specified, -l will be ignored\n"
"\n"
"    -c, --cores CORES\n"
"        Use CORES of processors for the calculation (default: 1)\n"
"\n"
"    -p, --precision\n"
"        Set precision to given value (default: %g)\n"
"\n"
"    -N, --order\n"
"        Order of Gauss-Laguerre integrateion (default: %d)\n"
"\n"
"    -h,--help\n"
"        Show this help\n",
    LFAC, PRECISION, ORDER);
}

double integrand(double xi, double LbyR, int lmax, double precision)
{
    casimir_t casimir;
    integration_perf_t int_perf;
    double v = 0, v0 = 0;
    int m = 0, use_trace = 0;
    const double Q = 1/(1+LbyR);

    casimir_init(&casimir, Q, xi);
    casimir_set_lmax(&casimir, lmax);

    casimir_integrate_perf_init(&int_perf, 1*casimir.T, casimir.lmax);

    while(m <= lmax)
    {
        double v_m;

        if(use_trace)
            v_m = -casimir_trM(&casimir, 1, m, &int_perf);
        else
        {
            v_m = casimir_logdetD(&casimir, 1, m, &int_perf);
            if(fabs(v_m) < 1e-8)
                use_trace = 1;
        }

        if(m == 0)
        {
            v0 = v_m;
            v_m /= 2;
        }

        v += v_m;

        if(fabs(v_m/v0) < precision)
            break;

        m++;
    }

    casimir_integrate_perf_free(&int_perf);
    casimir_free(&casimir);

    return v;
}


int main(int argc, char *argv[])
{
    int order = ORDER, cores = 1, lmax = 0, k;
    double F0, alpha, LbyR = -1, lfac = LFAC, precision = PRECISION;
    edouble integral = 0, *xk, *ln_wk;

    /* parse command line options */
    while (1)
    {
        int c;
        struct option long_options[] =
        {
            { "help",      no_argument,       0, 'h' },
            { "LbyR",      required_argument, 0, 'x' },
            { "lmax",      required_argument, 0, 'L' },
            { "order",     required_argument, 0, 'N' },
            { "lscale",    required_argument, 0, 'l' },
            { "cores",     required_argument, 0, 'c' },
            { "precision", required_argument, 0, 'p' },
            { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "x:L:N:l:c:p:h", long_options, &option_index);

        /* Detect the end of the options. */
        if(c == -1)
            break;

        switch (c)
        {
            case 0:
                /* If this option sets a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
            case 'x':
                LbyR = atof(optarg);
                break;
            case 'L':
                lmax = atoi(optarg);
                break;
            case 'c':
                cores = atoi(optarg);
                break;
            case 'l':
                lfac = atof(optarg);
                break;
            case 'p':
                precision = atof(optarg);
                break;
            case 'N':
                order = atoi(optarg);
                break;
            case 'h':
                usage(stdout);
                exit(0);

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort();
        }
    }

    if(LbyR < 0)
    {
        fprintf(stderr, "LbyR must be positive.\n\n");
        usage(stderr);
        return 1;
    }
    if(order < 1)
    {
        fprintf(stderr, "order must be positive.\n\n");
        usage(stderr);
        return 1;
    }
    if(precision <= 0)
    {
        fprintf(stderr, "precision must be positive.\n\n");
        usage(stderr);
        return 1;
    }
    if(cores <= 0)
    {
        fprintf(stderr, "cores must be positive.\n\n");
        usage(stderr);
        return 1;
    }
    if(lmax <= 0)
    {
        if(lfac > 0)
        {
            lmax = lfac/LbyR;
        }
        else if(lfac < 0)
        {
            fprintf(stderr, "lfac must be positive\n\n");
            usage(stderr);
            return 1;
        }
        else
        {
            fprintf(stderr, "lmax must be positive\n\n");
            usage(stderr);
            return 1;
        }
    }

    /* disable buffering */
    fflush(stdin);
    fflush(stderr);
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    /* estimate, cf. eq. (6.33) */
    alpha = 2*LbyR/(1+LbyR);
    order = gausslaguerre_nodes_weights(order, &xk, &ln_wk);

    printf("# LbyR  = %.15g\n", LbyR);
    printf("# alpha = %.15g\n", alpha);
    printf("# order = %d\n", order);
    printf("# prec  = %g\n", precision);
    printf("# lmax  = %d\n", lmax);
    printf("# cores = %d\n", cores);
    printf("#\n");

    /* do Gauss-Legendre quadrature */
    for(k = 0; k < order; k++)
    {
        const edouble x    = xk[k];
        const edouble ln_w = ln_wk[k];
        const double f = integrand(x/alpha, LbyR, lmax, precision);

        integral += expe(ln_w+x)*f;

        printf("# k=%d, x=%.15g, logdetD(xi = x/alpha)=%.15g\n", k, (double)x, f);
    }

    F0 = (double)(integral/alpha/M_PI);

    printf("#\n");
    printf("# L/R, order, alpha, F(T=0)\n");
    printf("%.15g, %d, %.15g, %.15g\n", LbyR, order, alpha, F0);

    return 0;
}
