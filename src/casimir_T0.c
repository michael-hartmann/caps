#include <stdio.h>
#include <stdlib.h>

#include "integration_perf.h"
#include "gausslaguerre.h"
#include "libcasimir.h"
#include "sfunc.h"
#include "edouble.h"

#define PRECISION 1e-9
#define ORDER 50

double integrand(double LbyR, double xi);

double integrand(double LbyR, double xi)
{
    casimir_t casimir;
    integration_perf_t int_perf;
    double v = 0;
    int m = 0;
    double v0 = 0;
    const double Q = 1/(1+LbyR);
    const int lmax = MAX(20,6./LbyR);

    casimir_init(&casimir, Q, xi);
    casimir_set_lmax(&casimir, lmax);

    casimir_integrate_perf_init(&int_perf, 1*casimir.T, casimir.lmax);

    while(m <= lmax)
    {
        double v_m = casimir_logdetD(&casimir, 1, m, &int_perf);

        if(m == 0)
        {
            v0 = v_m;
            v_m /= 2;
        }

        v += v_m;

        if(fabs(v_m/v0) < PRECISION)
            break;

        m++;
    }

    casimir_integrate_perf_free(&int_perf);
    casimir_free(&casimir);

    return v;
}


int main(int argc, char *argv[])
{
    int k;
    double F0;
    const double LbyR = 0.1;
    const double alpha = 2*LbyR/(1+LbyR);
    edouble integral = 0;
    edouble *xk, *ln_wk;

    const int order = gausslaguerre_nodes_weights(ORDER, &xk, &ln_wk);

    printf("# LbyR  = %.15g\n", LbyR);
    printf("# alpha = %.15g\n", alpha);
    printf("# order = %d\n", order);
    printf("#\n");

    for(k = 0; k < order; k++)
    {
        const edouble x    = xk[k];
        const edouble ln_w = ln_wk[k];
        const double f = integrand(LbyR, x/alpha);

        integral += expe(ln_w+x)*f;

        printf("# k=%d, x=%.15g, f(x/alpha)=%.15g\n", k, (double)x, f);

        if(f == 0)
            break;
    }

    F0 = (double)(integral/alpha/M_PI);

    printf("#\n");
    printf("# L/R, order, alpha, F(T=0)\n");
    printf("%.15g, %d, %.15g, %.15g\n", LbyR, order, alpha, F0);

    return 0;
}
