#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "bessel.h"
#include "plm.h"
#include "logfac.h"
#include "matrix.h"
#include "quadpack.h"

#include "casimir_scalar.h"

static double log_integrand(double x, int l1, int l2, int m, double tau)
{
    const double log_sinhx = log(sinh(x));
    const double coshx = cosh(x);

    if(l1 == l2)
        return log_sinhx -tau*coshx + 2*Plm(l1,m,coshx);
    else
        return log_sinhx -tau*coshx + Plm(l1,m,coshx) + Plm(l2,m,coshx);
}

static double integrand(double x, void *args_)
{
    integrand_t *args = (integrand_t *)args_;
    return exp(log_integrand(x,args->l1,args->l2,args->m,args->tau)-args->max);
}

static double integrate(int l1, int l2, int m, double tau, double epsrel)
{
    integrand_t args;

    double xmax  = asinh((l1+l2+1)/tau);
    double width = 1/sqrt(tau*cosh(xmax));

    args.l1  = l1;
    args.l2  = l2;
    args.m   = m;
    args.tau = tau;
    args.max = log_integrand(xmax, l1, l2, m, tau);

    double a = fmax(0, xmax-8*width);
    double b = xmax+8*width;

    /* perform integrations in interval [a,b] */
    int neval = 0, ier = 0;
    double abserr = 0;

    /* I2: [a,b] */
    double I = dqags(integrand, a, b, 0, epsrel, &abserr, &neval, &ier, &args);

    #if 0
    printf("l1=%d, l2=%d, m=%d, tau=%g\n", l1, l2, m, tau);
    printf("width=%g\n", width);
    printf("xmax=%g, max=%g\n", xmax, args.max);
    printf("a=%g, b=%g\n", a, b);
    printf("I=%g\n", I);
    #endif

    return args.max+log(I);
}

static double casimir_rs(int l, double chi, char Y)
{
    if(Y == 'N')
        abort();
    else
    {
        double Inu, Knu;
        bessel_lnInuKnu(l, chi, &Inu, &Knu);

        return log(M_PI/2) + Inu-Knu;
    }
}

static double kernel(int i, int j, void *args_)
{
    double rp, rs, I;
    args_t *args = (args_t *)args_;

    const double LbyR = args->LbyR;
    const int m = args->m;
    const double xi_ = args->xi_;
    const int l1 = i+m, l2 = j+m;

    /* N_l1^m * N_l2^m */
    const double C = (logi(2*l1+1)+logi(2*l2+1)-log(4)+lfac(l1-m)+lfac(l2-m)-lfac(l1+m)-lfac(l2+m))/2;

    if(args->X == 'N')
        rp = 1;
    else
        rp = -1;

    if(isnan(args->cache_rs[i]))
        args->cache_rs[i] = casimir_rs(l1, xi_/(1+LbyR), args->Y);
    if(isnan(args->cache_rs[j]))
        args->cache_rs[j] = casimir_rs(l2, xi_/(1+LbyR), args->Y);

    rs = (args->cache_rs[i]+args->cache_rs[j])/2;

    I = integrate(l1,l2,m,2*xi_,args->epsrel);

    double v = C + rs + I;
    return -2*rp*exp(v);
}

double logdetD(double LbyR, double xi_, int m, int ldim, char X, char Y, double epsrel)
{
    args_t args;

    args.LbyR   = LbyR;
    args.xi_    = xi_;
    args.m      = m;
    args.X      = X;
    args.Y      = Y;
    args.epsrel = epsrel;

    args.cache_rs = xmalloc(ldim*sizeof(double));
    for(int i = 0; i < ldim; i++)
        args.cache_rs[i] = NAN;

    double v = kernel_logdet(ldim, kernel, &args, 1, DETALG_HODLR);

    xfree(args.cache_rs);

    return v;
}

int main(int argc, char *argv[])
{
    const double LbyR = 0.1;
    const int ldim = 6/LbyR;
    const char X = 'D', Y = 'D';
    double epsrel = 1e-6;
    double xi_ = 1;
    int m = 1;

    fflush(stdin);
    fflush(stderr);
    setvbuf(stdout, NULL, _IONBF, 0); 
    setvbuf(stderr, NULL, _IONBF, 0);

    double v = logdetD(LbyR, xi_, m, ldim, X, Y, epsrel);
    printf("%g\n", v);

    return 0;
}
