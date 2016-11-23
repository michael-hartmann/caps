/**
 * @file   integration.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   November, 2016
 * @brief  Perform integration for arbitrary materials
 */

#include <assert.h>
#include <math.h>

#include "quadpack.h"

#include "utils.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"
#include "hash-table.h"


typedef struct
{
    int l1,l2,m,p;
    double tau,zmax,normalization;
} integrand_t;

static double f_bisect(double left, double right, int N, int nu, double tau, double log_c);
void I_estimate_width(int l1, int l2, double tau, double eps, double *a, double *b);
double I_zmax(int l1, int l2, double tau);
double I_integrand(double z, void *args_);
double I_log_integrand_max(int l1, int l2, int m, double tau);
double I(int l1, int l2, int m, int p, double tau, double epsrel);

static double f_bisect(double left, double right, int N, int nu, double tau, double log_c)
{
    for(int i = 0; i < N; i++)
    {
        double middle = (right+left)/2;

        double fl = nu*log(left)  -tau*left  -log_c;
        double fm = nu*log(middle)-tau*middle-log_c;
        //double fr = nu*log(right) -tau*right -log_c;

        if(fl*fm < 0)
            right = middle;
        else
            left = middle;
    }

    return (left+right)/2;
}

void I_estimate_width(int l1, int l2, double tau, double eps, double *a, double *b)
{
    double nu = l1+l2-2;
    double delta = 1e-3;

    /* f(z) = z^(l1+l2-2)*exp(-tau*z) */
    double zmax = (l1+l2-2)/tau;
    double log_c = log(eps)+nu*log(zmax)-tau*zmax;

    /* a */
    if(zmax <= 5)
        *a = 0;
    else
    {
        double left = delta;
        double right = zmax;
        int N = ceil((log(right-left)-log(delta))/M_LOG2);
        *a = f_bisect(left, right, N, nu, tau, log_c);
    }

    /* b */
    for(int i = 1;; i++)
    {
        double left = i*zmax;
        double right = (i+1)*zmax;

        double fl = nu*log(left) -tau*left -log_c;
        double fr = nu*log(right)-tau*right-log_c;

        if(fl*fr <= 0)
        {
            int N = ceil((log(right-left)-log(delta))/M_LOG2);
            *b = f_bisect(left, right, N, nu, tau, log_c);
            return;
        }
    }
}

double I_zmax(int l1, int l2, double tau)
{
    return (l1+l2-2)/tau;
}

double I_integrand(double z, void *args_)
{
    if(z == 0)
        return 0;

    integrand_t *args = (integrand_t *)args_;

    int l1 = args->l1;
    int l2 = args->l2;

    double exp_function  = exp(args->tau*(z-args->zmax)/(l1+l2));
    double factor = args->normalization * exp_function;

    double Pl1m, Pl2m;
    Pl1mPl2m(l1,l2, args->m,1+z,factor, &Pl1m, &Pl2m);

    return Pl1m*Pl2m/(z*(z+2)); /* rp is missing */
}


double I_log_integrand_max(int l1, int l2, int m, double tau)
{
    double zmax = I_zmax(l1,l2,tau);
    return lfac(2*l1)+lfac(2*l2)-(l1+l2)*M_LOG2-lfac(l1)-lfac(l1-m)-lfac(l2)-lfac(l2-m) + (l1+l2-2)*log(zmax)-tau*zmax;
}

double I(int l1, int l2, int m, int p, double tau, double epsrel)
{
    const double zmax = I_zmax(l1,l2,tau);
    const int lmax = MAX(l1,l2);

    const double normalization = exp( (Plm_estimate(lmax,m,1+zmax) - log(zmax))/lmax );

    integrand_t args = {
        .l1   = l1,
        .l2   = l2,
        .m    = m,
        .p    = p,
        .tau  = tau,
        .zmax = zmax,
        .normalization = normalization
    };

    double a,b;
    I_estimate_width(l1, l2, tau, 1e-5, &a, &b);

    double abserr, result1, result2, result3;
    int neval1, neval2, neval3, ier1, ier2, ier3;
    result1 = dqags(I_integrand, a, b, 0, epsrel, &abserr, &neval1, &ier1, &args); /* [0,a] */
    result2 = dqags(I_integrand, 0, a, 0, epsrel, &abserr, &neval2, &ier2, &args); /* [a,b] */
    result3 = dqagi(I_integrand, b, 1, 0, epsrel, &abserr, &neval3, &ier3, &args); /* [b,∞] */

    assert(ier1 == 0 && ier2 == 0 && ier3 == 0);

    return result1+result2+result3;
}

integration_t *casimir_integrate_init(casimir_t *self, int n, int m)
{
    return NULL;
}

void casimir_integrate_free(integration_t *integration)
{

}

int casimir_integrate(integration_t *int_obj, int l1, int l2, double v[8])
{
    return 0;
}

double casimir_integrate_A(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    return 0;
}

double casimir_integrate_B(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    return 0;
}

double casimir_integrate_C(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    return 0;
}

double casimir_integrate_D(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    return 0;
}

double casimir_integrate_I(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    /*
    int m = self->m;
    double tau = self->tau;
    */

    return 0;
}
