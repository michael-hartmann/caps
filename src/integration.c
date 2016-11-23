/**
 * @file   integration.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   November, 2016
 * @brief  Perform integration for arbitrary materials
 */

#include <math.h>

#include "quadpack.h"

#include "utils.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"
//#include "hash-table.h"

#define ZMAX(l1,l2,tau) (((l1)+(l2)-2.)/(tau))

typedef struct
{
    int l1,l2,m;
    polarization_t p;
    double tau,zmax,normalization;
    casimir_t *casimir;
} integrand_t;

static double f_bisect(double left, double right, int N, int nu, double tau, double log_c);
static void I_estimate_width(int l1, int l2, double tau, double eps, double *a, double *b);
static double I_integrand(double z, void *args_);

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

static void I_estimate_width(int l1, int l2, double tau, double eps, double *a, double *b)
{
    double nu = l1+l2-2;
    double delta = 1e-3;

    /* f(z) = z^(l1+l2-2)*exp(-tau*z) */
    double zmax = ZMAX(l1,l2,tau);
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

double I_integrand(double z, void *args_)
{
    if(z == 0)
        return 0;

    integrand_t *args = (integrand_t *)args_;

    const int l1 = args->l1, l2 = args->l2;

    const double tau = args->tau;
    const double k = tau/2*sqrt(pow_2(z)+2*z);

    double rTE, rTM;
    casimir_rp(args->casimir, tau/2, k, &rTE, &rTM);

    const double exp_function = exp(tau*(z-args->zmax)/(l1+l2));
    const double factor = args->normalization * exp_function;

    double Pl1m, Pl2m;
    Pl1mPl2m(l1, l2, args->m, 1+z, factor, &Pl1m, &Pl2m);

    if(args->p == TE)
        return rTE*Pl1m*Pl2m/(z*(z+2));
    else
        return rTM*Pl1m*Pl2m/(z*(z+2));
}


double casimir_integrate_I(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    const int lmax = MAX(l1,l2);
    const int m = self->m;
    const double tau = self->tau;
    const double zmax = ZMAX(l1,l2,tau);

    const double normalization = exp( (Plm_estimate(lmax,m,1+zmax) - log(zmax))/lmax );
    const double epsrel = self->epsrel;

    integrand_t args = {
        .l1   = l1,
        .l2   = l2,
        .m    = m,
        .p    = p,
        .tau  = tau,
        .zmax = zmax,
        .normalization = normalization,
        .casimir = self->casimir
    };

    double a,b;
    I_estimate_width(l1, l2, tau, 1e-5, &a, &b);

    double abserr, result1, result2, result3;
    int neval1, neval2, neval3, ier1, ier2, ier3;
    result1 = dqags(I_integrand, a, b, 0, epsrel, &abserr, &neval1, &ier1, &args); /* [0,a] */
    result2 = dqags(I_integrand, 0, a, 0, epsrel, &abserr, &neval2, &ier2, &args); /* [a,b] */
    result3 = dqagi(I_integrand, b, 1, 0, epsrel, &abserr, &neval3, &ier3, &args); /* [b,∞] */

    TERMINATE(ier1 != 0 || ier2 != 0 || ier3 != 0, "ier1=%d, ier2=%d, ier3=%d\nl1=%d, l2=%d, m=%d, tau=%g\na=%g, b=%g", ier1, ier2, ier3, l1,l2,m,tau,a,b);

    *prefactor = lfac(2*l1)+lfac(2*l2)-(l1+l2)*M_LOG2-lfac(l1)-lfac(l1-m)-lfac(l2)-lfac(l2-m) + (l1+l2-2)*log(zmax)-tau*zmax;

    return result1+result2+result3;
}

integration_t *casimir_integrate_init(casimir_t *casimir, int n, int m, double epsrel)
{
    integration_t *self = (integration_t *)xmalloc(sizeof(integration_t));

    self->casimir = casimir;
    self->n = n;
    self->m = m;
    self->tau = 2*n*casimir->T;

    self->epsrel = epsrel;

    return self;
}

void casimir_integrate_free(integration_t *integration)
{
    xfree(integration);
}


double casimir_integrate_A(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    const int m = self->m;
    const double A0 = -MPOW(l2+m)*pow_2(self->m);

    const double I = casimir_integrate_I(self, l1, l2, p, prefactor);
    *prefactor += casimir_lnLambda(l1, l2, m);

    return A0*I;
}

double casimir_integrate_B(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    double I, prefactor1, prefactor2, prefactor3, prefactor4;
    const int m = self->m;
    const double B0 = -MPOW(l2+m+1);

    const double I1 = casimir_integrate_I(self, l1-1, l2-1, p, &prefactor1);
    const double I2 = casimir_integrate_I(self, l1+1, l2-1, p, &prefactor2);
    const double I3 = casimir_integrate_I(self, l1-1, l2+1, p, &prefactor3);
    const double I4 = casimir_integrate_I(self, l1+1, l2+1, p, &prefactor4);

    const double denom = (2*l1+1)*(2*l2+1);
    I  = (l1+1)*(l1+m)*(l2+1)*(l2+m)/denom*I1*exp(prefactor1-prefactor4);
    I -=   l1*(l1-m+1)*(l2+1)*(l2+m)/denom*I2*exp(prefactor2-prefactor4);
    I -=   (l1+1)*(l1+m)*l2*(l2-m+1)/denom*I3*exp(prefactor3-prefactor4);
    I +=     l1*(l1-m+1)*l2*(l2-m+1)/denom*I4;

    *prefactor = prefactor4+casimir_lnLambda(l1,l2,m);

    return B0*I;
}

double casimir_integrate_C(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    double I, prefactor1, prefactor2;
    const int m = self->m;
    const double C0 = -MPOW(l2+m);

    const double I1 = casimir_integrate_I(self, l1, l2-1, p, &prefactor1);
    const double I2 = casimir_integrate_I(self, l1, l2+1, p, &prefactor2);

    const double denom = 2*l2+1;
    I  = -(l2+1)*(l2+m)/denom*I1*exp(prefactor1-prefactor2);
    I += l2*(l2-m+1)/denom*I2;
    *prefactor = prefactor2+casimir_lnLambda(l1,l2,m);

    return C0*I;
}

double casimir_integrate_D(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    return MPOW(l1+l2+1)*casimir_integrate_C(self, l2, l1, p, prefactor);
}

int casimir_integrate(integration_t *int_obj, int l1, int l2, double v[8])
{
    return 0;
}
