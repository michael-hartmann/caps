/**
 * @file   integration.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   November, 2016
 * @brief  Perform integration for arbitrary materials
 */

#include <math.h>
#include <stdint.h>
#include <inttypes.h>

#include "quadpack.h"

#include "utils.h"
#include "gaunt.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"
#include "hash-table.h"

/* arguments for integrand in function I_integrand */
typedef struct
{
    int nu,m;
    polarization_t p; /* TE or TM */
    double tau,zmax,normalization;
    casimir_t *casimir;
} integrand_t;

/* entry of cache constisting of the value of the integral; exp(prefactor)*I */
typedef struct
{
    double v;
    sign_t sign;
} cache_entry_t;

/* allocate memory and create new cache entry */
static cache_entry_t *cache_entry_create(double v, sign_t sign)
{
    cache_entry_t *entry = xmalloc(sizeof(cache_entry_t));
    entry->v = v;
    entry->sign = sign;
    return entry;
}

/* free memory for cache entry */
static void cache_entry_destroy(void *entry)
{
    xfree(entry);
}

/* Create a hash from l1, l2 and p. The hash function breaks if l1 or l2 >
 * 2^31. This, however, exceeds the ressources available by orders of
 * magnitude. Other things will break earlier...
 */
static uint64_t hash(int l1, int l2, polarization_t p)
{
    uint64_t l1_ = l1;
    uint64_t l2_ = l2;
    uint64_t p_  = (p == TM) ? 1 : 0;

    return (l1_ << 32) | (l2_ << 1) | p_;
}

static double ZMAX(int nu, int m, double tau)
{
    return nu/tau;
}

static double f_bisect(double left, double right, int N, int nu, double tau, double log_c)
{
    /* We want to solve
     *      x^ν*exp(-τx) = c, c>0.
     * Taking the logarithm of both sides yields
     *      ν*log(x)-τx-log(c) = 0.
     * We find the solution using bisection method (N times).
     */
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

static void K_estimate_width(int nu, int m, double tau, double eps, double *a, double *b)
{
    double delta = 1e-2;

    /* f(z) = z^ν*exp(-τz) */
    double zmax = ZMAX(nu,m,tau);
    double log_c = log(eps)+nu*log(zmax)-tau*zmax;

    /* a */
    {
        double left = 0, right = zmax;
        int N = ceil((log(right-left)-log(delta))/M_LOG2);
        *a = f_bisect(left, right, N, nu, tau, log_c);
    }

    /* b */
    {
        if(nu == 1)
        {
            /* k(z) = rp*exp(-τz) */
            *b = -log(eps)/tau;
            return;
        }

        for(int i = 1;; i++)
        {
            double left = i*zmax, right = (i+1)*zmax;

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
}

static double K_integrand(double z, void *args_)
{
    if(z <= 0)
        return 0;

    integrand_t *args = (integrand_t *)args_;

    const int nu = args->nu, m = args->m;

    const double tau = args->tau;
    const double z2p2z = z*(z+2); /* z²+2z */
    const double k = tau/2*sqrt(z2p2z);

    double rTE, rTM;
    casimir_rp(args->casimir, tau/2, k, &rTE, &rTM);

    const double exp_function = exp(tau*(z-args->zmax)/nu);
    const double factor = args->normalization * exp_function;

    const double v = Plm(nu,2*m,1+z,factor)/z2p2z;

    //fprintf(stderr, "%g %g\n", z, v);

    if(args->p == TE)
        return rTE*v;
    else
        return rTM*v;
}

static double _casimir_integrate_K(integration_t *self, int nu, polarization_t p, sign_t *sign)
{
    const int m = self->m;

    const double tau = self->tau;
    const double zmax = ZMAX(nu,m,tau);
    const double epsrel = self->epsrel;

    double log_normalization;

    if(zmax > 0)
        log_normalization = (Plm_estimate(nu,2*m,1+zmax)-log(zmax*(zmax+2)))/nu;
    else
        /* nu == 1 => I(z) = rp*exp(-τz) */
        log_normalization = 0;

    integrand_t args = {
        .nu   = nu,
        .m    = m,
        .p    = p,
        .tau  = tau,
        .zmax = zmax,
        .normalization = exp(log_normalization),
        .casimir = self->casimir
    };

    double a,b;
    K_estimate_width(nu, m, tau, 1e-6, &a, &b);

    double abserr1 = 0, abserr2, abserr3, result1, result2, result3;
    int neval1, neval2, neval3, ier1, ier2, ier3;
    result2 = dqags(K_integrand, a, b, 0, epsrel, &abserr1, &neval2, &ier2, &args); /* [a,b] */
    abserr2 = abserr3 = abserr1;
    result1 = dqags(K_integrand, 0, a, 0, 1, &abserr2, &neval1, &ier1, &args); /* [0,a] */
    result3 = dqagi(K_integrand, b, 1, 0, 1, &abserr3, &neval3, &ier3, &args); /* [b,∞] */

    //printf("ier1=%d, ier2=%d, ier3=%d, nu=%d, m=%d, tau=%g, a=%g, b=%g\n", ier1, ier2, ier3, nu,m,tau,a,b);
    TERMINATE(ier1 != 0 || ier2 != 0 || ier3 != 0, "ier1=%d, ier2=%d, ier3=%d, nu=%d, m=%d, tau=%g, a=%g, b=%g", ier1, ier2, ier3, nu,m,tau,a,b);

    double sum = result1+result2+result3;
    *sign = copysign(1, sum);
    return log(fabs(sum))-tau*zmax + log_normalization*nu;
}


double casimir_integrate_K(integration_t *self, int nu, polarization_t p, sign_t *sign)
{
    HashTable *hash_table = self->hash_table_K;
    const uint64_t key = hash(0,nu,p);
    cache_entry_t *entry = hash_table_lookup(hash_table, key);

    if(entry)
    {
        /* lookup successful */
        *sign = entry->sign;
        return entry->v;
    }
    else
    {
        /* compute and save integral */
        double K = _casimir_integrate_K(self, nu, p, sign);
        entry = cache_entry_create(K, *sign);
        hash_table_insert(hash_table, key, entry);
        return K;
    }
}

static double _casimir_integrate_I(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign)
{
    const int m = self->m;
    const int qmax = gaunt_qmax(l1,l2,m);
    const int nu = l1+l2;
    double a_tilde[qmax+1];
    double v[qmax+1];
    sign_t s[qmax+1];

    gaunt(l1, l2, m, a_tilde);

    int q;
    for(q = 0; q <= qmax; q++)
    {
        double aq = a_tilde[q];
        double K = casimir_integrate_K(self, nu-2*q, p, &s[q]);
        //printf("q=%d, nu=%d, K=%.10e, a_tilde=%g, tau=%g\n", q, nu-2*q, K*exp(pre), a_tilde[q], self->tau);
        
        s[q] *= copysign(1, aq);
        v[q] = K+log(fabs(aq));

        if(v[q]-v[0] < -70)
            break;
    }

    return gaunt_log_a0(l1,l2,m)+logadd_ms(v, s, q, sign);
}

double casimir_integrate_I(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign)
{
    const int m = self->m;

    if(l1 < m || l2 < m)
    {
        *sign = 0;
        return -INFINITY;
    }

    /* simplification for perfect conductors */
    if(self->is_pc && p == TE)
    {
        const double v = casimir_integrate_I(self, l1, l2, TM, sign);
        *sign = -1;
        return v;
    }

    /* I(l1,l2) = I(l2,l1)
     * Make sure that l1 >= l2
     */
    if(l1 < l2)
        swap(&l1, &l2);

    HashTable *hash_table = self->hash_table_I;
    const uint64_t key = hash(l1,l2,p);
    cache_entry_t *entry = hash_table_lookup(hash_table, key);

    if(entry)
    {
        /* lookup successful */
        *sign = entry->sign;
        return entry->v;
    }
    else
    {
        /* compute and save integral */
        double I = _casimir_integrate_I(self, l1, l2, p, sign);
        entry = cache_entry_create(I, *sign);
        hash_table_insert(hash_table, key, entry);
        return I;
    }
}

integration_t *casimir_integrate_init(casimir_t *casimir, int n, int m, double epsrel)
{
    integration_t *self = (integration_t *)xmalloc(sizeof(integration_t));

    self->casimir = casimir;
    self->n = n;
    self->m = m;
    self->tau = 2*n*casimir->T;
    self->epsrel = epsrel;

    self->hash_table_I = hash_table_new(cache_entry_destroy);
    self->hash_table_K = hash_table_new(cache_entry_destroy);

    if(isinf(casimir_epsilonm1(casimir, INFINITY)))
        self->is_pc = true;
    else
        self->is_pc = false;

    return self;
}

void casimir_integrate_free(integration_t *integration)
{
    hash_table_free(integration->hash_table_I);
    hash_table_free(integration->hash_table_K);
    xfree(integration);
}


double casimir_integrate_A(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign)
{
    const int m = self->m;
    if(m == 0)
    {
        *sign = 0;
        return -INFINITY;
    }

    const double I1 = casimir_integrate_I(self, l1, l2, p, sign);
    const double log_A0 = 2*logi(m)+casimir_lnLambda(l1,l2,m)-self->tau;

    *sign *= -MPOW(l2);

    return log_A0+I1;
}

double casimir_integrate_B(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign)
{
    const int m = self->m;

    if(m == 0)
    {
        *sign = 0;
        return 0;
        /*
        I = casimir_integrate_I(self, l1, l2, p, &prefactor1);
        *prefactor = prefactor1+casimir_lnLambda(l1,l2,m)-self->tau;
        return B0*I;
        */
    }

    const double log_B0 = casimir_lnLambda(l1,l2,m)-self->tau;
    sign_t sign1, sign2, sign3, sign4;
    const double log_I1 = casimir_integrate_I(self, l1-1, l2-1, p, &sign1);
    const double log_I2 = casimir_integrate_I(self, l1+1, l2-1, p, &sign2);
    const double log_I3 = casimir_integrate_I(self, l1-1, l2+1, p, &sign3);
    const double log_I4 = casimir_integrate_I(self, l1+1, l2+1, p, &sign4);

    double I;
    const double denom = (2*l1+1.)*(2*l2+1.);
    I  = (l1+1.)*(l1+m)*(l2+1.)*(l2+m)/denom*sign1*exp(log_I1-log_I4);
    I -=   l1*(l1-m+1.)*(l2+1.)*(l2+m)/denom*sign2*exp(log_I2-log_I4);
    I -=   (l1+1.)*(l1+m)*l2*(l2-m+1.)/denom*sign3*exp(log_I3-log_I4);
    I +=     l1*(l1-m+1.)*l2*(l2-m+1.)/denom*sign4;

    *sign = -MPOW(l2+1)*copysign(1,I);

    return log_B0+log_I4+log(fabs(I));
}

double casimir_integrate_C(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign)
{
    const int m = self->m;
    if(m == 0)
    {
        *sign = 0;
        return -INFINITY;
    }

    const double log_C0 = logi(m)+casimir_lnLambda(l1,l2,m)-self->tau;

    sign_t sign1, sign2;
    const double log_I1 = casimir_integrate_I(self, l1, l2-1, p, &sign1);
    const double log_I2 = casimir_integrate_I(self, l1, l2+1, p, &sign2);

    const double denom = 2.*l2+1.;
    double I;
    I  = -(l2+1.)*(l2+m)/denom*sign1*exp(log_I1-log_I2);
    I += l2*(l2-m+1.)/denom*sign2;

    *sign = -MPOW(l2)*copysign(1,I);

    return log_C0 + log_I2+log(fabs(I));
}

double casimir_integrate_D(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign)
{
    double log_C = casimir_integrate_C(self, l2, l1, p, sign);
    *sign *= MPOW(l1+l2+1);
    return log_C;
}
