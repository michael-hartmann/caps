/**
 * @file   integration.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   April, 2017
 * @brief  Perform integration for arbitrary materials
 */

#include <math.h>
#include <stdint.h>
#include <inttypes.h>

#include "quadpack.h"

#include "utils.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"
#include "hash-table.h"

/* arguments for integrand in function K_integrand */
typedef struct
{
    int nu,m;
    polarization_t p; /* TE or TM */
    double factor,tau,log_normalization;
    casimir_t *casimir;
} integrand_t;

/* entry of cache consists of value of integral and sign: sign*exp(v) */
typedef struct
{
    double v;    /**< logarithmic value of integral */
    sign_t sign; /**< sign of integral */
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
static uint64_t hash(uint64_t l1, uint64_t l2, uint64_t p)
{
    return (l1 << 32) | (l2 << 1) | p;
}


/** @brief Bisect method for f(x) = c
 *
 * Apply the bisection method for
 *      f(x) = x^ν*exp(-τx) -c = 0,
 * where c = exp(log_c).
 *
 * @param [in] left left boundary
 * @param [in] right right boundary
 * @param [in] N number of times to apply bisection method
 * @param [in] nu parameter nu
 * @param [in] tau parameter tau
 * @param [in] log_c logarithm of parameter c
 * @retval x0 root of f(x)
 */
static double _f_bisect(double left, double right, int N, int nu, double tau, double log_c)
{
    for(int i = 0; i < N; i++)
    {
        const double middle = (right+left)/2;
        const double fl = nu*log(left)  -tau*left  -log_c;
        const double fm = nu*log(middle)-tau*middle-log_c;

        if(fl*fm < 0)
            right = middle;
        else
            left = middle;
    }

    return (left+right)/2;
}

/** @brief Analyze function f(x) = x^ν*exp(-τx)
 *
 * This function analyzes the function
 *      f(x) = x^ν*exp(-τx).
 *
 * At x=0 the function is zero,
 *      f(0) = 0,
 * then it increases up to the unique maximum at
 *      x_max = ν/τ,
 * and decreases for x > x_max.
 *
 * This function returns the logarithm of the maximum, i.e.
 *      log(f(x_max)).
 *
 * There exist two values x* with
 *      f(x*) = f(x_max)*eps.
 * The value x* in (0,x_max) will be stored in left, the value x* in
 * (x_max,infty) will be stored in right.
 *
 * @param [in] nu parameter nu
 * @param [in] tau parameter tau
 * @param [in] eps determines left and right boundaries (see description)
 * @param [in] tol accuracy of left and right
 * @param [out] left point in (0,x_max) with f(x_max)*eps = f(left)
 * @param [out] right point in (x_max,infty) with f(x_max)*eps = f(right)
 * @retval maximum logarithm of maximum, log(f(x_max))
 */
static double _f_estimate(int nu, double tau, double eps, double tol, double *left, double *right)
{
    /* f(x) = x^ν*exp(-τx) */
    const double x_max = nu/tau;
    const double log_f_max = nu*log(x_max)-tau*x_max; /* log(f(x_max)) */

    const int N = ceil(-log(tol/x_max)/M_LOG2);
    *left = _f_bisect(0, x_max, N, nu, tau, log_f_max+log(eps));

    for(int i = 1; true; i++)
    {
        double x = (i+1)*x_max;
        double log_f = nu*log(x) - tau*x;

        if((log_f-log_f_max) < log(eps))
        {
            *right = _f_bisect(i*x_max, (i+1)*x_max, N, nu, tau, log_f_max+log(eps));
            return log_f_max;
        }
    }
}

/* Estimate shape of integrand K assuming nu/tau >> 1.
 *
 * If nu/tau >> 1, we can approximate the associated Legendre polynomial
 *      Plm(ν,m,1+z) ≈ (2l)!/(2^l*l!*(l-m)!) z^ν.
 *
 * We are only interested in the rough shape of the integrand, so we also
 * ignore the Fresnel coefficient r_p. This means, we assume that the Fresnel
 * coefficient varies slowly compared to the other factors of the integrand.
 *
 * Thus the integrand looks like (m>0)
 *      k(z) ≈ (2l)!/(2^l*l!*(l-m)!) z^(ν-2) * exp(-τz).
 *
 * The maximum of the integrand is at
 *      zmax ≈ (ν-2)/τ.
 *
 * For m=0 we find
 *      k(z) = (2l)!/(2^l*l!*l!) z^ν * exp(-τz).
 * and
 *      zmax ≈ ν/τ.
 *
 * We estimate the position a < zmax where the integrand is approximately
 * zmax*eps with precision tol=1e-2 using the bisection method for the
 * simplified integrand k(z).
 *
 * We also estimate the position b > zmax where the integrand is approximately
 * zmax*eps with precision tol=1e-2 using bisection method.
 */
static double K_estimate_zlarge(int nu, int m, double tau, double eps, double *a, double *b, double *log_normalization)
{
    const double tol = 1e-2;
    double zmax;

    if(m > 0)
    {
        zmax = (nu-2)/tau;
        *log_normalization = Plm_estimate(nu,2*m,1+zmax)-log(zmax*(zmax+2))-tau*zmax;

        nu -= 2;
    }
    else
    {
        zmax = nu/tau;
        *log_normalization = Plm_estimate(nu,0,1+zmax)-tau*zmax;
    }

    /* k(z) = z^ν*exp(-τz) */
    _f_estimate(nu, tau, eps, tol, a, b);

    return zmax;
}

static double log_k(double z, int nu, int m, double tau)
{
    if(z == 0)
        return -INFINITY;

    if(m == 0)
        return Plm(nu,2,1+z) - tau*z;
    else
        return Plm(nu,2*m,1+z)-log(z*(z+2)) - tau*z;
}

static double K_estimate_zsmall(int nu, int m, double tau, double eps, double *a, double *b, double *log_normalization)
{
    double left, right, middle, log_k_zmax, zmax;
    const double tol = 1e-4;

    *a = 0;
    *b = 4;

    /* Now, we search a good approximation for zmax using Golden section
     * search.
     * The code for Golden section search is taken from Wikipedia:
     * https://en.wikipedia.org/wiki/Golden_section_search
     */
    double fc = NAN, fd = NAN;
    double c = *b-(*b-*a)/M_GM;
    double d = *a+(*b-*a)/M_GM;

    while(fabs(c-d) > tol)
    {
        /* we only have to caclulate fc or fd if it is NAN */
        if(isnan(fc))
            fc = log_k(c, nu, m, tau);
        if(isnan(fd))
            fd = log_k(d, nu, m, tau);

        if(fc > fd)
        {
            *b = d;
            fd = fc;
            fc = NAN;
        }
        else
        {
            *a = c;
            fc = fd;
            fd = NAN;
        }

        /* we recompute both c and d here to avoid loss of precision which may
         * lead to incorrect results or infinite loop */
        c = *b-(*b-*a)/M_GM;
        d = *a+(*b-*a)/M_GM;
    }

    /* now we have a good approximation for zmax */
    zmax = (*a+*b)/2;

    /* and we also know log(k(zmax)) which is approximately fc or fd */
    if(!isnan(fc))
        log_k_zmax = fc;
    else
        log_k_zmax = fd;

    /* thus we can set log_normalization */
    *log_normalization = log_k(zmax, nu,m,tau);

    /* now we are trying to estimate the interval [a,b] which gives the main
     * contributions to the integration; we are looking for a < zmax < b such
     * that a/zmax ≈ b/zmax ≈ eps
     *
     * We are using the bisection method and have to apply it N times.
     */
    const int N = ceil(-log(tol)/M_LOG2);

    /* a */
    left = 0;
    right = zmax;
    for(int i = 0; i < N; i++)
    {
        middle = (left+right)/2;
        double fm = log_k(middle, nu, m, tau);

        if((fm-log_k_zmax) < log(eps))
            left = middle;
        else
            right = middle;
    }
    *a = (left+right)/2;

    if(*a < 1e-5)
        *a = 0;

    /* b */
    left = zmax;
    right = 4;
    for(int i = 0; i < N; i++)
    {
        middle = (left+right)/2;
        double fm = log_k(middle, nu, m, tau);

        if((fm-log_k_zmax) < log(eps))
            right = middle;
        else
            left = middle;
    }
    *b = (left+right)/2;

    //printf("a=%g (%g), zmax=%g (%g), b=%g (%g) (%g)\n", *a, log_k(*a,nu,m,tau), zmax, log_k(zmax,nu,m,tau), *b, log_k(*b,nu,m,tau), *log_normalization*nu);
    //printf("log_k = %g\n", log_k(0.00258133,nu,m,tau));

    return zmax;
}

static double K_integrand(double z, void *args_)
{
    double v, rTE, rTM;
    integrand_t *args = (integrand_t *)args_;
    casimir_t *casimir = args->casimir;

    z *= args->factor;

    const int nu = args->nu, m = args->m;
    const double log_normalization = args->log_normalization;
    const double tau = args->tau;
    const double xi = tau/2;
    const double z2p2z = z*(z+2);

    if(m)
        v = exp(-log_normalization + Plm(nu,2*m,1+z)-tau*z)/z2p2z;
    else
        v = exp(-log_normalization + Plm(nu,2,1+z)-tau*z);

    casimir->rp(casimir, xi, xi*sqrt(z2p2z), &rTE, &rTM);

    TERMINATE(isnan(v) || isinf(v), "z=%g, nu=%d, m=%d, tau=%g, v=%g", z, nu, m, tau, v);

    if(args->p == TE)
        return rTE*v;
    else
        return rTM*v;
}

static double _casimir_integrate_K(integration_t *self, int nu, polarization_t p, sign_t *sign)
{
    double zmax,log_normalization,a,b;
    const int m = self->m;
    const double eps = 1e-6;
    const double tau = self->tau;
    const double epsrel = self->epsrel;

    integrand_t args = {
        .nu   = nu,
        .m    = m,
        .p    = p,
        .tau  = tau,
        .factor = 1,
        .casimir = self->casimir
    };

    /* estimate width and height of the peak of the integrand. */
    if(nu == 2 && m == 1)
    {
        /* Integrate
         *      r_p*Plm(2,2,1+z)/(z²+2z)*exp(-τz) = -3*r_p*exp(-τz)
         *
         * Maximum is at z=0.
         */
        zmax = 0;
        a = 0;
        b = -log(eps)/tau; /* exp(-14) =~ 1e-6 */
        log_normalization = log(3);
    }
    else if(nu/tau > 4)
        /* large z limit */
        zmax = K_estimate_zlarge(nu, m, tau, eps, &a, &b, &log_normalization);
    else
        /* small z limit */
        zmax = K_estimate_zsmall(nu, m, tau, eps, &a, &b, &log_normalization);

    /* log_normalization is the logarithm of the estimated maximum of the
     * integrand K
     */
    args.log_normalization = log_normalization;

    /* perform integrations in intervals [0,a], [a,b] and [b,∞] */
    int neval1 = 0, neval2 = 0, neval3 = 0, ier1 = 0, ier2 = 0, ier3 = 0;
    double abserr1 = 0, abserr2 = 0, abserr3 = 0, I1 = 0, I2 = 0, I3 = 0;

    /* I2: [a,b] */
    I2 = dqags(K_integrand, a, b, 0, epsrel, &abserr2, &neval2, &ier2, &args);

    /* I1: [0,a] */
    if(a > 0)
    {
        /* The contribution of this integral should be small, so use
         * Gauss-Kronrod G_K 7-15 as integration rule.
         */
        int limit = 200;
        I1 = dqage(K_integrand, 0, a, abserr2, epsrel, GK_7_15, &abserr1, &neval1, &ier1, &limit, &args);
    }

    /* I3: [b,∞]
     * Make a substitution for the integrand, i.e., we don't integrate over z
     * but over t=z*τ. The integrand exponentially decays as exp(-z*τ), so
     * after the substitution it decays as exp(-t). This makes life easier for
     * the quadrature routine.
     */
    args.factor = 1/tau;
    I3 = dqagi(K_integrand, b*tau, 1, abserr2*tau, epsrel, &abserr3, &neval3, &ier3, &args)/tau;

    const double sum = I1+I2+I3;

    bool warn = ier1 != 0 || ier2 != 0 || ier3 != 0 || isnan(sum) || sum == 0;
    WARN(warn, "ier1=%d, ier2=%d, ier3=%d, nu=%d, m=%d, tau=%.20g, zmax=%g, a=%g, b=%g, I1=%g, I2=%g, I3=%g", ier1, ier2, ier3, nu,m,tau,zmax,a,b, I1, I2, I3);

    *sign = SGN(sum);
    return log(fabs(sum)) + log_normalization;
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

/* eq. (3) */
static double alpha(double p, double n, double nu)
{
    return (((pow_2(p)-pow_2(n+nu+1))*(pow_2(p)-pow_2(n-nu)))/(4*pow_2(p)-1));
}

static double _casimir_integrate_I(integration_t *self, int l1, int l2, polarization_t p_, sign_t *sign)
{
    sign_t s;
    double K;

    const int m_ = self->m > 0 ? self->m : 1;
    const double n  = l1;
    const double nu = l2;
    const double m  = m_;
    const double n4 = l1+l2-2*m_;
    const int l1pl2 = l1+l2;

    /* eq. (24) */
    const int qmax = MIN(MIN(n,nu), (l1pl2-2*m_)/2);

    TERMINATE(qmax < 0, "l1=%d, l2=%d, m=%d\n", l1, l2, (int)m);

    /* eq. (28) */
    const double Ap = -2*m*(n-nu)*(n+nu+1);

    /* eq. (20) */
    const double log_a0 = lfac(2*l1)-lfac(l1)+lfac(2*l2)-lfac(l2)+lfac(l1+l2)-lfac(2*l1pl2)+lfac(l1pl2-2*m_)-lfac(l1-m_)-lfac(l2-m_);

    double aq[qmax+1];
    log_t array[qmax+1];


    /* q = 0 */
    int q = 0;
    aq[q] = 1;
    K = casimir_integrate_K(self, l1pl2-2*q, p_, &s);
    array[q].s = s;
    array[q].v = K;

    if(qmax == 0)
        goto done;

    /* q = 1 */
    q = 1;
    /* eq. (29) */
    aq[q] = (n+nu-1.5)*(1-(2*n+2*nu-1)/(n4*(n4-1))*((m-n)*(m-n+1)/(2*n-1)+(m-nu)*(m-nu+1)/(2*nu-1)));
    K = casimir_integrate_K(self, l1pl2-2*q, p_, &s);
    array[q].s = SGN(aq[q])*s;
    array[q].v = K+log(fabs(aq[q]));
    if(qmax == 1)
        goto done;

    q = 2;
    /* eq. (35) */
    aq[q] = (2*n+2*nu-1)*(2*n+2*nu-7)/4*( (2*n+2*nu-3)/(n4*(n4-1)) * ( (2*n+2*nu-5)/(2*(n4-2)*(n4-3)) \
                * ( (m-n)*(m-n+1)*(m-n+2)*(m-n+3)/(2*n-1)/(2*n-3) \
                + 2*(m-n)*(m-n+1)*(m-nu)*(m-nu+1)/((2*n-1)*(2*nu-1)) \
                + (m-nu)*(m-nu+1)*(m-nu+2)*(m-nu+3)/(2*nu-1)/(2*nu-3) ) - (m-n)*(m-n+1)/(2*n-1) \
                - (m-nu)*(m-nu+1)/(2*nu-1) ) +0.5);

    K = casimir_integrate_K(self, l1pl2-2*q, p_, &s);
    array[q].s = SGN(aq[q])*s;
    array[q].v = K+log(fabs(aq[q]));

    if(qmax == 2)
        goto done;

    double log_scaling = 0;
    for(q = 3; q <= qmax; q++)
    {
        const double p = n+nu-2*q;
        const double p1 = p-2*m;
        const double p2 = p+2*m;

        if(Ap != 0)
        {
            /* eqs. (26), (27) */
            double c0 = (p+2)*(p+3)*(p1+1)*(p1+2)*Ap*alpha(p+1,n,nu);
            double c1 = Ap*(Ap*Ap \
               + (p+1)*(p+3)*(p1+2)*(p2+2)*alpha(p+2,n,nu) \
               + (p+2)*(p+4)*(p1+3)*(p2+3)*alpha(p+3,n,nu));
            double c2 = -(p+2)*(p+3)*(p2+3)*(p2+4)*Ap*alpha(p+4,n,nu);

            aq[q] = (c1*aq[q-1] + c2*aq[q-2])/c0;
        }
        else
            /* eq. (30) */
            aq[q] = (p+1)*(p2+2)*alpha(p+2,n,nu)*aq[q-1] / ((p+2)*(p1+1)*alpha(p+1,n,nu));

        if(fabs(aq[q]) > 1e100)
        {
            log_scaling += log(fabs(aq[q]));
            aq[q-1] /= fabs(aq[q]);
            aq[q] = SGN(aq[q]);
        }
        else if(fabs(aq[q]) < 1e-100 && fabs(aq[q]) > 0)
        {
            log_scaling -= log(fabs(aq[q]));
            aq[q-1] *= fabs(aq[q]);
            aq[q] = SGN(aq[q]);
        }

        K = casimir_integrate_K(self, l1pl2-2*q, p_, &s);
        array[q].s = SGN(aq[q])*s;
        array[q].v = log_scaling+K+log(fabs(aq[q]));

        if((array[q].v - array[0].v) < -75)
            break;
    }

    double log_I;
    done:
    log_I = log_a0+logadd_ms(array, MIN(q,qmax)+1, sign);
    TERMINATE(!isfinite(log_I), "l1=%d, l2=%d, p=%d, log_I=%g", l1, l2, p_, log_I);
    return log_I;
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


/** @brief Initialize integration
 *
 * Initialize integration for (scaled) Matsubara frequency xi and quantum
 * number m. The aspect ratio L/R and the dielectric function of the metals is
 * taken from casimir. The integration is performed to a relative accuracy of
 * epsrel.
 *
 * This function returns an object in order to compute the actual integrals.
 * This memory of this object has to be freed after use by a call to \ref
 * casimir_integrate_free.
 *
 * @param [in] casimir Casimir object
 * @param [in] xi (scaled) Matsubara frequency
 * @param [in] m quantum number
 * @param [in] epsrel relative accuracy of integration
 * @retval integration object
 */
integration_t *casimir_integrate_init(casimir_t *casimir, double xi, int m, double epsrel)
{
    if(xi < 0 || m < 0 || epsrel <= 0)
        return NULL;

    integration_t *self = (integration_t *)xmalloc(sizeof(integration_t));

    self->casimir = casimir;
    self->m = m;
    self->tau = 2*xi;
    self->epsrel = epsrel;

    self->hash_table_I = hash_table_new(cache_entry_destroy);
    self->hash_table_K = hash_table_new(cache_entry_destroy);

    if(isinf(casimir_epsilonm1(casimir, INFINITY)))
        self->is_pc = true;
    else
        self->is_pc = false;

    return self;
}


/** @brief Free integration object
 *
 * @param [in,out] integration integration object
 */
void casimir_integrate_free(integration_t *integration)
{
    if(integration != NULL)
    {
        hash_table_free(integration->hash_table_I);
        hash_table_free(integration->hash_table_K);
        xfree(integration);
    }
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
    const double A0 = 2*logi(m)+casimir_lnLambda(l1,l2,m)-self->tau;

    *sign *= -MPOW(l2);

    const double A = A0+I1;
    TERMINATE(!isfinite(A), "l1=%d, l2=%d, m=%d, p=%d, I1=%g, A0=%g, A=%g", l1,l2,m,p, I1, A0, A);
    return A;
}

double casimir_integrate_B(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign)
{
    const int m = self->m;

    const double B0 = casimir_lnLambda(l1,l2,m)-self->tau;

    if(m == 0)
    {
        const double I = casimir_integrate_I(self, l1, l2, p, sign);
        const double B = B0+I;

        TERMINATE(!isfinite(B), "l1=%d, l2=%d, m=%d, p=%d, I=%g, B0=%g, B=%g", l1,l2,m,p, I, B0, B);

        *sign *= -MPOW(l2+1);

        return B;
    }

    sign_t sign1, sign2, sign3, sign4;
    const double I1 = casimir_integrate_I(self, l1-1, l2-1, p, &sign1);
    const double I2 = casimir_integrate_I(self, l1+1, l2-1, p, &sign2);
    const double I3 = casimir_integrate_I(self, l1-1, l2+1, p, &sign3);
    const double I4 = casimir_integrate_I(self, l1+1, l2+1, p, &sign4);

    double I;
    const double denom = (2*l1+1.)*(2*l2+1.);
    I  = (l1+1.)*(l1+m)*(l2+1.)*(l2+m)/denom*sign1*exp(I1-I4);
    I -=   l1*(l1-m+1.)*(l2+1.)*(l2+m)/denom*sign2*exp(I2-I4);
    I -=   (l1+1.)*(l1+m)*l2*(l2-m+1.)/denom*sign3*exp(I3-I4);
    I +=     l1*(l1-m+1.)*l2*(l2-m+1.)/denom*sign4;

    *sign = -MPOW(l2+1)*SGN(I);

    const double B = B0+I4+log(fabs(I));
    TERMINATE(!isfinite(B), "l1=%d, l2=%d, m=%d, p=%d, I1=%g, I2=%g, I3=%g, I4=%g, B0=%g, B=%g", l1,l2,m,p, I1,I2,I3,I4, B0, B);
    return B;
}

double casimir_integrate_C(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign)
{
    const int m = self->m;
    if(m == 0)
    {
        *sign = 0;
        return -INFINITY;
    }

    const double C0 = logi(m)+casimir_lnLambda(l1,l2,m)-self->tau;

    sign_t sign1, sign2;
    const double I1 = casimir_integrate_I(self, l1, l2-1, p, &sign1);
    const double I2 = casimir_integrate_I(self, l1, l2+1, p, &sign2);

    const double denom = 2*l2+1;
    double I;
    I  = -(l2+1.)*(l2+m)/denom*sign1*exp(I1-I2);
    I += l2*(l2-m+1.)/denom*sign2;

    *sign = -MPOW(l2)*SGN(I);

    const double C = C0+I2+log(fabs(I));
    TERMINATE(!isfinite(C), "l1=%d, l2=%d, m=%d, p=%d, I1=%g, I2=%g, C0=%g, C=%g", l1,l2,m,p, I1,I2, C0, C);
    return C;
}

double casimir_integrate_D(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign)
{
    double log_C = casimir_integrate_C(self, l2, l1, p, sign);
    *sign *= MPOW(l1+l2+1);
    return log_C;
}

integration_plasma_t *casimir_integrate_plasma_init(double omegap, double epsrel)
{
    integration_plasma_t *self;
    self = (integration_plasma_t *)xmalloc(sizeof(integration_plasma_t));
    self->omegap = omegap;
    self->epsrel = epsrel;
    self->cache = hash_table_new(cache_entry_destroy);

    return self;
}

typedef struct {
    int nu;
    double omegap, log_prefactor;
} integrand_plasma_t;

static double _integrand_plasma(double z, void *args_)
{
    integrand_plasma_t *args = (integrand_plasma_t *)args_;

    const int nu               = args->nu;
    const double omegap        = args->omegap;
    const double log_prefactor = args->log_prefactor;

    const double k = 0.5*z;
    const double beta = sqrt(1+pow_2(omegap/k));
    const double rTE = (1-beta)/(1+beta);

    return -rTE * exp(log_prefactor -z+nu*log(z));
}

double casimir_integrate_plasma(integration_plasma_t *self, double LbyR, int l1, int l2, int m)
{
    const int nu = l1+l2;
    const double y = -M_LOG2-log1p(LbyR);
    double I;
    cache_entry_t *entry = hash_table_lookup(self->cache, nu);

    if(entry == NULL)
    {
        const double epsrel = self->epsrel;

        /* compute integral
         *      1/prefactor * int_0^infty dz r_TE e^(-z) z^nu
         * with
         *      prefactor = nu!
         */
        const double log_prefactor = -lfac(nu);

        /* find left and right boundaries */
        double a,b;
        const double tol = 1e-2;
        const double eps = 1e-6;
        _f_estimate(nu, 1, eps, tol, &a, &b);

        /* perform integrations in intervals [0,a], [a,b] and [b,∞] */
        integrand_plasma_t args = { .nu = nu, .omegap = self->omegap, .log_prefactor = log_prefactor };

        int neval1 = 0, neval2 = 0, neval3 = 0, ier1 = 0, ier2 = 0, ier3 = 0;
        double abserr1 = 0, abserr2 = 0, abserr3 = 0, I1 = 0, I2 = 0, I3 = 0;

        /* I2: [a,b] */
        I2 = dqags(_integrand_plasma, a, b, 0, epsrel, &abserr2, &neval2, &ier2, &args);

        /* I1: [0,a] */
        if(a > 0)
        {
            /* The contribution of this integral should be small, so use
             * Gauss-Kronrod G_K 7-15 as integration rule.
             */
            int limit = 200;
            I1 = dqage(_integrand_plasma, 0, a, abserr2, epsrel, GK_7_15, &abserr1, &neval1, &ier1, &limit, &args);
        }

        /* I3: [b,∞] */
        I3 = dqagi(_integrand_plasma, b, 1, abserr2, epsrel, &abserr3, &neval3, &ier3, &args);

        I = I1+I2+I3;
        bool warn = ier1 != 0 || ier2 != 0 || ier3 != 0 || isnan(I) || I == 0;
        WARN(warn, "ier1=%d, ier2=%d, ier3=%d, nu=%d, m=%d, a=%g, b=%g, I1=%g, I2=%g, I3=%g", ier1, ier2, ier3, nu,m,a,b, I1, I2, I3);

        entry = cache_entry_create(I, +1);
        hash_table_insert(self->cache, nu, entry);
    }
    else
        I = entry->v;

    return exp( lfac(nu)-0.5*(lfac(l1+m)+lfac(l1-m)+lfac(l2+m)+lfac(l2-m)) + (nu+1)*y )*sqrt((l1*l2)/((l1+1.)*(l2+1.)))*I;
}

void casimir_integrate_plasma_free(integration_plasma_t *self)
{
    hash_table_free(self->cache);
    xfree(self);
}
