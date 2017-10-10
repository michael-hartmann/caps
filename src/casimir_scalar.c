#include <omp.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "bessel.h"
#include "constants.h"
#include "plm.h"
#include "libcasimir.h"
#include "logfac.h"
#include "matrix.h"
#include "misc.h"
#include "quadpack.h"

#include "casimir_scalar.h"

/** @brief Initialize integration
 *
 * We want to compute the integral exp(-2*xi_*x)*Plm(l,2*m,x) from 1..oo.
 *
 * @param [in] m      magnetic quantum number
 * @param [in] xi_    xi*(R+L)/c
 * @param [in] ldim   size of vector space
 * @param [in] epsrel relative accuracy of integration
 * @retval self integration object
 */
integration_scalar_t *integration_init(int m, double xi_, int ldim, double epsrel)
{
    integration_scalar_t *self = xmalloc(sizeof(integration_scalar_t));

    self->xi_ = xi_;
    self->m   = m;
    self->K   = xmalloc(3*ldim*sizeof(double));
    self->epsrel = epsrel;

    /* initialize values with NANs */
    for(int i = 0; i < 3*ldim; i++)
        self->K[i] = NAN;

    return self;
}

/** @brief Free integration object
 *
 * @param [in] self integration object
 */
void integration_free(integration_scalar_t *self)
{
    if(self != NULL)
    {
        xfree(self->K);
        xfree(self);
    }
}

/** @brief logarithm of integrand
 *
 * The function returns log( exp(-2*xi_*x)*Plm(l,m,x) ) which is the logarithm
 * of the value of the integrand for parameters l,m,xi_ at x.
 *
 * @param [in] x position
 * @param [in] l azimuthal quantum number
 * @param [in] m magnetic quantum number
 * @param [in] xi_ xi*(R+L)/c
 * @retval logarithm of integrand
 */
static double log_integrand_K(double x, int l, int m, double xi_)
{
    return -2*xi_*x + Plm(l,m,x);
}

/** @brief Compute first and second derivative of integrand
 *
 * The integrand is given by exp(-2*xi_*x)*Plm(l,m,x). This function computes
 * the first and second derivative of the integrand.
 *
 * @param [in] x position
 * @param [in] l azimuthal quantum number
 * @param [in] m magnetic quantum number
 * @param [in] xi_ xi*(R+L)/c
 * @param [out] df first derivative of integrand
 * @param [out] d2f second derivative of integrand
 */
static void derivs(double x, int l, int m, double xi_, double *df, double *d2f)
{
    double p = 0, q = 0;
    const double c2 = (x+1)*(x-1);
    const double c  = sqrt(c2);

    if(m+1 <= l)
        p = 1/plm_continued_fraction(l,m+1,x);
    if(m+2 <= l)
        q = p/plm_continued_fraction(l,m+2,x);

    *df = -2*xi_ + p/c + m*x/c2;
    *d2f = (q - p*p - m*(x*x+1)/c2)/c2;
}

/** @brief Estimate position and width of maximum of integrand
 *
 * This function uses Newton's method to find the position of the maximum of
 * the integrand. The width is estimated using the second derivative at the
 * maximum. The value of the integral is estimated by approximating the maximum
 * with a Gaussian function and evaluating the integral analytical.
 *
 * @param [in] l azimuthal quantum number
 * @param [in] m magnetic quantum number
 * @param [in] xi_ xi*(R+L)/c
 * @param [out] Delta width of maximum
 * @param [out] I estimated value of integral
 * @retval xmax position of maximum
 */
static double estimate_integrand(int l, int m, double xi_, double *Delta, double *I)
{
    double df, d2f;
    const int maxiter = 20;
    double x = fmax(1.001,cosh(asinh((l+1)/(2*xi_)))-1);

    /* Newton's method */
    for(int i = 0; i < maxiter; i++)
    {
        const double xlast = x ;
        derivs(x,l,m,xi_, &df, &d2f);
        x = x - df/d2f;

        double w = 2/sqrt(-d2f);
        if(x < 1)
            x = 1+(xlast-1)/2;
        else if(fabs(x-xlast)/w < 1e-3)
            break;
    }

    derivs(x,l,m,xi_, &df, &d2f);
    const double f = log_integrand_K(x,l,m,xi_);
    *Delta = 2/sqrt(-d2f);
    *I = log(2*M_PI/-d2f)/2+f;

    return x;
}

/** @brief Integrand
 *
 * This is the function the integration function calls. It returns the value of
 * the integrand at position x, multiplies is by an factor of exp(-args->max)
 * and returns the result.
 *
 * @param [in] x position
 * @param [in] args_ pointer to struct containing parameters
 * @retval Kx scaled value of integrand
 */
static double integrand_K(double x, void *args_)
{
    integrand_t *args = (integrand_t *)args_;
    return exp(log_integrand_K(x,args->l,args->m,args->xi_)-args->max);
}

/** @brief Compute integral
 *
 * Compute the integral exp(-2*xi_*x)*Plm(l,2m,x) for x=1..oo.
 *
 * @param [in] l azimuthal quantum number
 * @param [in] m magnetic quantum number
 * @param [in] xi_ xi*(R+L)/c
 * @param [in] epsrel relative accuracy of integration
 * @retval K value of integral
 */
static double _integrate_K(int l, int m, double xi_, double epsrel)
{
    double width, Ie;
    const int nmax = 100;
    integrand_t args;

    /* exact solution for m = 0; 7.141 of Gradshteyn and Ryzhik */
    if(m == 0)
    {
        /* l = 0 is extremely simple */
        if(l == 0)
            /* int_0^infty dx sinh(x) exp(-2*xi_*cosh(x)) = exp(-2*xi_)/(2*xi_) */
            return -2*xi_-log(2*xi_);

        return -log(xi_*M_PI)/2 +bessel_lnKnu(l, 2*xi_);
    }

    /* position of maximum */
    const double xmax  = estimate_integrand(l, m, xi_, &width, &Ie);
    const double fxmax = log_integrand_K(xmax, l, m, xi_);

    args.l   = l;
    args.m   = m;
    args.xi_ = xi_;
    args.max = Ie;

    /* left and right boundary of integration */
    double a = fmax(1, xmax-8*width);
    double b = xmax+8*width;

    /* check left border */
    if(a > 1)
    {
        int i;
        for(i = 0; i < nmax; i++)
        {
            const double fa = log_integrand_K(a, l, m, xi_);
            if(exp(fa-fxmax) < epsrel)
                break;

            a = 1+0.5*(a-1);
        }

        TERMINATE(i == nmax, "l=%d, m=%d, xi_=%g, xmax=%g, f(xmax)=%g, a=%g", l, m, xi_, xmax, fxmax, a);
    }

    /* check right border */
    {
        int i;
        for(i = 0; i < nmax; i++)
        {
            const double fb = log_integrand_K(b, l, m, xi_);
            if(exp(fb-fxmax) < epsrel)
                break;

            b = 1+2*(b-1);
        }

        TERMINATE(i == nmax, "l=%d, m=%d, xi_=%g, xmax=%g, f(xmax)=%g, b=%g", l, m, xi_, xmax, fxmax, b);
    }

    /* perform integrations in interval [a,b] */
    int neval = 0, ier = 0;
    double abserr = 0;

    /* I: [a,b] */
    double I = dqags(integrand_K, a, b, 0, epsrel, &abserr, &neval, &ier, &args);

    TERMINATE(ier != 0, "ier=%d, l=%d, m=%d, xi_=%g, a=%g, b=%g, xmax=%.10g, width=%.10g, epsrel=%g, abserr=%g", ier, l, m, xi_, a, b, xmax, width, epsrel, abserr);

    return args.max+log(I);
}

/** @brief Get value of integral
 *
 * Get the value of the integral exp(-2*xi_*x)*Plm(l,m,x) for x=1..oo. This
 * function uses a cache to avoid recomputing integrals.
 *
 * @param [in] l azimuthal quantum number
 * @param [in] m magnetic quantum number
 * @param [in] xi_ xi*(R+L)/c
 * @param [in] epsrel relative accuracy of integration
 * @retval K value of integral
 */
static double integrate_K(integration_scalar_t *self, int l)
{
    int m = self->m;
    int index = l-m;
    double v = self->K[index];

    if(isnan(v))
        v = self->K[index] = _integrate_K(l,2*m,self->xi_,self->epsrel);

    return v;
}


/* eq. (3) */
static double alpha(double p, double n, double nu)
{
    return (((pow_2(p)-pow_2(n+nu+1))*(pow_2(p)-pow_2(n-nu)))/(4*pow_2(p)-1));
}

/** @brief Get value of integral
 *
 * Get the value of the integral exp(-2*xi_*x)*Plm(l1,m,x)*Plm(l2,m,x) for
 * x=1..oo.
 *
 * @param [in] self integration object
 * @param [in] l1 first azimuthal quantum number
 * @param [in] l2 second azimuthal quantum number
 * @retval I value of integral
 */
static double integrate_I(integration_scalar_t *self, int l1, int l2)
{
    double K;

    const int m_ = self->m;
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
    K = integrate_K(self, l1pl2-2*q);
    array[q].s = 1;
    array[q].v = K;

    if(qmax == 0)
        goto done;

    /* q = 1 */
    q = 1;
    /* eq. (29) */
    aq[q] = (n+nu-1.5)*(1-(2*n+2*nu-1)/(n4*(n4-1))*((m-n)*(m-n+1)/(2*n-1)+(m-nu)*(m-nu+1)/(2*nu-1)));
    K = integrate_K(self, l1pl2-2*q);
    array[q].s = SGN(aq[q]);
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

    K = integrate_K(self, l1pl2-2*q);
    array[q].s = SGN(aq[q]);
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

        K = integrate_K(self, l1pl2-2*q);
        array[q].s = SGN(aq[q]);
        array[q].v = log_scaling+K+log(fabs(aq[q]));

        if((array[q].v - array[0].v) < -75)
            break;
    }

    double log_I;
    sign_t sign;
    done:
    log_I = log_a0+logadd_ms(array, MIN(q,qmax)+1, &sign);
    TERMINATE(!isfinite(log_I), "l1=%d, l2=%d, m=%d, log_I=%g", l1, l2, m_, log_I);
    return log_I;
}

/** @brief Get reflection coefficient
 *
 * Get the reflection coefficient at the sphere for boundary condition Y.
 *
 * @param [in] casimir casimir object
 * @param [in] l azimuthal quantum number
 * @param [in] xi_ xi*(R+L)/c
 * @param [in] Y boundary condition (either D,N or R)
 * @retval rs reflection coefficient
 */
static double _rs(casimir_t *casimir, int l, double xi_, char Y)
{
    double lna, lnb;
    casimir_lnab_perf(casimir, xi_, l, &lna, &lnb);

    if(Y == 'N')
        return lna;
    else
        return lnb;
}

/** @brief Kernel of scalar Casimir effect
 *
 * Return matrix element of round-trip operator for scalar Casimir effect for
 * fixed magnetic quantum number m, ratio R/L and boundary conditions X and Y.
 *
 * @param i [in] azimuthal quantum number (row)
 * @param j [in] azimuthal quantum number (column)
 * @param args [in] pointer to structure containing the parameters
 * @param Mij matrix element
 */
static double kernel(int i, int j, void *args_)
{
    const args_t *args = (args_t *)args_;
    const int m = args->m;
    const int l1 = i+m, l2 = j+m;

    /* N_l1^m * N_l2^m */
    const double C = (logi(2*l1+1)+logi(2*l2+1)+lfac(l1-m)+lfac(l2-m)-lfac(l1+m)-lfac(l2+m))/2;

    if(isnan(args->rs[i]))
        args->rs[i] = _rs(args->casimir, l1, args->xi_, args->Y);
    if(isnan(args->rs[j]))
        args->rs[j] = _rs(args->casimir, l2, args->xi_, args->Y);

    const double rs = (args->rs[i]+args->rs[j])/2;
    const double I = integrate_I(args->integration, l1, l2);
    const double rp = args->X == 'D' ? 1 : -1;

    return rp*exp(C+rs+I);
}

/** @brief Compute logdetD for fixed magnetic quantum number m
 *
 * @param [in] LbyR L/R
 * @param [in] m magnetic quantum number
 * @param [in] xi_ xi*(L+R)/c
 * @param [in] ldim size of vector space
 * @param [in] X boundary condition at plate (D, N or R)
 * @param [in] Y boundary condition at sphere (D, N or R)
 * @param [in] epsrel relative accuracy of logdet value
 * @param logdetD
 */
double logdetD_m(double LbyR, double xi_, int m, int ldim, char X, char Y, double epsrel)
{
    args_t args;

    args.LbyR   = LbyR;
    args.xi_    = xi_;
    args.m      = m;
    args.X      = X;
    args.Y      = Y;
    args.epsrel = epsrel;
    args.integration = integration_init(m, xi_, 4*ldim, epsrel);
    args.casimir = casimir_init(LbyR);

    args.rs = xmalloc(ldim*sizeof(double));
    for(int i = 0; i < ldim; i++)
        args.rs[i] = NAN;

    double v = kernel_logdet(ldim, kernel, &args, 1, DETALG_HODLR);

    xfree(args.rs);
    integration_free(args.integration);
    casimir_free(args.casimir);

    return v;
}

/** @brief Compute logdetD
 *
 * @param [in] LbyR L/R
 * @param [in] xi_ xi*(L+R)/c
 * @param [in] ldim size of vector space
 * @param [in] X boundary condition at plate (D, N or R)
 * @param [in] Y boundary condition at sphere (D, N or R)
 * @param [in] epsrel relative accuracy of logdet value
 * @param logdetD
 */
double logdetD(double LbyR, double xi_, int ldim, char X, char Y, double epsrel)
{
    double values[1024] = { 0 };
    bool done = false;

    values[0] = logdetD_m(LbyR, xi_, 0, ldim, X, Y, epsrel);
    if(values[0] == 0)
        return 0;

    #pragma omp parallel for shared(done) schedule(static,1)
    for(size_t m = 1; m < sizeof(values)/sizeof(values[0]); m++)
    {
        if(!done)
        {
            double v = values[m] = logdetD_m(LbyR, xi_, m, ldim, X, Y, epsrel);
            if(v/values[0] < epsrel)
                done = true;
        }
    }

    double sum = values[0];
    for(size_t m = 1; m < sizeof(values)/sizeof(values[0]); m++)
        sum += 2*values[m];

     return sum;
}

/** @brief Integrand of xi
 *
 * @param [in] x position
 * @param [in] args_ pointer to structure with parameters
 */
static double integrand_xi(double x, void *args_)
{
    double xi_, v;
    integrand_xi_t *args = args_;

    xi_ = x*args->alpha;
    v = logdetD(args->LbyR, xi_, args->ldim, args->X, args->Y, args->epsrel);

    printf("# xi*R/c=%.10g, logdetD=%.10g\n", xi_, v);

    return v;
}

/** @brief Compute energy
 *
 * Compute energy for T=0, aspect ratio L/R, and boundary conditions X and Y.
 *
 * @param [in] LbyR L/R
 * @param [in] X boundary condition at plate (D, N or R)
 * @param [in] Y boundary condition at sphere (D, N or R)
 * @param [in] ldim dimension of vector space
 * @param [in] epsrel relative accuracy
 * @retval E energy
 */
double casimir_E(double LbyR, char X, char Y, int ldim, double epsrel)
{
    double abserr, I;
    int neval, ier;
    integrand_xi_t args;

    args.LbyR   = LbyR;
    args.ldim   = ldim;
    args.X      = X;
    args.Y      = Y;
    args.epsrel = epsrel/100;

    args.alpha = (1+LbyR)/(2*LbyR);

    I = dqagi(integrand_xi, 0, 1, 0, epsrel, &abserr, &neval, &ier, &args);

    TERMINATE(ier != 0, "ier=%d, LbyR=%g, X=%c, Y=%c, ldim=%d, epsrel=%g", ier, LbyR, X, Y, ldim, epsrel);
    return I*args.alpha/(2*M_PI);
}

/** @brief Print usage to stream
 *
 * @param self [in] self string with name of the program, i.e., argv[0]
 * @param stream [in] stream print message to stream
 */
static void usage(char *self, FILE *stream)
{
    fprintf(stream, "Usage: %s LbyR [bc ldim epsrel]\n\n", self);

    fprintf(stream, "Compute the Casimir energy for a scalar field with perfect conductors at T=0.\n\n");

    fprintf(stream, "Options:\n");
    fprintf(stream, "  LbyR:   aspect ratio L/R\n");
    fprintf(stream, "  bc:     boundary condition on plate and sphere (DD, DN, ND or NN)\n");
    fprintf(stream, "  ldim:   dimension of vector space\n");
    fprintf(stream, "  epsrel: relative accuracy of integration\n\n");

    fprintf(stream, "The number of threads can be set using OMP_NUM_THREADS.\n\n");

    casimir_build(stream, NULL);
}

int main(int argc, char *argv[])
{
    double LbyR, epsrel = 5e-7;
    int eta = 10, ldim = -1;
    char X = 'D', Y = 'D';

    disable_buffering();

    if(argc < 2)
    {
        usage(argv[0], stderr);
        return 1;
    }

    LbyR = atof(argv[1]);
    if(LbyR <= 0)
    {
        fprintf(stderr, "LbyR must be positive\n\n");
        usage(argv[0], stderr);
        return 1;
    }

    if(argc > 2)
    {
        for(size_t i = 0; i < strlen(argv[2]); i++)
            argv[2][i] = toupper(argv[2][i]);

        if(strcmp(argv[2], "DD") != 0 && strcmp(argv[2], "DN") && strcmp(argv[2], "ND") && strcmp(argv[2], "NN"))
        {
            fprintf(stderr, "boundary condition must be either DD, DN, ND or NN\n\n");
            usage(argv[0], stderr);
            return 1;
        }
        X = toupper(argv[2][0]);
        Y = toupper(argv[2][1]);
    }

    if(argc > 3)
    {
        ldim = atoi(argv[3]);
        if(ldim <= 0)
        {
            fprintf(stderr, "ldim must be positive\n\n");
            usage(argv[0], stderr);
            return 1;
        }
    }
    else
        ldim = eta/LbyR;

    if(argc > 4)
    {
        epsrel = atof(argv[4]);
        if(epsrel <= 0)
        {
            fprintf(stderr, "epsrel must be positive\n\n");
            usage(argv[0], stderr);
            return 1;
        }
    }

    printf("# LbyR   = %g\n", LbyR);
    printf("# X      = %c\n", X);
    printf("# Y      = %c\n", Y);
    printf("# ldim   = %d\n", ldim);
    printf("# epsrel = %g\n", epsrel);
    printf("# cores  = %d\n", omp_get_max_threads());
    printf("#\n");

    const double E = casimir_E(LbyR, X, Y, ldim, epsrel);

    printf("\n");
    printf("# L/R, X, Y, ldim, E*(L+R)/(hbar*c)\n");
    printf("%.15g, %c, %c, %d, %.15g\n", LbyR, X, Y, ldim, E);

    return 0;
}
