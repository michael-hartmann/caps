#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "bessel.h"
#include "plm.h"
#include "logfac.h"
#include "matrix.h"
#include "misc.h"
#include "quadpack.h"

#include "casimir_scalar.h"

integration_t *integration_init(int m, double xi_, int ldim, double epsrel)
{
    integration_t *self = xmalloc(sizeof(integration_t));

    self->xi_ = xi_;
    self->m   = m;
    self->epsrel = epsrel;
    self->K  = xmalloc(3*ldim*sizeof(double));

    for(int i = 0; i < 3*ldim; i++)
        self->K[i] = NAN;

    return self;
}

void integration_free(integration_t *self)
{
    if(self != NULL)
    {
        xfree(self->K);
        xfree(self);
    }
}

static double log_integrand_K(double x, int l, int m, double xi_)
{
    const double log_sinhx = log(sinh(x));
    const double coshx = cosh(x);

    return log_sinhx -(2*xi_)*coshx + Plm(l,m,coshx);
}

static double integrand_K(double x, void *args_)
{
    integrand_t *args = (integrand_t *)args_;
    return exp(log_integrand_K(x,args->l,args->m,args->xi_)-args->max);
}

static double _integrate_K(int l, int m, double xi_, double epsrel)
{
    integrand_t args;

    double xmax  = asinh((l+1)/(2*xi_));
    double width = 1/sqrt(2*xi_*cosh(xmax));

    args.l   = l;
    args.m   = m;
    args.xi_ = xi_;
    args.max = log_integrand_K(xmax, l, m, xi_);

    double a = fmax(0, xmax-8*width);
    double b = xmax+8*width;

    /* perform integrations in interval [a,b] */
    int neval = 0, ier = 0;
    double abserr = 0;

    /* I2: [a,b] */
    double I = dqags(integrand_K, a, b, 0, epsrel, &abserr, &neval, &ier, &args);

    #if 0
    printf("l1=%d, l2=%d, m=%d, tau=%g\n", l1, l2, m, tau);
    printf("width=%g\n", width);
    printf("xmax=%g, max=%g\n", xmax, args.max);
    printf("a=%g, b=%g\n", a, b);
    printf("I=%g\n", I);
    #endif

    return args.max+log(I);
}

static double integrate_K(integration_t *self, int l)
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

static double integrate_I(integration_t *self, int l1, int l2)
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

static double _rs(int l, double chi, char Y)
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

    if(isnan(args->rs[i]))
        args->rs[i] = _rs(l1, xi_/(1+LbyR), args->Y);
    if(isnan(args->rs[j]))
        args->rs[j] = _rs(l2, xi_/(1+LbyR), args->Y);

    rs = (args->rs[i]+args->rs[j])/2;

    I = integrate_I(args->integration, l1, l2);

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
    args.integration = integration_init(m, xi_, 4*ldim, epsrel);

    args.rs = xmalloc(ldim*sizeof(double));
    for(int i = 0; i < ldim; i++)
        args.rs[i] = NAN;

    double v = kernel_logdet(ldim, kernel, &args, 1, DETALG_HODLR);

    xfree(args.rs);
    integration_free(args.integration);

    return v;
}

int main(int argc, char *argv[])
{
    const double epsrel = 1e-6;

    if(argc < 3)
    {
        fprintf(stderr, "usage: %s LbyR xi_ m\n", argv[0]);
        return 1;
    }

    const double LbyR = atof(argv[1]);
    const int m       = atoi(argv[2]);
    const double xi_  = atof(argv[3]);

    double v = logdetD(LbyR, xi_, m, 6/LbyR, 'D', 'D', epsrel);

    printf("LbyR    = %g\n", LbyR);
    printf("m       = %d\n", m);
    printf("xi_     = %g\n", xi_);
    printf("logdetD = %.10g\n", v);

    return 0;
}
