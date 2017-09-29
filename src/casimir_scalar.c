#include <omp.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "constants.h"
#include "plm.h"
#include "libcasimir.h"
#include "logfac.h"
#include "matrix.h"
#include "misc.h"
#include "quadpack.h"

#include "casimir_scalar.h"

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

void integration_free(integration_scalar_t *self)
{
    if(self != NULL)
    {
        xfree(self->K);
        xfree(self);
    }
}

static double log_integrand_K(double x, int l, int m, double xi_)
{
    const double log_sinhx = log(0.5)+x+log1p(-exp(-2*x)); /* log(sinh(x)) */
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

    /* position of maximum */
    const double xmax  = asinh((l+1)/(2*xi_));
    /* width of maximum */
    const double width = 1/sqrt(2*xi_*cosh(xmax));

    args.l   = l;
    args.m   = m;
    args.xi_ = xi_;
    args.max = log_integrand_K(xmax, l, m, xi_);

    /* left and right boundary of integration */
    const double a = fmax(0, xmax-8*width);
    const double b = xmax+8*width;

    /* perform integrations in interval [a,b] */
    int neval = 0, ier = 0;
    double abserr = 0;

    /* I2: [a,b] */
    double I = dqags(integrand_K, a, b, 0, epsrel, &abserr, &neval, &ier, &args);

    TERMINATE(ier != 0, "ier=%d, l=%d, m=%d, xi_=%g, epsrel=%g", ier, l, m, xi_, epsrel);

    return args.max+log(I);
}

/* Integrate
 *      sinh(x) * exp(-2*xi_*cosh(x)) * P(l,2m,cosh(x))
 * from x=0...oo
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

/* Integrate
 *      sinh(x) * exp(-2*xi_*cosh(x)) * P(l1,m,cosh(x)) * P(l2,m,cosh(x))
 * from x=0...oo
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

static double _rs(casimir_t *casimir, int l, double xi_, char Y)
{
    double lna, lnb;
    casimir_lnab_perf(casimir, xi_, l, &lna, &lnb);

    if(Y == 'N')
        return lna;
    else
        return lnb;
}

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

static double integrand_xi(double x, void *args_)
{
    double xi_, v;
    integrand_xi_t *args = args_;

    xi_ = x*args->alpha;
    v = logdetD(args->LbyR, xi_, args->ldim, args->X, args->Y, args->epsrel);

    printf("# xi*R/c=%.10g, logdetD=%.10g\n", xi_, v);

    return v;
}

double casimir_E(double LbyR, char X, char Y, int ldim, double epsrel)
{
    double abserr;
    int neval, ier;
    integrand_xi_t args;

    args.LbyR   = LbyR;
    args.ldim   = ldim;
    args.X      = X;
    args.Y      = Y;
    args.epsrel = epsrel/100;

    args.alpha = (1+LbyR)/(2*LbyR);

    double I = dqagi(integrand_xi, 0, 1, 0, epsrel, &abserr, &neval, &ier, &args)*args.alpha/(2*M_PI);

    TERMINATE(ier != 0, "ier=%d, LbyR=%g, X=%c, Y=%c, ldim=%d, epsrel=%g", ier, LbyR, X, Y, ldim, epsrel);
    return I;
}

static void usage(char *self, FILE *stream)
{
    fprintf(stream, "Usage: %s LbyR [bc ldim epsrel]\n\n", self);

    fprintf(stream, "Compute the Casimir energy for a scalar field with perfect conductors at T=0.\n\n");

    fprintf(stream, "Options:\n");
    fprintf(stream, "  LbyR:   aspect ratio L/R\n");
    fprintf(stream, "  bc:     boundary condition on plate and sphere (DD, DN, ND or NN)\n");
    fprintf(stream, "  ldim:   dimension of vector space\n");
    fprintf(stream, "  epsrel: relative accuracy of integration\n\n");

    fprintf(stream, "The number of threads can be set using OMP_NUM_THREADS.\n");
}

int main(int argc, char *argv[])
{
    double E, LbyR, epsrel = 5e-7;
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

    E = casimir_E(LbyR, X, Y, ldim, epsrel);

    printf("\n");
    printf("# L/R, X, Y, ldim, E*(L+R)/(hbar*c)\n");
    printf("%.15g, %c, %c, %d, %.15g\n", LbyR, X, Y, ldim, E);

    return 0;
}
