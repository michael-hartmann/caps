/* Implementation of the integrand for the linear correction to PFA for
 * arbitrary materials at zero temperature.
 *
 * Refernce: Teo, Material dependence of Casimir interaction between a sphere
 * and a plate: First analytic correction beyond proximity force approximation,
 * https://doi.org/10.1103/PhysRevD.88.045019
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double integrand(const double xi_, const double x, const double eps, const double e);

#define BUF_ELEMS 16777216L

/**
 * @brief Sum elements in array
 *
 * Compute sum of the elements of the array input. The function uses Kahan
 * summation algorithm to reduce numerical error.
 *
 * The algorithm is taken from Wikipedia, see
 * https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 *
 * @param [in] input array
 * @param [in] N length of array
 * @return sum sum of array elements
 */
static double __attribute__((optimize("O0"))) kahan_sum(double input[], size_t N)
{
    double sum = 0;
    double c = 0; /* running compensation for lost low-order bits */

    for(size_t i = 0; i < N; i++)
    {
        double y = input[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}

/* see also theory/teo_prd.pdf */
double integrand(const double xi_, const double x, const double eps, const double e)
{
    double *terms = malloc(BUF_ELEMS*sizeof(double));
    if(terms == NULL)
    {
        fprintf(stderr, "Cannot allocate memory.\n");
        abort();
    }

    /* in the following, we assume that both sphere and plate have the same
     * dielectric properties */

    const double f = xi_/x;
    const double tau2 = (1+f)*(1-f); /* τ² */
    const double tau = sqrt(tau2);   /* τ */
    const double tau4 = tau2*tau2;   /* τ^4 */
    const double tau6 = tau4*tau2;   /* τ^6 */
    const double l    = tau*x/e;

    const double eps2 = eps*eps; /* ε² */

    const double alpha = f*f;           /* 1-τ² */
    const double alpha2 = alpha*alpha;  /* (1-τ²)² */
    const double alpha3 = alpha2*alpha; /* (1-τ²)³ */

    const double beta2 = f*f*(eps-1)+1; /* β² = ε*(1-τ²)+τ² */
    const double beta  = sqrt(beta2);   /* β */
    const double beta3 = beta2*beta;    /* β³ */

    const double rTE = (beta-1)/(beta+1);     /* T0_TE */
    const double rTM = (eps-beta)/(eps+beta); /* T0_TM */
    const double rTE2 = rTE*rTE; /* r_TE² */
    const double rTM2 = rTM*rTM; /* r_TM² */

    const double K1_TE = -2*tau/beta;
    const double K1_TM = 2*eps*tau*alpha/(beta*(eps+tau2));
    const double W1_TE = -4*tau/beta; /* 2*K1_TE */
    const double W1_TM = 4*eps*tau*alpha/(beta*(eps+tau2)); /* 2*K1_TM */

    const double K2_TE = -eps*alpha/beta3 + 2*tau2/beta2;
    const double K2_TM = eps2*alpha2/(beta3*(eps+tau2)) -tau2*(-eps2*tau2+eps2+eps+1)/(beta2*(eps+tau2)) + (tau2*(eps*beta+1)*(eps*beta+1))/(beta2*(beta+eps)*(beta+eps));

    const double W2_TE = (8*tau2+4*tau4+4*eps-4*eps*tau4)/beta3 + 4*alpha2*(eps+beta)*(eps+beta)/(tau2*beta2*(beta+1)*(beta+1)) - 4*alpha*(tau2+eps)/(tau2*beta2);
    const double W2_TM = -eps*alpha*(8*tau2+4*tau4+4*eps-4*eps*tau4)/((eps+tau2)*beta3) + 4*alpha2*eps2*(1+beta)*(1+beta)/(tau2*beta2*(beta+eps)*(beta+eps)) - 4*eps2*alpha3/(tau2*(tau2+eps)*beta2);
    const double Y2_TE = -tau/beta - (8*eps*tau2-3*eps-5*eps*tau4+9*tau2+5*tau4)/(12*beta2);
    const double Y2_TM = eps*tau*alpha/((eps+tau2)*beta) - (7*eps2*tau4-4*eps2*tau2-3*eps2-5*eps*tau6+13*eps*tau4-18*eps*tau2+5*tau6-3*tau4)/(12*(eps+tau2)*beta2);

    double qs = 1; /* q^s */
    const double q = exp(-2*e*l/tau);

    double rTE_2s = 1; /* r_TE^(2s) */
    double rTM_2s = 1; /* r_TM^(2s) */

    for(size_t i = 0; i < BUF_ELEMS; i++)
    {
        qs *= q;

        rTE_2s *= rTE2; /* rTE^(2s) */
        rTM_2s *= rTM2; /* rTM^(2s) */

        const double s = i+1;   /* s */
        const double s2 = s*s;  /* s² */
        const double s3 = s2*s; /* s³ */

        const double A = e*e*l*tau/3*(s3+2*s) + e/3*((tau2-2)*s2-3*tau*s+2*tau2-1) + (tau4+tau2-12)/(12*l*tau)*s + (1+tau)*alpha/(2*l*tau) - tau*alpha/(3*l*s);

        const double B = alpha/(l*tau);
        const double X_num   = (rTE*rTM+rTM2)*rTE_2s - (rTE*rTM+rTE2)*rTM_2s;
        const double X_denom = (rTE+rTM)*(rTE-rTM);
        const double X = X_num/X_denom;

        const double CV = -e*tau/3*(s3+2*s)+alpha/(6*l)*s2 + tau/(2*l)*s + (1-4*tau2)/(12*l);
        const double CJ = -e*tau/6*(s3-s) + (s2-1)/(12*l);
        const double C_TE = CV*K1_TE+CJ*W1_TE;
        const double C_TM = CV*K1_TM+CJ*W1_TM;

        const double DVV = tau/(12*l)*(s3-2*s2+2*s-1);
        const double DVJ = tau/(12*l)*(s3-s);
        const double DJJ = tau/(48*l)*(s3-2*s2-s+2);
        const double DV  = tau/(6*l) *(2*s2-3*s+1);
        const double DJ  = tau/(12*l)*(s2-1);
        const double D_TE = DVV*K1_TE*K1_TE + DVJ*K1_TE*W1_TE + DJJ*W1_TE*W1_TE + (s*tau/(2*l)+DV)*K2_TE + DJ*W2_TE + s*tau/l*Y2_TE;
        const double D_TM = DVV*K1_TM*K1_TM + DVJ*K1_TM*W1_TM + DJJ*W1_TM*W1_TM + (s*tau/(2*l)+DV)*K2_TM + DJ*W2_TM + s*tau/l*Y2_TM;

        const double v = qs/s2*( rTE_2s*(A+C_TE+D_TE) + rTM_2s*(A+C_TM+D_TM) + X*B );

        /* only PFA result: const double v = qs/s2*( rTE_2s + rTM_2s ); */

        if(isnan(v))
        {
            printf("v=nan, i=%zu, xi_=%.16g, x=%.16g, tau=%.16g, l=%.16g\n", i, xi_, x, tau, l);
            printf("%g %g %g %g %g %g %g %g %g %g %g\n", qs, s2, rTE_2s, rTM_2s, A, C_TE, D_TE, C_TM, D_TM, X, B);
            abort();
        }

        terms[i] = v;

        if(v == 0 || fabs(v/terms[0]) < 1e-20)
        {
            const double sum = kahan_sum(terms, i+1);
            free(terms);
            return sum;
        }
    }

    fprintf(stderr, "sum over s did not converge: tau=%g, l=%g, eps=%g, e=%g, terms[0]=%g\n", tau, l, eps, e, terms[0]);
    abort();
    return NAN;
}
