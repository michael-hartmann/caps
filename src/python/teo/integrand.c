#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BUF_SIZE 134217728

static double kahan_sum(double input[], size_t N)
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

double integrand(const double tau, const double l, const double eps, const double e)
{
	double *terms = malloc(BUF_SIZE*sizeof(double));
    if(terms == NULL)
    {
        fprintf(stderr, "Cannot allocate memory.\n");
        abort();
    }

    //const double eps2 = eps*eps; /* ε² */

    const double tau2 = tau*tau;   /* τ² */
    //const double tau4 = tau2*tau2; /* τ^4 */
    //const double tau6 = tau4*tau2; /* τ^6 */

    const double alpha = (1+tau)*(1-tau); /* 1-τ² */
    //const double alpha2 = alpha*alpha;    /* (1-τ²)² */
    //const double alpha3 = alpha2*alpha;   /* (1-τ²)³ */

    const double beta2 = eps*alpha+tau2; /* β² = ε*(1-τ²)+τ² */
    const double beta  = sqrt(beta2);    /* β */
    //const double beta3 = beta2*beta;     /* β³ */

    const double T0_TE = (beta-1)/(beta+1);
    const double T0_TM = (eps-beta)/(eps+beta);
    const double T0_TE2 = T0_TE*T0_TE;
    const double T0_TM2 = T0_TM*T0_TM;

    #if 0
    const double K1_TE = -2*tau/beta;
    const double K1_TM = 2*eps*tau*alpha/(beta*(eps+tau2));
    const double W1_TE = -4*tau/beta; /* 2*K1_TE */
    const double W1_TM = 4*eps*tau*alpha/(beta*(eps+tau2)); /* 2*K1_TM */

    const double K2_TE = -eps*alpha/beta3 + 2*tau2/beta2;
    const double K2_TM = eps2*alpha2/(beta3*(eps+tau2)) -tau2*(-eps2*tau2+eps2+eps+1)/(beta2*(eps+tau2));
    const double W2_TE = (8*tau2+4*tau4+4*eps-4*eps*tau4)/beta3 + 4*alpha2*(eps+beta)*(eps+beta)/(tau2*beta2*(beta+1)*(beta+1)) - 4*alpha*(tau2+eps)/(tau2*beta2);
    const double W2_TM = -eps*alpha*(8*tau2+4*tau4+4*eps-4*eps*tau4)/((eps+tau2)*beta3) + 4*alpha2*eps2*(1+beta)*(1+beta)/(tau2*beta2*(beta+eps)*(beta+eps)) - 4*eps2*alpha3/(tau2*(tau2+eps)*beta2);
    const double Y2_TE = -tau/beta - (8*eps*tau2-3*eps-5*eps*tau4+9*tau2+5*tau4)/(12*beta2);
    const double Y2_TM = eps*tau*alpha/((eps+tau2)*beta) - (7*eps2*tau4-4*eps2*tau2-3*eps2-5*eps*tau6+13*eps*tau4-18*eps*tau2+5*tau6-3*tau4)/(12*(eps+tau2)*beta2);
    #endif

    double qs = 1; /* q^s */
    const double q = exp(-2*e*l/tau);

    double T0_TE_2s = 1;
    double T0_TM_2s = 1;

    for(int i = 0; i < BUF_SIZE; i++)
    {
        qs *= q;

        T0_TE_2s *= T0_TE2; /* T0_TE^(2s) */
        T0_TM_2s *= T0_TM2; /* T0_TM^(2s) */

        const double s = i+1;   /* s */
        const double s2 = s*s;  /* s² */
        #if 0
        const double s3 = s2*s; /* s³ */

        const double A = e*e*l*tau/3*(s3+2*s) + e/3*((tau2-2)*s2-3*tau*s+2*tau2-1) + (tau4+tau2-12)/(12*l*tau)*s + (1+tau)*alpha/(2*l*tau) - tau*alpha/(3*l*s);
        const double B = alpha/(2*l*tau);

        const double X = 2*T0_TE*T0_TM*(T0_TE_2s-T0_TM_2s)/(T0_TE2-T0_TM2) + 2*(T0_TE2*T0_TM2)*(T0_TE_2s/T0_TE2-T0_TM_2s/T0_TM2)/(T0_TE2-T0_TM2);

        const double CV = -e*tau/3*(s3+2*s)+alpha/(6*l)*s2 + tau/(2*l)*s + (1-4*tau2)/(12*l);
        const double CJ = -e*tau/6*(s3-s) + (s2-1)/(12*l);
        const double C_TE = CV*K1_TE+CJ*W1_TE;
        const double C_TM = CV*K1_TM+CJ*W1_TM;

        const double DVV = tau/(12*l)*(s3-2*s2+2*s-1);
        const double DVJ = tau/(12*l)*(s3-s);
        const double DJJ = tau/(48*l)*(s3-2*s2-s+2);
        const double DV  = tau/(6*l)*(2*s2-3*s+1);
        const double DJ  = tau/(12*l)*(s2-1);
        const double D_TE = DVV*K1_TE*K1_TE + DVJ*K1_TE*W1_TE + DJJ*W1_TE*W1_TE + (s*tau/(2*l)+DV)*K2_TE + DJ*W2_TE + s*tau/l*Y2_TE;
        const double D_TM = DVV*K1_TM*K1_TM + DVJ*K1_TM*W1_TM + DJJ*W1_TM*W1_TM + (s*tau/(2*l)+DV)*K2_TM + DJ*W2_TM + s*tau/l*Y2_TM;
        #endif

        //const double v = qs/s2*( T0_TE_2s*(A+C_TE+D_TE) + T0_TM_2s*(A+C_TM+D_TM) + X*B );
        const double v = qs/s2*( T0_TE_2s + T0_TM_2s );

        terms[i] = v;

        //printf("%d, %e\n", i, v/terms[0]);
        if(v == 0 || fabs(v/terms[0]) < 1e-12)
		{
			const double sum = kahan_sum(terms, i);
            free(terms);
            return sum;
		}
    }

    fprintf(stderr, "sum over s did not converge: tau=%g, l=%g, eps=%g, e=%g\n", tau, l, eps, e);
	abort();
    return NAN;
}
