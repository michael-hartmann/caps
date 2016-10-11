/**
 * @file   integration.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   October, 2016
 * @brief  Perform integration for Drude planes
 */

#include <math.h>

#include "utils.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"

/* nodes and weights for Gauß-Kronrod (G7,K15); make sure that the first 7
 * points are the Gauß nodes! */
static double gausskronrod[15][3] = {
    /* nodes,             weight Gauß,       weight Kronrod */
    { +0.949107912342759, 0.129484966168870, 0.063092092629979},
    { -0.949107912342759, 0.129484966168870, 0.063092092629979},
    { +0.741531185599394, 0.279705391489277, 0.140653259715525},
    { -0.741531185599394, 0.279705391489277, 0.140653259715525},
    { +0.405845151377397, 0.381830050505119, 0.190350578064785},
    { -0.405845151377397, 0.381830050505119, 0.190350578064785},
    {  0.000000000000000, 0.417959183673469, 0.209482141084728},

    { +0.991455371120813, 0, 0.022935322010529},
    { -0.991455371120813, 0, 0.022935322010529},
    { +0.864864423359769, 0, 0.104790010322250},
    { -0.864864423359769, 0, 0.104790010322250},
    { +0.586087235467691, 0, 0.169004726639267},
    { -0.586087235467691, 0, 0.169004726639267},
    { +0.207784955007898, 0, 0.204432940075298},
    { -0.207784955007898, 0, 0.204432940075298}
};


/**
 * @brief Evaluate integrands A,B,C,D
 *
 * Evaluate the integrands at t in [0,1] for (l1,l2,m) and save integrands in
 * array v. v contains values for A_TE, A_TM, ..., D_TE, D_TM.
 *
 * We want to integrate:
 * A = A_0 Λ(l1,l2,m) √(|a_l1|*|a_l2|) exp(-τ) m²  τ  ∫ dz r_p exp(-z) P_l1^m  P_l2^m  / (z²+2τz)
 * B = B_0 Λ(l1,l2,m) √(|a_l1|*|a_l2|) exp(-τ)   1/τ³ ∫ dz r_p exp(-z) P'_l1^m P'_l2^m   (z²+2τz)
 * C = C_0 Λ(l1,l2,m) √(|a_l1|*|a_l2|) exp(-τ) m 1/τ  ∫ dz r_p exp(-z) P_l1^m  P'_l2^m
 * D = D_0 Λ(l1,l2,m) √(|a_l1|*|a_l2|) exp(-τ) m 1/τ  ∫ dz r_p exp(-z) P'_l1^m P_l2^m
 *         \---------- prefactor ------------/
 *
 * The integration is from 0 to ∞, the associated Legendre polynomials and
 * their derivatives are evaluated at the argument arg=1+x/τ.
 *
 * Here, we make the substitution z=t/(1-t), dz = dt/(1-t)². The range of
 * integration is now t=0...1.
 *
 * This function returns
 *  v[A_TE] = prefactor m²  τ  1/(1-t)² r_TE exp(-z) P_l1^m  P_l2^m  / (z²+2τz)
 *  v[A_TM] = prefactor m²  τ  1/(1-t)² r_TM exp(-z) P_l1^m  P_l2^m  / (z²+2τz)
 *
 *  v[B_TE] = prefactor   1/τ³ 1/(1-t)² r_TE exp(-z) P'_l1^m P'_l2^m   (z²+2τz)
 *  v[B_TM] = prefactor   1/τ³ 1/(1-t)² r_TM exp(-z) P'_l1^m P'_l2^m   (z²+2τz)
 *
 *  v[C_TE] = prefactor m 1/τ  1/(1-t)² r_TE exp(-z) P_l1^m  P'_l2^m
 *  v[C_TM] = prefactor m 1/τ  1/(1-t)² r_TM exp(-z) P_l1^m  P'_l2^m
 *
 *  v[D_TE] = prefactor m 1/τ  1/(1-t)² r_TE exp(-z) P'_l1^m P_l2^m
 *  v[D_TM] = prefactor m 1/τ  1/(1-t)² r_TM exp(-z) P'_l1^m P_l2^m
 * where
 *      z = t/(1-t)
 *      τ = 2nT
 *      r_p: Fresnel coefficient for p=TE,TM.
 * The prefactor is supposed to be given as
 *      prefactor = Λ(l1,l2,m) √(|a_l1|*|a_l2|) exp(-τ).
 *
 * @param [in]  int_obj integration object, see \ref casimir_integrate_init
 * @param [in]  t position of integrand, t in [0,1]
 * @param [in]  l1 parameter l1
 * @param [in]  l2 parameter l2
 * @param [out] v output array v of log(integrand) for A_TE,A_TM,...,D_TE,D_TM
 */
void casimir_integrate_integrands(integration_t *int_obj, double t, int l1, int l2, double log_prefactor, double v[8])
{
    /* 0 <= t <= 1 */
    TERMINATE(t < 0 || t > 1, "0 <= t <= 1, invalid value for t: t=%g", t);

    /* l1 > 0, l2 > 0 */
    TERMINATE(l1 <= 0 || l2 <= 0, "l1,l2 > 0, l1=%d, l2=%d\n", l1,l2);

    /* initialize to +0 */
    for(int i = 0; i < 8; i++)
    {
        v[i] = -INFINITY;
        s[i] = +1;
    }

    /* t=1 corresponds to z=∞; the integrands vanish for z→∞ */
    if(t == 1)
        return;

    int m = int_obj->m;
    double tau = int_obj->tau;
    double log_tau = log(tau);
    double z = t/(1-t);
    double log_dz = -2*log1p(-t); /* 1/(1-t)² */
    double log_term = log(pow_2(z)+2*tau*z); /* z²+2τz */

    /* calculate Fresnel coefficients */
    double log_rTE, log_rTM;
    {
        double rTE,rTM;
        double k = sqrt(pow_2(z)+2*tau*z)/2; /* sqrt(z²+2τz)/2 */
        casimir_rp(int_obj->casimir, int_obj->nT, k, &rTE, &rTM);
        log_rTE = log(-rTE);
        log_rTM = log(+rTM);
    }

    /* calculate products of associated Legendre polynomials: Pl1m*Pl2m,
     * dPl1m*dPl2m, dPl1m*Pl2m, Pl1m*dPl2m */
    plm_combination_t comb;
    plm_PlmPlm(l1, l2, m, 1+z/tau, &comb);

    /* prefactor 1/τ³ 1/(1-t)² r_p exp(-z) P'_l1^m P'_l2^m   (z²+2τz) */
    double log_B = log_prefactor -3*log_tau -z +log_term +comb.lndPl1mdPl2m +log_dz;
    v[B_TE] = -comb.sign_dPl1mdPl2m*exp(log_rTE + log_B); /* B, TE */
    v[B_TM] = +comb.sign_dPl1mdPl2m*exp(log_rTM + log_B); /* B, TM */

    if(m > 0)
    {
        int m2 = pow_2(m); /* m² */

        /* prefactor m² τ 1/(1-t)² r_p exp(-z) P_l1^m  P_l2^m  / (z²+2τz) */
        double log_A = log_prefactor +log_tau -z -log_term +comb.lnPl1mPl2m +log_dz;
        v[A_TE] = -comb.sign_Pl1mPl2m*m2*exp(log_rTE + log_A); /* A, TE */
        v[A_TM] = +comb.sign_Pl1mPl2m*m2*exp(log_rTM + log_A); /* A, TM */

        /* prefactor m 1/τ  1/(1-t)² r_p exp(-z) P_l1^m  P'_l2^m */
        double C_ = log_prefactor -log_tau -z +comb.lnPl1mdPl2m +log_dz;
        v[C_TE] = -comb.sign_dPl1mPl2m*m*exp(log_rTE + log_C); /* C, TE */
        v[C_TM] = +comb.sign_dPl1mPl2m*m*exp(log_rTM + log_C); /* C, TM */

        /* prefactor m 1/τ 1/(1-t)² r_p exp(-z) P'_l1^m P_l2^m */
        double D_ = log_prefactor -log_tau -z +comb.lndPl1mPl2m +log_dz;
        v[D_TE] = -comb.sign_dPl1mPl2m*m*exp(log_rTE + log_D); /* D, TE */
        v[D_TM] = +comb.sign_dPl1mPl2m*m*exp(log_rTM + log_D); /* D, TM */
    }
}


static void integrate_gauss_kronrod(integration_t *int_obj, int l1, int l2, double a, double b, interval_t *interval);

/**
 * @brief Integrate integrands from a to b
 *
 * @param [in] int_obj integration object, see \ref casimir_integrate_init
 * @param [in] l1 parameter l1
 * @param [in] l2 parameter l2
 * @param [in] a left border of integration
 * @param [in] b right border of integration, 0 <= a < b <= 1
 * @param [out] interval result of integration
 */
static void integrate_gauss_kronrod(integration_t *int_obj, int l1, int l2, double a, double b, interval_t *interval)
{
    const double dx = (b-a)/2;
    double points_G7[8][7];
    double points_K15[8][15];
    sign_t signs[8][15];

    TERMINATE(l1 <= 0 || l2 <= 0, "l1,l2 > 0, l1=%d, l2=%d\n", l1,l2);
    TERMINATE(a < 0 || a >= b || b > 1, "0 <= a < b <= 1, b=%g, a=%g", b, a);

    interval->a = a;
    interval->b = b;

    /* prefactor */
    double log_prefactor;
    {
        /* prefactor = Λ(l1,l2,m) √(|a_l1|*|a_l2|) exp(-τ) dx */
        casimir_t *casimir = int_obj->casimir;
        int n = int_obj->n;
        double log_dx = log(dx);
        double log_Lambda = casimir_lnLambda(l1, l2, int_obj->m);
        double tau = log(int_obj->tau);

        sign_t dummy1,dummy2;
        double log_al1, log_al2, log_bl1, log_bl2;

        casimir_mie_cache_get(casimir, l1, n, &log_al1, &dummy1, &log_bl1, &dummy2);
        casimir_mie_cache_get(casimir, l2, n, &log_al2, &dummy1, &log_bl2, &dummy2);

        log_prefactor = log_Lambda + (log_al1+log_al2)/2 -tau +log_dx;
    }

    /* calculate integrands at nodes and save results in G7 for Gauß, and K15
     * for Kronrod */
    double G7[8] = 0; double K15[8] = 0;
    for(int i = 0; i < 15; i++)
    {
        double xi = gausskronrod[i][0]; /* node Kronrod */
        double zi = (xi+1)*dx+a; /* transform node to interval [a,b] */

        double wiG = gausskronrod[i][1]; /* weight Gauss */
        double wiK = gausskronrod[i][2]; /* weight Kronrod */

        /* calculate integrands A_TE, A_TM, ..., D_TE, D_TM at node zi */
        double v[8];
        casimir_integrate_integrands(int_obj, zi, l1, l2, log_prefactor, v);

        if(i < 7)
        {
            for(int j = 0; j < 8; j++)
                G7[i] += wiG*v[j];
        }

        for(int j = 0; j < 8; j++)
            K15[i] = wiK*v[j];
    }

    for(int i = 0; i < 8; i++)
    {
        double G7_i = G7[i];
        double K15_i = K15[i];

        interval->K15[i] = K15_i
        interval->err[i] = fabs(K15_i-G7_i);
    }
}


/**
 * @brief Initialize integration
 *
 * Initialize integration for xi=nT and m.
 *
 * @param [in]  casimir Casimir object
 * @param [out] int_obj integration object
 * @param [in] nT Matsubara frequency
 * @param [in] m magnetic number
 * @retval 0
 */
int casimir_integrate_init(casimir_t *casimir, integration_t *int_obj, int n, int m)
{
    double nT = n*casimir-T;

    TERMINATE(m < 0, "m > 0, m=%d", m);
    TERMINATE(nT <= 0, "nT > 0, nT=%g", nT);

    int_obj->casimir = casimir;
    int_obj->m   = m;
    int_obj->n   = n;
    int_obj->nT  = nT;
    int_obj->tau = 2*nT;

    return 0;
}


/**
 * @brief Free integration object
 *
 * At the moment this function does nothing.
 *
 * @param [in,out] self integration object
 * @retval 0
 */
int casimir_integrate_free(__attribute__((unused)) integration_t *self)
{
    /* NOP at the moment */
    return 0;
}


static double estimate_error(interval_t intervals[], int N, double v[8], double relerror[8], int *index);

static double estimate_error(interval_t intervals[], int N, double v[8], double relerror[8], int *index)
{
    double err2[8][N];  /* log of squared error */
    double K15[8][N];   /* log of integral according to Kronrod */
    sign_t signs[8][N]; /* signs corresponding to K15 */

    /* copy value of integral, signs and error of every subinterval to arrays
     * K15, signs, err2 */
    for(int i = 0; i < N; i++)
    {
        interval_t *interval = &intervals[i];

        for(int j = 0; j < 8; j++)
        {
            err2[j][i]  = 2*interval->err[j];
            K15[j][i]   = interval->K15[j];
            signs[j][i] = interval->signs[j];
        }
    }

    /* We assume that errors in subintervals are independent, then:
     *      err = sqrt( 1/(N*(N-1)) \sum_j err_j^2 )
     */
    for(int j = 0; j < 8; j++)
    {
        v[j] = logadd_ms(K15[j], signs[j], N, &s[j]);
        if(v[j] == -INFINITY)
            relerror[j] = -INFINITY;
        else
            relerror[j] = (logadd_m(err2[j], N) - log(N*(N-1.)))/2 - v[j];
    }

    /* find integral with largest error; i.e. A_TE, A_TM, ..., C_TE or C_TM */
    int jmax = 0;
    for(int j = 1; j < 8; j++)
    {
        if(relerror[j] > relerror[jmax])
            jmax = j;
    }

    /* find index of subinterval with largest error */
    int imax = 0;
    for(int i = 1; i < N; i++)
    {
        if(err2[jmax][i] > err2[jmax][imax])
            imax = i;
    }

    *index = imax;

    return max(relerror, 8);
}


int casimir_integrate(integration_t *self, int l1, int l2, casimir_integrals_t *cint)
{
    #define NMIN  10 /* minimum number of intervals */
    #define NMAX 150 /* maximum number of intervals */
    #define MAXERROR -20 /* XXX */

    int N = NMIN;
    double v[8];
    sign_t s[8];
    double relerror[8];
    int index;

    /* Use fixed size arrays on stack instead dynamic allocation with malloc.
     * If we need more than NMAX intervals, we have a serious problem and
     * should terminate with an error. This also avoids possible infinity
     * loops.
     * Using the stack also makes execution faster and code simpler.
     */
    interval_t intervals[NMAX];

    /* Split integral in NMIN subintervals and calculate integral in every
     * subinterval using (G7,K15) */
    for(int i = 0; i < N; i++)
    {
        double a = (double)i/N;
        double b = (i+1.)/N;

        integrate_gauss_kronrod(self, l1, l2, a, b, &intervals[i]);
    }

    while(N < NMAX)
    {
        /* sum integrals of every subinterval and estimate error */
        double error = estimate_error(intervals, N, v, s, relerror, &index);

        /* if error is sufficiently small, we're done */
        if(error < MAXERROR)
        {
            /* copy results to cint */
            cint->lnA_TE = v[A_TE];
            cint->lnA_TM = v[A_TM];

            cint->lnB_TE = v[B_TE];
            cint->lnB_TM = v[B_TM];

            cint->lnC_TE = v[C_TE];
            cint->lnC_TM = v[C_TM];

            cint->lnD_TE = v[D_TE];
            cint->lnD_TM = v[D_TM];

            int m = self->m;
            cint->signA_TE = A0(l1,l2,m) * s[A_TE];
            cint->signA_TM = A0(l1,l2,m) * s[A_TM];
            cint->signB_TE = B0(l1,l2,m) * s[B_TE];
            cint->signB_TM = B0(l1,l2,m) * s[B_TM];
            cint->signC_TE = C0(l1,l2,m) * s[C_TE];
            cint->signC_TM = C0(l1,l2,m) * s[C_TM];
            cint->signD_TE = D0(l1,l2,m) * s[D_TE];
            cint->signD_TM = D0(l1,l2,m) * s[D_TM];

            return 0;
        }

        /* accuracy is not high enough; split the interval with largest error
         * in half */
        double a = intervals[index].a; /* left */
        double b = intervals[index].b; /* right */
        double m = (a+b)/2; /* middle */

        integrate_gauss_kronrod(self, l1, l2, a, m, &intervals[index]);
        integrate_gauss_kronrod(self, l1, l2, m, b, &intervals[N]);

        /* number of subintervals has increased */
        N++;
    }

    /* give up :( */
    TERMINATE(1, "Integral did not converge. l1=%d, l2=%d, m=%d, nT=%g", l1,l2,self->m,self->nT);
    return 1;
}
