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

/* nodes and weights for Gauß-Kronrod (G7,K15) */
static double gausskronrod[15][3] =
{
    /* nodes,             weight Gauß,       weight Kronrod */
    { +0.949107912342759, 0.129484966168870, 0.063092092629979 },
    { -0.949107912342759, 0.129484966168870, 0.063092092629979 },
    { +0.741531185599394, 0.279705391489277, 0.140653259715525 },
    { -0.741531185599394, 0.279705391489277, 0.140653259715525 },
    { +0.405845151377397, 0.381830050505119, 0.190350578064785 },
    { -0.405845151377397, 0.381830050505119, 0.190350578064785 },
    {  0.000000000000000, 0.417959183673469, 0.209482141084728 },

    { +0.991455371120813, 0, 0.022935322010529 },
    { -0.991455371120813, 0, 0.022935322010529 },
    { +0.864864423359769, 0, 0.104790010322250 },
    { -0.864864423359769, 0, 0.104790010322250 },
    { +0.586087235467691, 0, 0.169004726639267 },
    { -0.586087235467691, 0, 0.169004726639267 },
    { +0.207784955007898, 0, 0.204432940075298 },
    { -0.207784955007898, 0, 0.204432940075298 }
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
 *         \---------- prefactor ------------/ \------ what we calculate in this function ------/
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
        v[i] = 0;

    /* t=1 corresponds to z=∞; the integrands vanish for z→∞ */
    if(t == 1)
        return;

    double tau = int_obj->tau;
    double log_tau = log(tau);
    double z = t/(1-t);
    double log_dz = -2*log1p(-t); /* 1/(1-t)² */
    double log_term = log(pow_2(z)+2*tau*z); /* z²+2τz */

    /* calculate Fresnel coefficients */
    double rTE,rTM;
    {
        double k = sqrt(pow_2(z)+2*tau*z)/2; /* sqrt(z²+2τz)/2 */
        casimir_rp(int_obj->casimir, int_obj->nT, k, &rTE, &rTM);
    }

    /* calculate products of associated Legendre polynomials: Pl1m*Pl2m,
     * dPl1m*dPl2m, dPl1m*Pl2m, Pl1m*dPl2m */
    int m = int_obj->m;
    plm_combination_t comb;
    plm_PlmPlm(l1, l2, m, 1+z/tau, &comb);

    /* prefactor 1/τ³ 1/(1-t)² r_p exp(-z) P'_l1^m P'_l2^m   (z²+2τz) */
    double B = comb.sign_dPl1mdPl2m*exp(log_prefactor -3*log_tau -z +log_term +comb.lndPl1mdPl2m +log_dz);
    v[B_TE] = rTE*B; /* B, TE */
    v[B_TM] = rTM*B; /* B, TM */

    if(m > 0)
    {
        /* prefactor m² τ 1/(1-t)² r_p exp(-z) P_l1^m  P_l2^m  / (z²+2τz) */
        double A = comb.sign_Pl1mPl2m*pow_2(m)*exp(log_prefactor +log_tau -z -log_term +comb.lnPl1mPl2m +log_dz);
        v[A_TE] = rTE*A; /* A, TE */
        v[A_TM] = rTM*A; /* A, TM */

        /* prefactor m 1/τ  1/(1-t)² r_p exp(-z) P_l1^m  P'_l2^m */
        double C = comb.sign_dPl1mPl2m*m*exp(log_prefactor -log_tau -z +comb.lnPl1mdPl2m +log_dz);
        v[C_TE] = rTE*C; /* C, TE */
        v[C_TM] = rTM*C; /* C, TM */

        /* prefactor m 1/τ 1/(1-t)² r_p exp(-z) P'_l1^m P_l2^m */
        double D = comb.sign_dPl1mPl2m*m*exp(log_prefactor -log_tau -z +comb.lndPl1mPl2m +log_dz);
        v[D_TE] = rTE*D; /* D, TE */
        v[D_TM] = rTM*D; /* D, TM */
    }
}


static void integrate_gauss_kronrod(integration_t *int_obj, int l1, int l2, int k, int N, double log_prefactor, interval_t *interval);

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
static void integrate_gauss_kronrod(integration_t *int_obj, int l1, int l2, int k, int N, double log_prefactor, interval_t *interval)
{
    const double a = (k+0.0)/N;
    //const double b = (k+1.0)/N;
    const double dx = 1./(2*N);

    TERMINATE(l1 <= 0 || l2 <= 0, "l1,l2 > 0, l1=%d, l2=%d\n", l1,l2);

    interval->k = k;
    interval->N = N;

    log_prefactor += log(dx);

    /* calculate integrands at nodes and save results in G7 for Gauß, and K15
     * for Kronrod */
    double G7[8]  = { 0 };
    double K15[8] = { 0 };

    double v[8];
    for(int i = 0; i < 15; i++)
    {
        double xi = gausskronrod[i][0]; /* node Kronrod */
        double zi = (xi+1)*dx+a; /* transform node to interval [a,b] */

        double wiG = gausskronrod[i][1]; /* weight Gauss */
        double wiK = gausskronrod[i][2]; /* weight Kronrod */

        /* calculate integrands A_TE, A_TM, ..., D_TE, D_TM at node zi */
        casimir_integrate_integrands(int_obj, zi, l1, l2, log_prefactor, v);

        for(int j = 0; j < 8; j++)
        {
            double vj = v[j];
            G7[j]  += wiG*vj;
            K15[j] += wiK*vj;
        }
    }

    /* copy to struct */
    for(int j = 0; j < 8; j++)
    {
        interval->K15[j] = K15[j];
        interval->err[j] = fabs(K15[j]-G7[j]);
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
    double nT = n*casimir->T;

    TERMINATE(m < 0, "m >= 0, m=%d", m);
    TERMINATE(n <= 0, "n > 0, n=%d", n);

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


static double estimate_error(interval_t intervals[], int N, double v[8], int *index);

static double estimate_error(interval_t intervals[], int N, double v[8], int *index)
{
    double sum_err2[8] = { 0 }; /* squared error */
    double maxerr = 0;
    *index = 0;

    for(int j = 0; j < 8; j++)
        v[j] = 0;

    /* copy value of integral and error of every subinterval to arrays K15 and
     * err2 */
    for(int i = 0; i < N; i++)
    {
        interval_t *interval = &intervals[i];

        for(int j = 0; j < 8; j++)
        {
            double err = interval->err[j];
            if(err > maxerr)
            {
                maxerr = err;
                *index = i;
            }

            sum_err2[j] += pow_2(err);
            v[j] += interval->K15[j];
        }
    }

    /* We assume that errors in subintervals are independent, then:
     *      err = sqrt( 1/(N*(N-1)) \sum_j err_j^2 )
     */
    double relerror[8];
    for(int j = 0; j < 8; j++)
    {
        if(sum_err2[j] == 0)
            relerror[j] = 0;
        else
            relerror[j] = sum_err2[j]/(v[j]*sqrt(N*(N-1.)));
    }

    return max(relerror, 8);
}


int casimir_integrate(integration_t *self, int l1, int l2, double v[8])
{
    #define NMIN  4 /* minimum number of intervals */
    #define PRECISION 1e-20 /* XXX */

    /* Use fixed size arrays on stack instead dynamic allocation with malloc.
     * If we need more than NMAX intervals, we have a serious problem and
     * should terminate with an error. This also avoids infinity loops.
     * Using the stack also makes execution faster and code simpler.
     */
    interval_t intervals[INTEGRATE_INTERVALS_MAX];

    /* calculate prefactor */
    double log_prefactor;
    {
        casimir_t *casimir = self->casimir;
        double log_Lambda = casimir_lnLambda(l1, l2, self->m);

        int n = self->n;
        sign_t dummy1,dummy2;
        double log_al1, log_al2, log_bl1, log_bl2;
        casimir_mie_cache_get(casimir, l1, n, &log_al1, &dummy1, &log_bl1, &dummy2);
        casimir_mie_cache_get(casimir, l2, n, &log_al2, &dummy1, &log_bl2, &dummy2);

        /* prefactor = Λ(l1,l2,m) √(|a_l1|*|a_l2|) exp(-τ) */
        log_prefactor = log_Lambda + (log_al1+log_al2)/2 -self->tau;
    }

    /* we label the left and right border of every interval by k and N:
     *      a = k/N   and   b = (k+1)/N
     * For convenience, N will always be a power of 2.
     */

    /* we start with NMIN subintervals */
    int subintervals = NMIN;

    /* Calculate integral of every subinterval using (G7,K15) */
    for(int k = 0; k < subintervals; k++)
        integrate_gauss_kronrod(self, l1, l2, k, subintervals, log_prefactor, &intervals[k]);

    while(subintervals < INTEGRATE_INTERVALS_MAX)
    {
        int index;

        /* sum integrals of every subinterval and estimate error */
        double error = estimate_error(intervals, subintervals, v, &index);

        /* if error is sufficiently small, we're done */
        if(error < PRECISION)
        {
            const int m = self->m;
            /* signs of integrals A,B,C and D */
            const sign_t a0 = A0(l1,l2,m);
            const sign_t b0 = B0(l1,l2,m);
            const sign_t c0 = C0(l1,l2,m);
            const sign_t d0 = D0(l1,l2,m);

            v[A_TE] = a0*v[A_TE];
            v[A_TM] = a0*v[A_TM];

            v[B_TE] = b0*v[B_TE];
            v[B_TM] = b0*v[B_TM];

            v[C_TE] = c0*v[C_TE];
            v[C_TM] = c0*v[C_TM];

            v[D_TE] = d0*v[D_TE];
            v[D_TM] = d0*v[D_TM];

            return 0;
        }

        /* accuracy is not high enough; split interval with largest error */
        int k = intervals[index].k;
        int N = intervals[index].N;

        integrate_gauss_kronrod(self, l1, l2, 2*k,   2*N, log_prefactor, &intervals[index]);
        integrate_gauss_kronrod(self, l1, l2, 2*k+1, 2*N, log_prefactor, &intervals[subintervals]);

        /* number of subintervals has increased */
        subintervals++;
    }

    /* give up :( */
    TERMINATE(1, "Integral did not converge. l1=%d, l2=%d, m=%d, nT=%g", l1,l2,self->m,self->nT);
    return 1;
}
