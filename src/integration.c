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
static double log_gausskronrod[15][3] =
{
    /* nodes                 weight Gauß           weight Kronrod */
    {  0.00000000000000000, -0.87237149793931956, -1.5631167886725021 },
    { +0.40584515137739717, -0.96277966333590204, -1.6588877593062039 },
    { -0.40584515137739717, -0.96277966333590204, -1.6588877593062039 },
    { +0.74153118559939444, -1.2740184029883301,  -1.9614575682357287 },
    { -0.74153118559939444, -1.2740184029883301,  -1.9614575682357287 },
    { +0.94910791234275852, -2.0441904959417476,  -2.7631598321848663 },
    { -0.94910791234275852, -2.0441904959417476,  -2.7631598321848663 },

    { +0.20778495500789847, -INFINITY, -1.5875152786694389 },
    { -0.20778495500789847, -INFINITY, -1.5875152786694389 },
    { +0.58608723546769113, -INFINITY, -1.7778285961704769 },
    { -0.58608723546769113, -INFINITY, -1.7778285961704769 },
    { +0.86486442335976907, -INFINITY, -2.2557968329911389 },
    { -0.86486442335976907, -INFINITY, -2.2557968329911389 },
    { +0.99145537112081264, -INFINITY, -3.7750771108951247 },
    { -0.99145537112081264, -INFINITY, -3.7750771108951247 }
};



/* We want to integrate:
 * A = A_0 Λ(l1,l2,m) m² exp(-τ)*τ  ∫ dz r_p exp(-z) P_l1^m  P_l2^m  / (z²+2τz)
 * B = B_0 Λ(l1,l2,m)    exp(-τ)/τ³ ∫ dz r_p exp(-z) P'_l1^m P'_l2^m   (z²+2τz)
 * C = C_0 Λ(l1,l2,m) m  exp(-τ)/τ  ∫ dz r_p exp(-z) P_l1^m  P'_l2^m
 * D = D_0 Λ(l1,l2,m) m  exp(-τ)/τ  ∫ dz r_p exp(-z) P'_l1^m P_l2^m
 *
 * The integration is from 0 to ∞, the associated Legendre polynomials and
 * their derivatives are evaluated at the argument arg=1+x/τ.
 *
 * Here, we make the substitution z=t/(1-t), dz = dt/(1-t)². The range of
 * integration is now t=0...1.
 *
 * This function returns
 *  v[A_TE] = log(| Λ m² exp(-τ)*τ  1/(1-t)² r_TE exp(-z) P_l1^m  P_l2^m  / (z²+2τz) |)
 *  v[A_TM] = log(| Λ m² exp(-τ)*τ  1/(1-t)² r_TM exp(-z) P_l1^m  P_l2^m  / (z²+2τz) |)
 *
 *  v[B_TE] = log(| Λ    exp(-τ)/τ³ 1/(1-t)² r_TE exp(-z) P'_l1^m P'_l2^m   (z²+2τz) |)
 *  v[B_TM] = log(| Λ    exp(-τ)/τ³ 1/(1-t)² r_TM exp(-z) P'_l1^m P'_l2^m   (z²+2τz) |)
 *
 *  v[C_TE] = log(| Λ m  exp(-τ)/τ  1/(1-t)² r_TE exp(-z) P_l1^m  P'_l2^m            |)
 *  v[C_TM] = log(| Λ m  exp(-τ)/τ  1/(1-t)² r_TM exp(-z) P_l1^m  P'_l2^m            |)
 *
 *  v[D_TE] = log(| Λ m  exp(-τ)/τ  1/(1-t)² r_TE exp(-z) P'_l1^m P_l2^m             |)
 *  v[D_TM] = log(| Λ m  exp(-τ)/τ  1/(1-t)² r_TM exp(-z) P'_l1^m P_l2^m             |)
 * where
 *      z = t/(1-t)
 *      τ = 2nT
 *      Λ = Λ(l1,l2,m): prefactor
 *      r_p: Fresnel coefficient for p=TE,TM.
 *
 * The signs are stored in s[A_TE], s[A_TM], ..., s[D_TM].
 */
void casimir_integrate_integrands(integration_t *int_obj, double t, int l1, int l2, double v[8], sign_t s[8])
{
    plm_combination_t comb;

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

    double lnLambda = casimir_lnLambda(l1, l2, m);

    plm_PlmPlm(l1, l2, m, 1+z/tau, &comb);

    /* Λ exp(-τ)/τ³ 1/(1-t)² r_p exp(-z) dPl1m dPl2m (z²+2τz)*/
    double B_ = lnLambda -tau -3*log_tau -z + log_term + comb.lndPl1mdPl2m +log_dz;
    v[B_TE] = log_rTE + B_; /* B, TE */
    v[B_TM] = log_rTM + B_; /* B, TM */
    s[B_TE] = -comb.sign_dPl1mdPl2m;
    s[B_TM] = +comb.sign_dPl1mdPl2m;

    if(m > 0)
    {
        double log_m = log(m);

        /* Λ m² exp(-τ)*τ 1/(1-t)² r_p exp(-z) Pl1m Pl2m / (z²+2τz) */
        double A_ = lnLambda +2*log_m +log_tau -tau -z - log_term +comb.lnPl1mPl2m +log_dz;
        v[A_TE] = log_rTE + A_; /* A, TE */
        v[A_TM] = log_rTM + A_; /* A, TM */
        s[A_TE] = -comb.sign_Pl1mPl2m;
        s[A_TM] = +comb.sign_Pl1mPl2m;

        /* Λ m exp(-τ)/τ 1/(1-t)² r_p exp(-z) Pl1m dPl2m */
        double C_ = lnLambda +log_m - tau - log_tau -z +comb.lnPl1mdPl2m +log_dz;
        v[C_TE] = log_rTE + C_; /* C, TE */
        v[C_TM] = log_rTM + C_; /* C, TM */
        s[C_TE] = -comb.sign_Pl1mdPl2m;
        s[C_TM] = +comb.sign_Pl1mdPl2m;

        /* Λ m exp(-τ)/τ 1/(1-t)² r_p exp(-z) dPl1m Pl2m */
        double D_ = lnLambda +log_m -tau -log_tau -z +comb.lndPl1mPl2m +log_dz;
        v[D_TE] = log_rTE + D_; /* D, TE */
        v[D_TM] = log_rTM + D_; /* D, TM */
        s[D_TE] = -comb.sign_dPl1mPl2m;
        s[D_TM] = +comb.sign_dPl1mPl2m;
    }
}


static void integrate_gauss_kronrod(integration_t *int_obj, int l1, int l2, double a, double b, interval_t *interval);

/* a: left border of integration
 * b: right border of integration
 */
static void integrate_gauss_kronrod(integration_t *int_obj, int l1, int l2, double a, double b, interval_t *interval)
{
    const double dx = (b-a)/2;
    const double log_dx = log(dx);
    double points_G7[8][7];
    double points_K15[8][15];
    sign_t signs[8][15];

    TERMINATE(l1 <= 0 || l2 <= 0, "l1,l2 > 0, l1=%d, l2=%d\n", l1,l2);
    TERMINATE(b <= a, "b > a, b=%g, a=%g", b, a);

    interval->a = a;
    interval->b = b;

    for(int i = 0; i < 15; i++)
    {
        double xi = log_gausskronrod[i][0]; /* node Kronrod */
        double zi = (xi+1)*dx+a; /* corresponding node in interval [a,b] */

        double log_wiG = log_gausskronrod[i][1]; /* weight Gauss */
        double log_wiK = log_gausskronrod[i][2]; /* weight Kronrod */

        /* calculate integrands A_TE, A_TM, ..., D_TE, D_TM at node zi */
        double v[8];
        sign_t s[8];
        casimir_integrate_integrands(int_obj, zi, l1, l2, v, s);

        if(i < 7)
        {
            for(int j = 0; j < 8; j++)
                points_G7[j][i] = log_dx+log_wiG+v[j];
        }

        for(int j = 0; j < 8; j++)
        {
            points_K15[j][i] = log_dx+log_wiK+v[j];
            signs[j][i] = s[j];
        }
    }

    for(int i = 0; i < 8; i++)
    {
        sign_t sign_G7, sign_K15, dummy;

        double G7  = logadd_ms(points_G7[i],  signs[i], 7,  &sign_G7);
        double K15 = logadd_ms(points_K15[i], signs[i], 15, &sign_K15);

        /* 200|G7-K15| */
        const double log200 = 5.298317366548036; /* log(200) */
        double err = log200 + logadd_s(G7, +1, K15, -1, &dummy);

        interval->K15[i]   = K15;
        interval->signs[i] = sign_K15;
        interval->err[i]   = err;
    }

    interval->maxerr = max(interval->err, 8);
}

int casimir_integrate_init(casimir_t *casimir, integration_t *int_obj, double nT, int m)
{
    TERMINATE(m < 0, "m > 0, m=%d", m);
    TERMINATE(nT <= 0, "nT > 0, nT=%g", nT);

    int_obj->casimir = casimir;
    int_obj->m   = m;
    int_obj->nT  = nT;
    int_obj->tau = 2*nT;

    return 0;
}

int casimir_integrate_free(integration_t *self)
{
    /* NOP at the moment */
    return 0;
}


static double estimate_error(interval_t intervals[], int N, double v[8], sign_t s[8], double relerror[8], int *index);

static double estimate_error(interval_t intervals[], int N, double v[8], sign_t s[8], double relerror[8], int *index)
{
    double err2[8][N];
    double K15[8][N];
    sign_t signs[8][N];

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

    for(int j = 0; j < 8; j++)
    {
        v[j] = logadd_ms(K15[j], signs[j], N, &s[j]);
        if(v[j] == -INFINITY)
            relerror[j] = -INFINITY;
        else
            relerror[j] = (logadd_m(err2[j], N) - log(N*(N-1.)))/2 - v[j];
    }

    int jmax = 0;
    for(int j = 1; j < 8; j++)
    {
        if(relerror[j] > relerror[jmax])
            jmax = j;
    }

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
    #define MAXERROR -22

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

    /* NMIN intervals */
    for(int i = 0; i < N; i++)
    {
        double a = (double)i/N;
        double b = (i+1.)/N;

        integrate_gauss_kronrod(self, l1, l2, a, b, &intervals[i]);
    }

    while(true)
    {
        double error = estimate_error(intervals, N, v, s, relerror, &index);
        if(error < MAXERROR)
        {
            /* finished */
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

        if(N >= NMAX)
        {
            /* give up :( */
            TERMINATE(1, "Integral did not converge. l1=%d, l2=%d, m=%d, nT=%g", l1,l2,self->m,self->nT);
            return 1;
        }

        /* bisect interval with largest error */
        double a = intervals[index].a;
        double b = intervals[index].b;
        double m = (a+b)/2;

        integrate_gauss_kronrod(self, l1, l2, a, m, &intervals[index]);
        integrate_gauss_kronrod(self, l1, l2, m, b, &intervals[N]);

        N++;
    }

    return 0;
}
