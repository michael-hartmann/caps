/**
 * @file   integration.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   September, 2016
 * @brief  Perform integration for Drude planes
 */

#include <assert.h>
#include <math.h>

#include "utils.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"
#include "gausslaguerre.h"

/* nodes and weights vor Gauss-Kronrod */
static double gausskronrod[15][3] =
{
    /* node               weight Gauss       weight Kronrod */
    { +0.949107912342759, 0.129484966168870, 0.063092092629979 },
    { -0.949107912342759, 0.129484966168870, 0.063092092629979 },
    { +0.741531185599394, 0.279705391489277, 0.140653259715525 },
    { -0.741531185599394, 0.279705391489277, 0.140653259715525 },
    { +0.405845151377397, 0.381830050505119, 0.190350578064785 },
    { -0.405845151377397, 0.381830050505119, 0.190350578064785 },
    {  0.000000000000000, 0.417959183673469, 0.209482141084728 },

    { +0.991455371120813, 0.000000000000000, 0.022935322010529 },
    { -0.991455371120813, 0.000000000000000, 0.022935322010529 },
    { +0.864864423359769, 0.000000000000000, 0.104790010322250 },
    { -0.864864423359769, 0.000000000000000, 0.104790010322250 },
    { +0.586087235467691, 0.000000000000000, 0.169004726639267 },
    { -0.586087235467691, 0.000000000000000, 0.169004726639267 },
    { +0.207784955007898, 0.000000000000000, 0.204432940075298 },
    { -0.207784955007898, 0.000000000000000, 0.204432940075298 }
};

/* We want to integrate:
 * A = A_0 m² τ exp(-τ) ∫ dz r_p exp(-z) P_l1  P_l2  / (z²+2τz)
 * B = B_0 exp(-τ)/τ³   ∫ dz r_p exp(-z) P'_l1 P'_l2   (z²+2τz)
 * C = C_0 im exp(-τ)/τ ∫ dz r_p exp(-z) P_l1  P'_l2
 * D = D_0 im exp(-τ)/τ ∫ dz r_p exp(-z) P'_l1 P_l2
 *
 * The integration is from 0 to ∞, the associated Legendre polynomials and
 * their derivatives are evaluated at the argument arg=1+x/τ.
 *
 * Here, we make the substitution z=t/(1-t), dz = dt/(1-t)² and t=0...1.
 *
 * The prefactors are given by:
 * A_0 = (-1)^{\ell_2+m}
 * B_0 = (-1)^{\ell_2+m+1}
 * C_0 = (-1)^{\ell_2+m}
 * D_0 = (-1)^{\ell_1+\ell_2} C^m_{\ell_2\ell_1,p}
 */
void casimir_integrate_integrands(integration_t *int_obj, double t, int l1, int l2, double v[8], sign_t s[8])
{
    plm_combination_t comb;

    /* 0 <= t <= 1 */
    TERMINATE(t < 0, "invalid value for t: t=%g", t);
    TERMINATE(t > 1, "invalid value for t: t=%g", t);

    /* l1 > 0, l2 > 0 */
    TERMINATE(l1 <= 0, "l1 > 0, l1=%d\n", l1);
    TERMINATE(l2 <= 0, "l1 > 0, l1=%d\n", l2);

    /* t=1 corresponds to z=∞; the integrands vanish for z→∞ */
    if(t == 1)
    {
        /* return +0 */
        for(int i = 0; i < 8; i++)
        {
            v[i] = -INFINITY;
            s[i] = +1;
        }
        return;
    }

    int m = int_obj->m;
    double tau = int_obj->tau;
    double log_tau = log(tau);
    double z = t/(1-t);
    double log_dz = -2*log1p(-t); /* 1/(1-t)² */
    double log_term = log(pow_2(z)+2*tau*z); /* z²+2τz */
    double lnLambda = casimir_lnLambda(l1, l2, m);

    /* calculate Fresnel coefficients */
    double log_rTE, log_rTM;
    {
        double rTE,rTM;
        double k = sqrt(pow_2(z)+2*tau*z)/2; /* sqrt(z²+2τz)/2 */
        casimir_rp(int_obj->casimir, int_obj->nT, k, &rTE, &rTM);
        log_rTE = log(-rTE);
        log_rTM = log(+rTM);
    }

    plm_PlmPlm(l1, l2, m, 1+z/tau, &comb);

    /* exp(-τ)/τ³ r_p exp(-z) (z²+2τz) P'_l1*P'_l2 */
    double B_ = lnLambda -tau -3*log_tau -z + log_term + comb.lndPl1mdPl2m +log_dz;
    v[B_TE] = log_rTE + B_; /* B, TE */
    v[B_TM] = log_rTM + B_; /* B, TM */
    s[B_TE] = -comb.sign_dPl1mdPl2m;
    s[B_TM] = +comb.sign_dPl1mdPl2m;

    if(m == 0)
    {
        v[A_TE] = v[A_TM] = -INFINITY;
        s[A_TE] = s[A_TM] = 1;

        v[C_TE] = v[C_TM] = -INFINITY;
        s[C_TE] = s[C_TM] = 1;

        v[D_TE] = v[D_TM] = -INFINITY;
        s[D_TE] = s[D_TM] = 1;
    }
    else
    {
        double log_m = log(m);

        /* m² τ exp(-τ) ∫ dz r_p exp(-z) P_l1*P_l2 / (z²+2τz) */
        double A_ = lnLambda +2*log_m +log_tau -tau -z - log_term +comb.lnPl1mPl2m +log_dz;
        v[A_TE] = log_rTE + A_; /* A, TE */
        v[A_TM] = log_rTM + A_; /* A, TM */
        s[A_TE] = -comb.sign_Pl1mPl2m;
        s[A_TM] = +comb.sign_Pl1mPl2m;

        /* m exp(-τ)/τ ∫ dz r_p exp(-z) P_l1*P'_l2 */
        double C_ = lnLambda +log_m - tau - log_tau -z +comb.lnPl1mdPl2m +log_dz;
        v[C_TE] = log_rTE + C_; /* C, TE */
        v[C_TM] = log_rTM + C_; /* C, TM */
        s[C_TE] = -comb.sign_Pl1mdPl2m;
        s[C_TM] = +comb.sign_Pl1mdPl2m;

        /* m exp(-τ)/τ ∫ dz r_p exp(-z) P'_l1*P_l2 */
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

    interval->a = a;
    interval->b = b;

    assert(b > a);

    for(int i = 0; i < 15; i++)
    {
        double v[8];
        sign_t s[8];
        const double xi  = gausskronrod[i][0]; /* node Kronrod */
        const double wiG = gausskronrod[i][1]; /* weight Gauss */
        const double wiK = gausskronrod[i][2]; /* weight Kronrod */
        const double zi  = (xi+1)*dx+a; /* corresponding node in interval [a,b] */

        /* calculate integrands A_TE, A_TM, ..., D_TE, D_TM at node zi */
        casimir_integrate_integrands(int_obj, zi, l1, l2, v, s);

        if(wiG > 0)
        {
            for(int j = 0; j < 8; j++)
                points_G7[j][i] = log_dx+log(wiG)+v[j];
        }

        for(int j = 0; j < 8; j++)
        {
            points_K15[j][i] = log_dx+log(wiK)+v[j];
            signs[j][i] = s[j];
        }
    }

    for(int i = 0; i < 8; i++)
    {
        double G7, K15, err;
        sign_t sign_G7, sign_K15;

        G7  = logadd_ms(points_G7[i],  signs[i], 7,  &sign_G7);
        K15 = logadd_ms(points_K15[i], signs[i], 15, &sign_K15);

        err = exp(K15/2) * pow(200*fabs(1-exp(G7-K15)),1.5);

        interval->K15[i] = K15;
        interval->signs[i] = sign_K15;
        interval->err[i] = err;
    }

    double maxerr = interval->err[0];
    for(int i = 1; i < 8; i++)
        maxerr = MAX(maxerr, interval->err[i]);
    interval->maxerr = maxerr;
}

int casimir_integrate_init(casimir_t *casimir, integration_t *int_obj, double nT, int m)
{
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

int casimir_integrate(integration_t *self, int l1, int l2, casimir_integrals_t *cint)
{
    int m = self->m;
    const int N = 50; /* XXX */
    interval_t intervals[N];
    double v[8][N];
    sign_t s[8][N];

    for(int i = 0; i < N; i++)
    {
        double a = (double)i/N;
        double b = (i+1.)/N;

        integrate_gauss_kronrod(self, l1, l2, a, b, &intervals[i]);
        //printf("%g %g %g\n", a,b, intervals[i].signs[B_TM]*exp(intervals[i].K15[B_TM]));

        for(int j = 0; j < 8; j++)
        {
            v[j][i] = intervals[i].K15[j];
            s[j][i] = intervals[i].signs[j];
        }
    }

    cint->lnA_TE = logadd_ms(v[A_TE], s[A_TE], N, &cint->signA_TE);
    cint->lnA_TM = logadd_ms(v[A_TM], s[A_TM], N, &cint->signA_TM);

    cint->lnB_TE = logadd_ms(v[B_TE], s[B_TE], N, &cint->signB_TE);
    cint->lnB_TM = logadd_ms(v[B_TM], s[B_TM], N, &cint->signB_TM);

    cint->lnC_TE = logadd_ms(v[C_TE], s[C_TE], N, &cint->signC_TE);
    cint->lnC_TM = logadd_ms(v[C_TM], s[C_TM], N, &cint->signC_TM);

    cint->lnD_TE = logadd_ms(v[D_TE], s[D_TE], N, &cint->signD_TE);
    cint->lnD_TM = logadd_ms(v[D_TM], s[D_TM], N, &cint->signD_TM);

    cint->signA_TE *= -A0(l1,l2,m);
    cint->signA_TM *= -A0(l1,l2,m);
    cint->signB_TE *= -B0(l1,l2,m);
    cint->signB_TM *= -B0(l1,l2,m);
    cint->signC_TE *= -C0(l1,l2,m);
    cint->signC_TM *= -C0(l1,l2,m);
    cint->signD_TE *= -D0(l1,l2,m);
    cint->signD_TM *= -D0(l1,l2,m);

    return 0;
}
