/**
 * @file   integration.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   September, 2016
 * @brief  Perform integration for Drude planes
 */

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

/*
 * A = A_0 * \int_0^\infty dz r_p e^{-(\tau+1)z} P_{\ell_1}^m(1+z) P_{\ell_2}^m(1+z) \frac{1}{z^2+2z}
 * B = B_0 * \int_0^\infty dz r_p e^{-(\tau+1)z} {P_{\ell_1}^m}^\prime(1+z) {P_{\ell_2}^m(1+z)}^\prime (z^2+2z)
 * C = C_0 * \int_0^\infty dz r_p e^{-(\tau+1)z} P_{\ell_1}^m(1+z) {P_{\ell_2}^m(1+z)}^\prime
 * D = D_0 * \int_0^\infty dz r_p e^{-(\tau+1)z} {P_{\ell_1}^m(1+z)}^\prime P_{\ell_2}^m(1+z)
 *
 * A_0 = (-1)^{\ell_2+m} m^2
 * B_0 = (-1)^{\ell_2+m+1}
 * C_0 = (-1)^{\ell_2+m} m
 * D_0 = (-1)^{\ell_1+\ell_2} C^m_{\ell_2\ell_1,p}
 */
static void calculate_integrands(integration_t *int_obj, double z, int l1, int l2, double v[8], sign_t s[8]);

static void calculate_integrands(integration_t *int_obj, double z, int l1, int l2, double v[8], sign_t s[8])
{
    plm_combination_t comb;
    const int m = int_obj->m;
    const double nT = int_obj->nT;
    const double tau = 2*nT;
    double log_z22z = log(pow_2(z)+2*z);

    /* Fresnel coefficients */
    double log_rTE, log_rTM;
    {
        double rTE,rTM;
        double k = tau/2*sqrt(pow_2(z)+2*z);
        casimir_rp(int_obj->casimir, nT, k, &rTE, &rTM);
        log_rTE = log(-rTE);
        log_rTM = log(+rTM);
    }

    plm_PlmPlm(l1, l2, m, 1+z, &comb);

    const double A = -(tau+1)*z -log_z22z +comb.lnPl1mPl2m;
    v[0] = log_rTE + A; /* A, TE */
    v[1] = log_rTM + A; /* A, TM */
    s[0] = -comb.sign_Pl1mPl2m;
    s[1] = +comb.sign_Pl1mPl2m;

    const double B = -(tau+1)*z +log_z22z +comb.lndPl1mdPl2m;
    v[2] = log_rTE + B; /* B, TE */
    v[3] = log_rTM + B; /* B, TM */
    s[2] = -comb.sign_dPl1mdPl2m;
    s[3] = +comb.sign_dPl1mdPl2m;

    const double C = -(tau+1)*z +comb.lnPl1mdPl2m;
    v[4] = log_rTE + C; /* C, TE */
    v[5] = log_rTM + C; /* C, TM */
    s[4] = -comb.sign_Pl1mdPl2m;
    s[5] = +comb.sign_Pl1mdPl2m;

    const double D = -(tau+1)*z +comb.lndPl1mPl2m;
    v[6] = log_rTE + D; /* D, TE */
    v[7] = log_rTM + D; /* D, TM */
    s[6] = -comb.sign_dPl1mPl2m;
    s[7] = +comb.sign_dPl1mPl2m;
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

    for(int i = 0; i < 15; i++)
    {
        double v[8];
        sign_t s[8];
        const double xi  = gausskronrod[i][0]; /* node Kronrod */
        const double wiG = gausskronrod[i][1]; /* weight Gauss */
        const double wiK = gausskronrod[i][2]; /* weight Kronrod */
        const double zi  = (xi+1)*dx+a; /* corresponding node in interval [a,b] */

        /* calculate integrands A_TE, A_TM, ..., D_TE, D_TM at node zi */
        calculate_integrands(int_obj, zi, l1, l2, v, s);

        if(wiG > 0)
        {
            for(int j = 0; j < 8; j++)
                points_G7[j][i] = log_dx+wiG+v[j];
        }

        for(int j = 0; j < 8; j++)
        {
            points_K15[j][i] = log_dx+wiK+v[j];
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
        interval->err[i] = err;
    }

    double maxerr = interval->err[0];
    for(int i = 1; i < 8; i++)
        maxerr = MAX(maxerr, interval->err[i]);
    interval->maxerr = maxerr;
}

int casimir_integrate_init(casimir_t *self, integration_t *int_obj, double nT, int m)
{
    int_obj->casimir = self;
    int_obj->m    = m;
    int_obj->nT   = nT;
    int_obj->lmax = self->lmax;

    return 0;
}

int casimir_integrate_free(integration_t *self)
{
    /* NOP at the moment */
    return 0;
}

int casimir_integrate(integration_t *int_drude, int l1, int l2, casimir_integrals_t *cint)
{
    #if 0
    const int N = 50; /* intervals */

    for(int i = 0; i < N; i++)
    {
        G7


        double a = ((double)i+0)/N; /* left */
        double b = ((double)i+1)/N; /* right */

        for(int j = 0; j < 15; j++)
        {
            double xi  = gausskronrod[i][0];
            double wiG = gausskronrod[i][1];
            double wiK = gausskronrod[i][2];

            const double zi  = (xi+1)*dx+a;
            const double fzi = f(zi, args);

            integral_G7  += wiG*fzi;
            integral_K15 += wiK*fzi;
        }

    }
    #endif

    /* XXX NOT IMPLEMENTED YET XXX */
    cint->lnA_TE = 0;
    cint->lnA_TM = 0;
    cint->lnB_TE = 0;
    cint->lnB_TM = 0;
    cint->lnC_TE = 0;
    cint->lnC_TM = 0;
    cint->lnD_TE = 0;
    cint->lnD_TM = 0;
    cint->signA_TE = 1;
    cint->signA_TM = 1;
    cint->signB_TE = 1;
    cint->signB_TM = 1;
    cint->signC_TE = 1;
    cint->signC_TM = 1;
    cint->signD_TE = 1;
    cint->signD_TM = 1;

    return 0;
}
