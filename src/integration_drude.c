/**
 * @file   integration_drude.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   September, 2016
 * @brief  Perform integration for Drude planes
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

#include "floattypes.h"
#include "utils.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration_drude.h"
#include "gausslaguerre.h"

void integrate_gauss_laguerre(casimir_t *self,
                              casimir_integrals_t *cint,
                              int l1, int l2, int m, double nT);

void integrands_drude(float64 x, integrands_drude_t *integrands,
                      struct integ_context* context,
                      unsigned int index,
                      unsigned int iteration);

/** @brief Calculate integrals A,B,C,D including prefactor Lambda vor Drude metals
 *
 * This function calculates
 *    Lambda(l1,l2,m)*A_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*B_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*C_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*D_(l1,l2)^(m)
 * for Drude metals.
 *
 * @param [in]  Drude integration object
 * @param [in]  l1   \f$\ell_1\f$
 * @param [in]  l2   \f$\ell_2\f$
 * @param [out] cint logarithms of values and signs of integrals
 */
void casimir_integrate_drude(integration_drude_t* int_drude, int l1, int l2, casimir_integrals_t* cint)
{
    integrate_gauss_laguerre(int_drude->casimir, cint, l1, l2, int_drude->m, int_drude->nT);
}


void casimir_integrate_drude_init(casimir_t* self, integration_drude_t* int_drude, double nT, int m, int lmax)
{
    int_drude->casimir = self;
    int_drude->nT      = nT;
    int_drude->m       = m;
    int_drude->lmax    = lmax;
}


void casimir_integrate_drude_free(integration_drude_t* int_drude)
{
}


/*
 * ******************** Gauss-Laguerre Integration ********************
 */
static void integrands_drude_laguerre(float64 x, integrands_drude_t *integrands,
                                      casimir_t *self, double nT, int l1, int l2, int m)
{
    plm_combination_t comb;
    const float64 tau        = 2*nT;
    const float64 k          = sqrt64(pow_2(x)/4 + nT*x);
    const float64 log_factor = log64(pow_2(x)+2*tau*x);
    float64 r_TE, r_TM;

    casimir_rp(self, nT, k, &r_TE, &r_TM);
    const float64 lnr_TE = log64(-r_TE);
    const float64 lnr_TM = log64(r_TM);

    plm_PlmPlm(l1, l2, m, 1+x/tau, &comb);

    const float64 A = comb.lnPl1mPl2m - log_factor;
    integrands->lnA_TE = lnr_TE + A;
    integrands->lnA_TM = lnr_TM + A;
    integrands->sign_A = comb.sign_Pl1mPl2m;

    const float64 B = comb.lndPl1mdPl2m + log_factor;
    integrands->lnB_TE = lnr_TE + B;
    integrands->lnB_TM = lnr_TM + B;
    integrands->sign_B = comb.sign_dPl1mdPl2m;

    const float64 C = comb.lnPl1mdPl2m;
    integrands->lnC_TE = lnr_TE + C;
    integrands->lnC_TM = lnr_TM + C;
    integrands->sign_C = comb.sign_Pl1mdPl2m;

    const float64 D = comb.lndPl1mPl2m;
    integrands->lnD_TE = lnr_TE + D;
    integrands->lnD_TM = lnr_TM + D;
    integrands->sign_D = comb.sign_dPl1mPl2m;
}


/** @brief Calculate integrals A,B,C,D including prefactor Lambda vor Drude metals
 *
 * This function calculates
 *    Lambda(l1,l2,m)*A_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*B_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*C_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*D_(l1,l2)^(m)
 * for Drude metals.
 *
 * @param [in]  self Casimir object
 * @param [out] cint logarithms of values and signs of integrals
 * @param [in]  l1   \f$\ell_1\f$
 * @param [in]  l2   \f$\ell_2\f$
 * @param [in]  m    \f$m\f$
 * @param [in]  nT   \f$nT\f$
 */
void integrate_gauss_laguerre(casimir_t *self, casimir_integrals_t *cint, int l1, int l2, int m, double nT)
{
    integrands_drude_t integrand;
    const float64 tau = 2*nT;
    const float64 ln_tau = log64(tau);
    const float64 ln_Lambda = casimir_lnLambda(l1, l2, m, NULL); /* sign: -1 */
    float64 prefactor;
    float64 *xk, *ln_wk;
    /* XXX use order 50 at the moment XXX */
    const int N = gausslaguerre_nodes_weights(50, &xk, &ln_wk);

    /* allocate space for signs_A, signs_B, signs_C, signs_D */
    sign_t *signs_ABCD = xmalloc(4*N*sizeof(sign_t));
    sign_t *signs_A = signs_ABCD;
    sign_t *signs_B = signs_A+1*N;
    sign_t *signs_C = signs_A+2*N;
    sign_t *signs_D = signs_A+3*N;

    /* allocate space for lnA_TE, lnA_TM, lnB_TE, lnB_TM, lnC_TE, lnC_TM,
     * lnD_TE, lnD_TM */
    float64 *ln_ABCD = xmalloc(4*2*N*sizeof(integrands_drude_t));
    float64 *lnA_TE  = ln_ABCD;
    float64 *lnA_TM  = ln_ABCD + 1*N;
    float64 *lnB_TE  = ln_ABCD + 2*N;
    float64 *lnB_TM  = ln_ABCD + 3*N;
    float64 *lnC_TE  = ln_ABCD + 4*N;
    float64 *lnC_TM  = ln_ABCD + 5*N;
    float64 *lnD_TE  = ln_ABCD + 6*N;
    float64 *lnD_TM  = ln_ABCD + 7*N;

    for(int i = 0; i < N; i++)
    {
        integrands_drude_laguerre(xk[i], &integrand, self, nT, l1, l2, m);

        lnA_TE[i]  = ln_wk[i] + integrand.lnA_TE;
        lnA_TM[i]  = ln_wk[i] + integrand.lnA_TM;
        signs_A[i] = integrand.sign_A;

        lnB_TE[i]  = ln_wk[i] + integrand.lnB_TE;
        lnB_TM[i]  = ln_wk[i] + integrand.lnB_TM;
        signs_B[i] = integrand.sign_B;

        lnC_TE[i]  = ln_wk[i] + integrand.lnC_TE;
        lnC_TM[i]  = ln_wk[i] + integrand.lnC_TM;
        signs_C[i] = integrand.sign_C;

        lnD_TE[i]  = ln_wk[i] + integrand.lnD_TE;
        lnD_TM[i]  = ln_wk[i] + integrand.lnD_TM;
        signs_D[i] = integrand.sign_D;
    }


    /* B */
    prefactor = ln_Lambda -tau-3*ln_tau; /* exp(-tau)/tau³ */
    cint->lnB_TE = prefactor + logadd_ms(lnB_TE, signs_B, N, &cint->signB_TE);
    cint->lnB_TM = prefactor + logadd_ms(lnB_TM, signs_B, N, &cint->signB_TM);

    cint->signB_TM = -MPOW(l2+m+1) * cint->signB_TM;
    cint->signB_TE = +MPOW(l2+m+1) * cint->signB_TE;


    if(m > 0)
    {
        const float64 log_m = log(m);

        /* A */
        prefactor = ln_Lambda + 2*log_m+ln_tau-tau; /* m²*tau*exp(-tau) */
        cint->lnA_TE = prefactor + logadd_ms(lnA_TE, signs_A, N, &cint->signA_TE);
        cint->lnA_TM = prefactor + logadd_ms(lnA_TM, signs_A, N, &cint->signA_TM);

        /* r_TE is negative, r_TM is positive and Lambda(l1,l2,m) is negative.
           => TM negative sign, TE positive sign */
        cint->signA_TM = -MPOW(l2+m) * cint->signA_TM;
        cint->signA_TE = +MPOW(l2+m) * cint->signA_TE;


        /* C */
        prefactor = ln_Lambda + log_m-tau-ln_tau; /* m*exp(-tau)/tau */
        cint->lnC_TE = prefactor + logadd_ms(lnC_TE, signs_C, N, &cint->signC_TE);
        cint->lnC_TM = prefactor + logadd_ms(lnC_TM, signs_C, N, &cint->signC_TM);

        cint->signC_TM = -MPOW(l2+m) * cint->signC_TM;
        cint->signC_TE = +MPOW(l2+m) * cint->signC_TE;


        /* D */
        /* prefactor is identical to C */
        cint->lnD_TE = prefactor + logadd_ms(lnD_TE, signs_D, N, &cint->signD_TE);
        cint->lnD_TM = prefactor + logadd_ms(lnD_TM, signs_D, N, &cint->signD_TM);

        cint->signD_TM = -MPOW(l2+m+1) * cint->signD_TM;
        cint->signD_TE = +MPOW(l2+m+1) * cint->signD_TE;
    }
    else
    {
        cint->lnA_TM = cint->lnA_TE = -INFINITY;
        cint->signA_TM = cint->signA_TE = +1;

        cint->lnC_TM = cint->lnC_TE = -INFINITY;
        cint->signC_TM = cint->signC_TE = +1;

        cint->lnD_TM = cint->lnD_TE = -INFINITY;
        cint->signD_TM = cint->signD_TE = +1;
    }

    xfree(ln_ABCD);
    xfree(signs_ABCD);
}
