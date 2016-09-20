/**
 * @file   integration_perf.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   November, 2015
 * @brief  Perform integration for perfect reflectors
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "floattypes.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration_perf.h"
#include "utils.h"

static void casimir_integrate_perf_m0(integration_perf_t *self, int l1, int l2, casimir_integrals_t *cint);


/** @brief Integrate polynomial x^offset*p*exp(-tau*x)
 *
 * This function evaluates the integral
 * \f[
 *      \inf_0^\infty \mathrm{d}x \, x^\mathrm{offset} p(x) \exp{(-\tau x)},
 * \f]
 * where \f$p(x)\f$ is a polynomial. The coefficients of p are given as
 * logarithms and are assumed to be positive.
 *
 * @param [in] p coefficients of polynomial
 * @param [in] len_p number of coefficients of polynomial p
 * @param [in] offset parameter
 * @param [in] tau parameter
 * @retval logarithm of value of integral
 */
static float80 polyintegrate(float80 p[], const int len_p, const int offset, const float80 tau)
{
    const float80 log_tau = log80(tau);
    float80 list[len_p];

    for(int k = offset; k < len_p+offset; k++)
        list[k-offset] = lgamma80(k+1)-(k+1)*log_tau+p[k-offset];

    return logadd_m(list, len_p);
}


/** @brief Multiply two polynomials
 *
 * Multiply the polynomials p1 (len_p1 coefficients) and p2 (len_p2
 * coefficients) and store result in p.
 *
 * @param [in]  p1 coefficients of polynomial p1
 * @param [in]  len_p1 number of coefficients of polynomial p1
 * @param [in]  p2 coefficients of polynomial p1
 * @param [in]  len_p2 number of coefficients of polynomial p2
 * @param [out] p polynomial p1*p2
 */
static void log_polymult(float80 p1[], const int len_p1, float80 p2[], const int len_p2, float80 p[])
{
    int len = len_p1+len_p2-1;
    float80 temp[len];

    for(int k = 0; k < len; k++)
    {
        int elems = 0;

        for(int i = 0; i < MIN(k+1, len_p1); i++)
        {
            int j = k-i;
            if(i < len_p1 && j < len_p2)
                temp[elems++] = p1[i]+p2[j];
        }

        if(elems)
            p[k] = logadd_m(temp, elems);
        else
            p[k] = -INFINITY;
    }
}


/** @brief Evaluate integral I
 *
 * This function evaluates the integral
 * \f[
 *      \mathcal{I}_\nu^{2m}(\tau=2nT) = (-1)^m \int_0^\infty \mathrm{d}z \, exp{(-\tau z)} (z^2+2z)^{-1} * \mathrm{Plm}_\nu^{2m}(1+z).
 * \f]
 *
 * The parameter m will be determined automatically. Usually m is the same
 * value of m that was given to casimir_integrate_perf_init, except for m=0,
 * where m=2 is used.
 *
 * The value of the integral is always positive.
 *
 * @param [in] self integration object
 * @param [in] nu
 * @retval logarithm of value of integral
 */
float80 casimir_integrate_perf_I(integration_perf_t *self, int nu)
{
    float80 v = self->cache_I[nu];

    if(isnan(v))
    {
        const double tau = 2*self->nT;
        const int m = (self->m != 0) ? self->m : 2;

        /* polynom (z+2)^(m-1) */
        float80 *p1 = xmalloc(m*sizeof(float80));
        /* polynom d^(2m)/dz^(2m) P_(nu)(1+z) */
        float80 *p2 = xmalloc((nu+1-2*m)*sizeof(float80));
        /* polynom p1*p2 */
        float80 *p = xmalloc((nu-m)*sizeof(float80));

        /* Every monom of both polynoms is positive. So we can save the
         * logarithms of the coefficients. */

        for(int k = 0; k <= m-1; k++)
            p1[k] = lgamma80(m)-lgamma80(k+1)-lgamma80(m-k)+(m-1-k)*M_LOG2;

        for(int k = 2*m; k <= nu; k++)
            p2[k-2*m] = lgamma80(k+nu+1)-lgamma80(k+1)-lgamma80(k-2*m+1)-lgamma80(-k+nu+1)-k*M_LOG2;

        log_polymult(p1, m, p2, nu+1-2*m, p); /* len: nu-m */

        v = self->cache_I[nu] = polyintegrate(p, -m+nu, m-1, tau);
        TERMINATE(!isfinite(v), "I=%Lg, nu=%d, m=%d\n", v, nu, m);

        xfree(p1);
        xfree(p2);
        xfree(p);
    }

    return v;
}


/** @brief Evaluate integral K
 *
 * This function evaluates the integral
 * \f[
 *      \mathcal{K}_{\ell_1,\ell_2}^m(\tau=2nT) = a_0 \sum_{q=0}^{q_\mathrm{max}} \tilde a_q \mathcal{I}_\nu^{2m}(\tau=2nT)
 * \f]
 *
 * The sign is stored in sign.
 *
 * @param [in] self integration object
 * @param [in] l1 parameter l1
 * @param [in] l2 parameter l2
 * @param [out] sign sign of integral

 * @retval logarithm of value of integral
 */
float80 casimir_integrate_K(integration_perf_t *self, const int l1, const int l2, sign_t *sign)
{
    const int index = l1*(self->lmax+2)+l2;
    float80 v = self->cache_K[index];

    if(isnan(v))
    {
        const int m = (self->m != 0) ? self->m : 2;
        const int qmax       = gaunt_qmax(l1,l2,m);
        const float80 log_a0 = gaunt_log_a0(l1,l2,m);
        const int elems = MAX(0,1+qmax);

        float80 *a    = xmalloc(elems*sizeof(float80));
        sign_t *signs = xmalloc(elems*sizeof(sign_t));

        gaunt(l1, l2, m, a);

        for(int q = 0; q <= qmax; q++)
        {
            signs[q] = copysign80(1, a[q]);
            a[q] = log80(fabs80(a[q])) + casimir_integrate_perf_I(self, l1+l2-2*q);
        }

        v = self->cache_K[index] = log_a0+logadd_ms(a, signs, elems, sign);
        self->cache_K_signs[index] = *sign;

        xfree(a);
        xfree(signs);

        TERMINATE(isnan(v), "casimir_integrate_K l1=%d, l2=%d, m=%d, %Lg\n", l1, l2, m, v);

        return v;
    }
    else
    {
        *sign = self->cache_K_signs[index];
        return v;
    }
}


/** @brief Initialize integration_perf_t object
 *
 * This function initializes the integration_perf_t object self for Matsubara frequency xi=nT.
 *
 * @param [out] self integration object
 * @param [in] nT Matsubara frequency, xi=nT
 * @param [in] m angular momentum z-axis
 * @param [out] lmax maximum value of l1,l2
 */
void casimir_integrate_perf_init(integration_perf_t *self, double nT, int m, int lmax)
{
    self->nT   = nT;
    self->lmax = lmax;
    self->m    = m;

    /* allocate memory and initialize chache for I */
    const int elems_I = 2*(lmax+2);
    self->cache_I = xmalloc(elems_I*sizeof(float80));
    for(int i = 0; i < elems_I; i++)
        self->cache_I[i] = NAN;

    /* allocate memory and initialize chache for K */
    const int elems_K = pow_2(lmax+2);
    self->cache_K_signs = xmalloc(elems_K*sizeof(sign_t));
    self->cache_K = xmalloc(elems_K*sizeof(float80));
    for(int i = 0; i < elems_K; i++)
        self->cache_K[i] = NAN;
}


/** @brief Free integration_perf_t object
 *
 * This function frees the integration_perf_t object self.
 *
 * @param [in,out] self integration object
 */
void casimir_integrate_perf_free(integration_perf_t *self)
{
    xfree(self->cache_I);
    xfree(self->cache_K);
    xfree(self->cache_K_signs);

    self->cache_I = NULL;
    self->cache_K = NULL;
    self->cache_K_signs = NULL;
}

/* integrate for m=0 */
static void casimir_integrate_perf_m0(integration_perf_t *self, int l1, int l2, casimir_integrals_t *cint)
{
    float80 log_B;
    sign_t sign_B;
    const float80 lnLambda = casimir_lnLambda(l1, l2, 0, NULL);

    sign_t signs_B[4] = {1,1,1,1};
    float80 log_Bn[4] = { -INFINITY, -INFINITY, -INFINITY, -INFINITY };

    if((l1-1) >= 2 && (l2-1) >= 2)
        log_Bn[0] = casimir_integrate_K(self, l1-1, l2-1, &signs_B[0]);

    if((l1-1) >= 2)
    {
        log_Bn[1] = casimir_integrate_K(self, l1-1, l2+1, &signs_B[1]);
        signs_B[1] *= -1;
    }

    if((l2-1) >= 2)
    {
        log_Bn[2] = casimir_integrate_K(self, l1+1, l2-1, &signs_B[2]);
        signs_B[2] *= -1;
    }

    log_Bn[3] = casimir_integrate_K(self, l1+1, l2+1, &signs_B[3]);

    log_B = logadd_ms(log_Bn, signs_B, 4, &sign_B) - log80(2*l1+1) - log80(2*l2+1);

    cint->lnB_TM = cint->lnB_TE = lnLambda+log_B0(m,self->nT)+log_B;
    cint->signB_TM = sign_B*sign_B0(l2,m,TM);
    cint->signB_TE = sign_B*sign_B0(l2,m,TE);

    /* set, A,C,D = +0 */
    cint->lnA_TE = cint->lnA_TM   = -INFINITY;
    cint->signA_TM = cint->signA_TE = 1;

    cint->lnC_TE = cint->lnC_TM   = -INFINITY;
    cint->signC_TM = cint->signC_TE = 1;

    cint->lnD_TE = cint->lnD_TM   = -INFINITY;
    cint->signD_TM = cint->signD_TE = 1;
}

/** @brief Evaluate integrals A,B,C,D
 *
 * This function evaluates the integral A,B,C,D for l1,l2. The results are stored in cint.
 *
 * Note that we also include the factor \f$\Lambda_{\ell_1,\ell_2}^{(m)}\f$.
 *
 * @param [in] self integration object
 * @param [in] l1 parameter l1
 * @param [in] l2 parameter l2
 * @param [out] cint values of integrals
 */
void casimir_integrate_perf(integration_perf_t *self, int l1, int l2, casimir_integrals_t *cint)
{
    const int m  = self->m;
    const float80 nT = self->nT;
    const float80 lnLambda = casimir_lnLambda(l1, l2, m, NULL);
    float80 value, v[4];
    sign_t sign, signs[4];

    if(m == 0)
    {
        casimir_integrate_perf_m0(self, l1, l2, cint);
        TERMINATE(!isfinite(cint->lnB_TM), "lnB=%Lg, l1=%d,l2=%d,m=%d,nT=%Lg", cint->lnB_TM,l1,l2,m,nT);

        return;
    }

    /* m != 0 */

    /* A */
    {
        value = casimir_integrate_K(self, l1, l2, &sign);

        cint->lnA_TM = cint->lnA_TE = lnLambda+log_A0(m,nT)+value;
        cint->signA_TM = sign*sign_A0(l2,m,TM); /* - because Lambda(l1,l2,m) is negative */
        cint->signA_TE = sign*sign_A0(l2,m,TE);
    }

    /* B */
    {
        v[0] = v[1] = v[2] = v[3] = -INFINITY;
        signs[0] = signs[1] = signs[2] = signs[3] = 1;

        if((l1-1) >= m && (l2-1) >= m)
            v[0] = log80(l1+1)+log80(l1+m)+log80(l2+1)+log80(l2+m) + casimir_integrate_K(self, l1-1, l2-1, &signs[0]);

        if((l1-1) >= m)
        {
            v[1] = log80(l1+1)+log80(l1+m)+log80(l2)+log80(l2-m+1) + casimir_integrate_K(self, l1-1, l2+1, &signs[1]);
            signs[1] *= -1;
        }

        if((l2-1) >= m)
        {
            v[2] = log80(l1)+log80(l1-m+1)+log80(l2+1)+log80(l2+m) + casimir_integrate_K(self, l1+1, l2-1, &signs[2]);
            signs[2] *= -1;
        }

        v[3] = log80(l1)+log80(l1-m+1)+log80(l2)+log80(l2-m+1) + casimir_integrate_K(self, l1+1, l2+1, &signs[3]);

        /* add */
        value = logadd_ms(v, signs, 4, &sign);

        cint->lnB_TM = cint->lnB_TE = lnLambda+log_B0(m,nT)-log80(2*l1+1)-log80(2*l2+1)+value;
        cint->signB_TM = sign*sign_B0(l2,m,TM);
        cint->signB_TE = sign*sign_B0(l2,m,TE);
    }

    /* C */
    {
        v[0] = v[1] = -INFINITY;
        signs[0] = signs[1] = 1;

        if((l2-1) >= m)
            v[0] = log80(l2+1)+log80(l2+m) + casimir_integrate_K(self, l1, l2-1, &signs[0]);

        v[1] = log80(l2)+log80(l2-m+1) + casimir_integrate_K(self, l1, l2+1, &signs[1]);
        signs[1] *= -1;

        /* add */
        value = logadd_ms(v, signs, 2, &sign);

        cint->lnC_TM = cint->lnC_TE = lnLambda+log_C0(m,nT)-log80(2*l2+1) + value;
        cint->signC_TM = sign*sign_C0(l2,m,TM);
        cint->signC_TE = sign*sign_C0(l2,m,TE);
    }

    /* D */
    {
        v[0] = v[1] = -INFINITY;
        signs[0] = signs[1] = 1;

        if((l1-1) >= m)
            v[0] = log80(l1+1)+log80(l1+m) + casimir_integrate_K(self, l1-1, l2, &signs[0]);

        v[1] = log80(l1)+log80(l1-m+1) + casimir_integrate_K(self, l1+1, l2, &signs[1]);
        signs[1] *= -1;

        /* add */
        value = logadd_ms(v, signs, 2, &sign);

        cint->lnD_TM = cint->lnD_TE = lnLambda+log_D0(m,nT)-log80(2*l1+1) + value;
        cint->signD_TM = sign*sign_D0(l2,m,TM);
        cint->signD_TE = sign*sign_D0(l2,m,TE);
    }

    TERMINATE(!isfinite(cint->lnA_TM), "lnA=%Lg, l1=%d,l2=%d,m=%d,nT=%Lg", cint->lnA_TM,l1,l2,m,nT);
    TERMINATE(!isfinite(cint->lnB_TM), "lnB=%Lg, l1=%d,l2=%d,m=%d,nT=%Lg", cint->lnB_TM,l1,l2,m,nT);
    TERMINATE(!isfinite(cint->lnC_TM), "lnC=%Lg, l1=%d,l2=%d,m=%d,nT=%Lg", cint->lnC_TM,l1,l2,m,nT);
    TERMINATE(!isfinite(cint->lnD_TM), "lnD=%Lg, l1=%d,l2=%d,m=%d,nT=%Lg", cint->lnD_TM,l1,l2,m,nT);
}
