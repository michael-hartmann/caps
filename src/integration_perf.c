/**
 * @file   integration_perf.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   November, 2015
 * @brief  Perform integration for perfect reflectors
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "edouble.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration_perf.h"
#include "utils.h"

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
static edouble polyintegrate(edouble p[], const int len_p, const int offset, const edouble tau)
{
    int k;
    const edouble log_tau = loge(tau);
    edouble list[len_p];

    for(k = offset; k < len_p+offset; k++)
        list[k-offset] = lgammae(k+1)-(k+1)*log_tau+p[k-offset];

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
static void log_polymult(edouble p1[], const int len_p1, edouble p2[], const int len_p2, edouble p[])
{
    int k, len = len_p1+len_p2-1;
    edouble temp[len];

    for(k = 0; k < len; k++)
    {
        int i, elems = 0;

        for(i = 0; i < MIN(k+1, len_p1); i++)
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
 *      \mathcal{I}_\nu^{2m}(tau) = (-1)^m \int_0^\infty \mathrm{d}z \, exp{(-\tau z)} (z^2+2z)^{-1} * \mathrm{Plm}_\nu^{2m}(1+z).
 * \f]
 * The value of the integral is always positive.
 *
 * @param [in] nu
 * @param [in] m2
 * @param [in] tau
 * @retval logarithm of value of integral
 */
edouble casimir_integrate_I(integration_perf_t *self, int nu)
{
    const int m = self->m;
    edouble v = self->cache_I[nu];

    if(isnan(v))
    {
        int k;
        double tau = 2*self->nT;

        edouble p1[m];        /* polynom (z+2)^(m-1) */
        edouble p2[nu+1-2*m]; /* polynom d^(2m)/dz^(2m) P_(nu)(1+z) */
        edouble p[-m+nu];     /* polynom p1*p2 */

        /* Every monom of both polynoms is positive. So we can save the
         * logarithms of the coefficients. */

        for(k = 0; k <= m-1; k++)
            p1[k] = lgammae(m)-lgammae(k+1)-lgammae(m-k)+(m-1-k)*LOG2;

        for(k = 2*m; k <= nu; k++)
            p2[k-2*m] = lgammae(k+nu+1)-lgammae(k+1)-lgammae(k-2*m+1)-lgammae(-k+nu+1)-k*LOG2;

        log_polymult(p1, m, p2, nu+1-2*m, p); /* len: nu-m */

        v = self->cache_I[nu] = polyintegrate(p, -m+nu, m-1, tau);
        TERMINATE(!isfinite(v) || !isnan(v), "I=%Lg, nu=%d, m=%d\n", v, nu, m);
    }

    return v;
}

edouble casimir_integrate_K(integration_perf_t *self, const int l1, const int l2, sign_t *sign)
{
    const int index = l1*(self->lmax+1)+l2;
    edouble v = self->cache_K[index];

    if(isnan(v))
    {
        int q;
        const int m = self->m;
        const int qmax       = gaunt_qmax(l1,l2,m);
        const edouble log_a0 = gaunt_log_a0(l1,l2,m);
        const int elems = MAX(0,1+qmax);
        edouble *a = NULL;
        sign_t *signs = NULL;

        a     = xmalloc(elems*sizeof(edouble));
        signs = xmalloc(elems*sizeof(sign_t));

        gaunt(l1, l2, m, a);
        
        for(q = 0; q <= qmax; q++)
        {
            signs[q] = copysigne(1, a[q]);
            a[q] = loge(a[q]) + casimir_integrate_I(self, l1+l2-2*q);
        }

        v = self->cache_K[index] = log_a0+logadd_ms(a, signs, elems, sign);

        xfree(a);
        xfree(signs);
    }

    return v;
}


void casimir_integrate_perf_init(integration_perf_t *self, double nT, int m, int lmax)
{
    int i, elems_I, elems_K;;

    self->tau  = 2*nT;
    self->nT   = nT;
    self->lmax = lmax;
    self->m    = m;

    elems_I = 2*(lmax+1);
    self->cache_I = xmalloc(elems_I*sizeof(edouble));
    for(i = 0; i < elems_I; i++)
        self->cache_I[i] = NAN;

    elems_K = pow_2(lmax+1);
    self->cache_K = xmalloc(elems_K*sizeof(edouble));
    for(i = 0; i < elems_K; i++)
        self->cache_K[i] = NAN;
}

void casimir_integrate_perf_free(integration_perf_t *self)
{
    if(self->cache_I != NULL)
    {
        xfree(self->cache_I);
        self->cache_I = NULL;
    }
    if(self->cache_K != NULL)
    {
        xfree(self->cache_K);
        self->cache_K = NULL;
    }
}

/** @brief Evaluate integrals A,B,C,D
 *
 * Note that we also include the factor \f$\Lambda_{\ell_1,\ell_2}^{(m)}\f$.
 */
void casimir_integrate_perf(integration_perf_t *self, int l1, int l2, casimir_integrals_t *cint)
{
    const int m  = self->m;
    const double nT = self->nT;
    const edouble lnLambda = casimir_lnLambda(l1, l2, m, NULL);
    edouble value, v[4];
    sign_t sign, signs[4];

    /* A */
    {
        value = casimir_integrate_K(self, l1, l2, &sign);

        cint->lnA_TM = cint->lnA_TE = lnLambda+log_A0(m,nT)+value;
        cint->signA_TM = -sign*sign_A0(l2,m,TM); /* - because Lambda(l1,l2,m) is negative */
        cint->signA_TE = -sign*sign_A0(l2,m,TE);
    }

    /* B */
    {
        v[0] = v[1] = v[2] = v[3] = -INFINITY;

        if((l1-1) >= m && (l2-1) >= m)
        {
            v[0] = loge(l1+1)+loge(l1+m)+loge(l2+1)+loge(l2+m) + casimir_integrate_K(self, l1-1, l2-1, &signs[0]);
        }

        if((l1-1) >= m)
        {
            v[1] = loge(l1+1)+loge(l1+m)+loge(l2)+loge(l2-m+1) + casimir_integrate_K(self, l1-1, l2+1, &signs[1]);
            signs[1] *= -1;
        }

        if((l2-1) >= m)
        {
            v[2] = loge(l1)+loge(l1-m+1)+loge(l2+1)+loge(l2+m) + casimir_integrate_K(self, l1+1, l2-1, &signs[2]);
            signs[2] *= -1;
        }

        v[3] = loge(l1)+loge(l1-m+1)+loge(l2)+loge(l2-m+1) + casimir_integrate_K(self, l1+1, l2+1, &signs[3]);

        /* add */
        value = logadd_ms(v, signs, 4, &sign);
        
        cint->lnB_TM = cint->lnB_TE = lnLambda+log_B0(m,nT)-loge(2*l1+1)-loge(2*l2+1)+value;
        cint->signB_TM = -sign*sign_B0(l2,m,TM);
        cint->signB_TE = -sign*sign_B0(l2,m,TE);
    }

    /* C */
    {
        v[0] = v[1] = -INFINITY;

        if((l2-1) >= m)
        {
            v[0] = loge(l2+1)+loge(l2+m) + casimir_integrate_K(self, l1, l2-1, &signs[0]);
        }

        v[1] = loge(l2)+loge(l2-m+1) + casimir_integrate_K(self, l1, l2+1, &signs[1]);
        signs[1] *= -1;

        /* add */
        value = logadd_ms(v, signs, 2, &sign);

        cint->lnC_TM = cint->lnD_TE = lnLambda+log_C0(m,nT)-loge(2*l2+1) + value;
        cint->signC_TM = -sign*sign_C0(l2,m,TM);
        cint->signC_TE = -sign*sign_C0(l2,m,TE);
    }

    /* D */
    {
        v[0] = v[1] = -INFINITY;

        if((l1-1) >= m)
        {
            v[0] = loge(l1+1)+loge(l1+m) + casimir_integrate_K(self, l1-1, l2, &signs[0]);
        }

        v[1] = loge(l1)+loge(l1-m+1) + casimir_integrate_K(self, l1+1, l2, &signs[1]);
        signs[1] *= -1;

        /* add */
        value = logadd_ms(v, signs, 2, &sign);

        cint->lnC_TM = cint->lnD_TE = lnLambda+log_D0(m,nT)-loge(2*l1+1) + value;
        cint->signC_TM = -sign*sign_D0(l2,m,TM);
        cint->signC_TE = -sign*sign_D0(l2,m,TE);
    }

    TERMINATE(!isfinite(cint->lnA_TM) || !isnan(cint->lnA_TM), "lnA=%Lg, l1=%d,l2=%d,m=%d,nT=%g", cint->lnA_TM,l1,l2,m,nT);
    TERMINATE(!isfinite(cint->lnB_TM) || !isnan(cint->lnB_TM), "lnB=%Lg, l1=%d,l2=%d,m=%d,nT=%g", cint->lnB_TM,l1,l2,m,nT);
    TERMINATE(!isfinite(cint->lnC_TM) || !isnan(cint->lnC_TM), "lnC=%Lg, l1=%d,l2=%d,m=%d,nT=%g", cint->lnC_TM,l1,l2,m,nT);
    TERMINATE(!isfinite(cint->lnD_TM) || !isnan(cint->lnD_TM), "lnD=%Lg, l1=%d,l2=%d,m=%d,nT=%g", cint->lnD_TM,l1,l2,m,nT);
}
