/**
 * @file   sfunc.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   April, 2015
 * @brief  various functions implementing mostly special functions
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "floattypes.h"
#include "sfunc.h"
#include "utils.h"

/**
* @name Add numbers given as logarithms
*/
/*@{*/


/**
 * @brief Sum len array elements of input
 *
 * Use Kahan summation algorithm to reduce the numerical error obtained by
 * adding all array elements of input.
 *
 * See: https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 *
 * @param [in] input array
 * @param [in] len length of array
 * @retval sum
 */
float80 kahan_sum(float80 input[], size_t len)
{
    float80 sum = 0;
    float80 c = 0; /* running compensation for lost low-order bits */

    for(size_t i = 0; i < len; i++)
    {
        float80 y = input[i] - c;
        float80 t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}



/**
 * @brief Calculate relative difference between two numbers
 *
 * This function computes \f$\log(|a-b|/\mathrm{max}(a,b))\f$, where
 * \f$\mathrm{log_a} = \log a\f$ and \f$\mathrm{log_b} = \log b\f$.
 *
 * @param [in] log_a number
 * @param [in] log_b number
 * @return log_diff_relative relative difference
 */
inline float80 logdiff_rel(const float80 log_a, const float80 log_b)
{
    if(isinf(log_a) && log_a < 0)
        return 0;
    else if(isinf(log_b) && log_b < 0)
        return 0;

    if(log_a > log_b)
        return log1p80(-exp80(log_b-log_a));
    else
        return log1p80(-exp80(log_a-log_b));
}

/**
 * @brief Add two numbers given by their logarithms.
 *
 * Both numbers are assumed to be nonnegative.
 *
 * @param [in] log_a number
 * @param [in] log_b number
 * @return log_sum \f$\log{\left[\exp{(\mathrm{log\_a})}+\exp{(log\_b)}\right]}\f$
 */
inline float80 logadd(const float80 log_a, const float80 log_b)
{
    if(isinf(log_a) && log_a < 0)
        return log_b;
    else if(isinf(log_b) && log_b < 0)
        return log_a;

    if(log_a > log_b)
        return log_a + log1p80(exp80(log_b-log_a));
    else
        return log_b + log1p80(exp80(log_a-log_b));
}


/**
 * @brief Add two numbers given by their logarithms, respecting signs
 *
 * @param [in]  log_a number
 * @param [in]  sign_a sign of a
 * @param [in]  log_b number
 * @param [in]  sign_b sign of b
 * @param [out] sign sign of result
 * @return log_sum \f$\log{\left( \mathrm{sign\_a}\cdot\exp{(\mathrm{log\_a})}+\mathrm{sign\_b} \cdot \exp{(log\_b)} \right)}\f$
 */
float80 logadd_s(const float80 log_a, const sign_t sign_a, const float80 log_b, const sign_t sign_b, sign_t *sign)
{
    if(isinf(log_a) && log_a < 0)
    {
        *sign = sign_b;
        return log_b;
    }
    else if(isinf(log_b) && log_b < 0)
    {
        *sign = sign_a;
        return log_a;
    }

    if(log_a > log_b)
    {
        *sign = sign_a;
        return log_a + log1p80(sign_a*sign_b*exp80(log_b-log_a));
    }
    else
    {
        *sign = sign_b;
        return log_b + log1p80(sign_a*sign_b*exp80(log_a-log_b));
    }
}


/**
 * @brief Add numbers given by their logarithms.
 *
 * len numbers in list will be added. The numbers are assumed to be positive.
 *
 * @param [in]  list array of numbers given by logarithm
 * @param [in]  len length of list
 * @return log_sum \f$\log{\sum_{i=1}^\mathrm{len} \mathrm{sign\_i}\cdot\exp{(\mathrm{log\_i})}}\f$
 */
inline float80 logadd_m(const float80 list[], const int len)
{
    float80 sum;
    float80 max = list[0];

    for(int i = 1; i < len; i++)
        if(list[i] > max)
            max = list[i];

    sum = exp80(list[0]-max);
    for(int i = 1; i < len; i++)
        sum += exp80(list[i]-max);

    return max + log80(sum);
}


inline float80 logadd_ms(const float80 list[], const sign_t signs[], const int len, sign_t *sign)
{
    float80 sum;
    float80 max = list[0];

    for(int i = 1; i < len; i++)
        if(list[i] > max)
            max = list[i];

    sum = signs[0]*exp80(list[0]-max);
    for(int i = 1; i < len; i++)
        sum += signs[i]*exp80(list[i]-max);

    *sign = copysign80(1, sum);
    return max + log80(fabs80(sum));
}

/*@}*/


/**
* @name Bessel functions
*/
/*@{*/

void bessel_lnInuKnu(int nu, const float80 x, float80 *lnInu_p, float80 *lnKnu_p)
{
    const float80 logx = log80(x);
    float80 lnKnu = 0, lnKnup = log1p80(1./x);

    /* calculate Knu, Knup */
    {
        const float80 prefactor = -x+0.5*(LOGPI-LOG2-logx);

        if(nu == 0)
        {
            lnKnu  = prefactor+lnKnu;
            lnKnup = prefactor+lnKnup;
        }
        else
        {
            for(int l = 2; l <= nu+1; l++)
            {
                float80 lnKn_new = logadd(log80(2*l-1)+lnKnup-logx, lnKnu);
                lnKnu  = lnKnup;
                lnKnup = lnKn_new;
            }

            lnKnup = prefactor+lnKnup;
            lnKnu  = prefactor+lnKnu;
        }

        TERMINATE(!isfinite(lnKnup), "Couldn't calculate Bessel functions, nu=%d, x=%Lg\n", nu, x);

        if(lnKnu_p != NULL)
            *lnKnu_p = lnKnu;
    }

    if(lnInu_p != NULL)
    {
        #define an(n,nu,x) (2*((nu)+0.5+(n))/(x))

        float80 num   = an(2,nu,x)+1/an(1,nu,x);
        float80 denom = an(2,nu,x);
        float80 ratio = (an(1,nu,x)*num)/denom;
        float80 ratio_last = 0;

        for(int l = 3; 1; l++)
        {
            num   = an(l,nu,x)+1/num;
            denom = an(l,nu,x)+1/denom;
            ratio *= num/denom;

            if(ratio_last != 0 && fabs80(1.L-ratio/ratio_last) < 1e-20L)
                break;

            ratio_last = ratio;
        }

        *lnInu_p = -logx-lnKnu-logadd(lnKnup-lnKnu, -log80(ratio));
        #undef an
    }
}


float80 bessel_lnKnu(const int nu, const float80 x)
{
    float80 Knu;
    bessel_lnInuKnu(nu, x, NULL, &Knu);
    return Knu;
}


float80 bessel_lnInu(const int nu, const float80 x)
{
    float80 Inu;
    bessel_lnInuKnu(nu, x, &Inu, NULL);
    return Inu;
}

/*@}*/

/**
 * @brief Calculate double factorial \f$n!!\f$
 *
 * @param n non-negative integer
 * @return doublefactorial \f$n!!\f$
 */
float80 ln_doublefact(int n)
{
    /* see e.g. http://en.wikipedia.org/wiki/Double_factorial */
    if(n == 0 || n == 1) /* 0!! = 1!! = 0 */
        return 0;

    if(n % 2 == 0) /* even */
    {
        int k = n/2;
        return k*LOG2 + lgamma80(1+k);
    }
    else /* odd */
    {
        int k = (n+1)/2;
        return lgamma80(1+2*k) - k*LOG2 - lgamma80(1+k);
    }
}


/**
* @name Associated Legendre polynomials
*/
/*@{*/

/**
 * @brief Calculate associated Legendre polynomials for argument \f$x \ge 1\f$
 *
 * This function calculates associated legendre functions and its derivatives
 * for \f$m \ge 0\f$ and \f$x \ge 1\f$.
 *
 * Associated Legendre polynomials are defined as follows:
 * \f[
 *     P_\ell^m(x) = (-1)^m (1-x^2)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_\ell(x)
 * \f]
 * where \f$P_\ell(x)\f$ denotes a Legendre polynomial.
 *
 * As \f$P_\ell(x)\f$ are ordinary polynomials, the only problem is the term
 * \f$(1-x^2)^{m/2}\f$ when extending the domain to values of \f$x > 1\f$. We will
 * use the convention \f$\sqrt{-x} \equiv +i \sqrt{x}\f$.
 *
 * Note: Products of associated legendre polynomials with common \f$m\f$ are
 * unambiguous, because
 * \f[
 *     i^2 = (-i)^2 = -1.
 * \f]
 *
 * This function calculates the associated Legendre polynomials for
 * \f$\ell=m,...,\ell_\mathrm{max}\f$.
 *
 * @param [in] lmax maximum value of \f$\ell\f$
 * @param [in] m order
 * @param [in] x position to evaluate associated Legendre polynomial
 * @param [out] lnplm array of logarithms of values
 * @param [out] sign corressponding signs
 */
inline void plm_lnPlm_array(int lmax, int m, float80 x, float80 lnplm[], sign_t sign[])
{
    float80 logx = log80(x);

    if(m == 0)
    {
        sign[0] = +1;
        lnplm[0] = 0; // log(1)
    }
    else
    {
        sign[0]  = MPOW((int)(m/2) + m%2);
        lnplm[0] = ln_doublefact(2*m-1) + m*0.5*log80(pow_2(x)-1); // l=m,m=m
    }

    if(lmax == m)
        return;

    sign[1]  = sign[0];
    lnplm[1] = lnplm[0]+logx+log80(2*m+1); // l=m+1, m=m

    for(int l = m+2; l <= lmax; l++)
    {
        lnplm[l-m] = logadd_s(log80(2*l-1)+logx+lnplm[l-m-1], sign[l-m-1], log80(l+m-1)+lnplm[l-m-2], -sign[l-m-2], &sign[l-m]);
        lnplm[l-m]-= log80(l-m);
    }
}

/**
 * @brief Calculate \f$P_\ell^m(x)\f$
 *
 * See \ref plm_lnPlm_array.
 *
 * @param l degree
 * @param m order
 * @param x position
 * @param sign sign of result
 * @return lnPlm \f$\log\left|P_\ell^m(x)\right|\f$
 */
float80 plm_lnPlm(int l, int m, float80 x, sign_t *sign)
{
    float80 plm[l-m+1];
    sign_t  signs[l-m+1];

    plm_lnPlm_array(l, m, x, plm, signs);
    *sign = signs[l-m];

    return plm[l-m];
}

/**
 * @brief Calculate \f$P_\ell^m(x)\f$
 *
 * See \ref plm_lnPlm_array.
 *
 * @param l degree
 * @param m order
 * @param x position
 * @return Plm \f$P_\ell^m(x)\f$
 */
float80 plm_Plm(int l, int m, float80 x)
{
    sign_t sign;
    float80 value = plm_lnPlm(l, m, x, &sign);
    return sign*exp80(value);
}

/**
 * @brief Calculate \f${P_\ell^m}^\prime (x)\f$
 *
 * See \ref plm_lnPlm_array.
 *
 * @param l degree
 * @param m order
 * @param x position
 * @param sign sign of result
 * @return lndPlm \f$\log\left|{P_\ell^m}^\prime(x)\right|\f$
 */
float80 plm_lndPlm(int l, int m, float80 x, sign_t *sign)
{
    const int lmax = l+1;
    float80 plm[lmax-m+1];
    sign_t signs[lmax-m+1];

    plm_lnPlm_array(lmax, m, x, plm, signs);

    return logadd_s(log80(l-m+1)+plm[l+1-m], signs[l+1-m], log80(l+1)+log80(x)+plm[l-m], -signs[l+1-m], sign) - log80(pow_2(x)-1);
}


/**
 * @brief Calculate \f${P_\ell^m}^\prime(x)\f$
 *
 * See \ref plm_lnPlm_array.
 *
 * @param l degree
 * @param m order
 * @param x position
 * @return dPlm \f${P_\ell^m}^\prime(x)\f$
 */
float80 plm_dPlm(int l, int m, float80 x)
{
    sign_t sign;
    float80 value = plm_lndPlm(l, m, x, &sign);
    return sign*exp80(value);
}

/**
 * @brief Calculate products of associated Legendre polynomials and its derivatives
 *
 * Calculate products \f$P_{\ell_1}^m(x) P_{\ell_2}^m(x)\f$,
 * \f${P_{\ell_1}^m}^\prime(x) P_{\ell_2}^m(x)\f$ \f$P_{\ell_1}^m(x)
 * {P_{\ell_2}^m}^\prime(x)\f$ \f${P_{\ell_1}^m}^\prime(x)
 * {P_{\ell_2}^m}^\prime(x)\f$.
 *
 * See \ref plm_lnPlm_array.
 *
 * @param [in]  l1 \f$\ell_1\f$
 * @param [in]  l2 \f$\ell_2\f$
 * @param [in]  m \f$m\f$
 * @param [in]  x \f$x\f$
 * @param [out] res save results in this struct
 */
void plm_PlmPlm(int l1, int l2, int m, float80 x, plm_combination_t *res)
{
    const int lmax = MAX(l1,l2)+1;
    float80 lnPlm[lmax-m+1];
    sign_t signs[lmax-m+1];
    float80 logx = log80(x);
    float80 logx2m1 = log80(pow_2(x)-1);
    float80 lnPl1m, lnPl2m, lndPl1m, lndPl2m;
    sign_t sign_Pl1m, sign_Pl2m, sign_dPl1m, sign_dPl2m;
    sign_t common_sign = MPOW(m%2);

    plm_lnPlm_array(lmax, m, x, lnPlm, signs);

    lnPl1m    = lnPlm[l1-m];
    sign_Pl1m = signs[l1-m];
    lnPl2m    = lnPlm[l2-m];
    sign_Pl2m = signs[l2-m];

    lndPl1m = logadd_s(log80(l1-m+1)+lnPlm[l1+1-m], signs[l1+1-m], log80(l1+1)+logx+lnPlm[l1-m], -signs[l1+1-m], &sign_dPl1m) - logx2m1;
    lndPl2m = logadd_s(log80(l2-m+1)+lnPlm[l2+1-m], signs[l2+1-m], log80(l2+1)+logx+lnPlm[l2-m], -signs[l2+1-m], &sign_dPl2m) - logx2m1;

    /* Pl1m*Pl2m */
    res->lnPl1mPl2m    = lnPl1m + lnPl2m;
    res->sign_Pl1mPl2m = common_sign * sign_Pl1m * sign_Pl2m;

    /* Pl1m*dPl2m */
    res->lnPl1mdPl2m    = lnPl1m + lndPl2m;
    res->sign_Pl1mdPl2m = common_sign * sign_Pl1m * sign_dPl2m;

    /* dPl1m*dPl2m */
    res->lndPl1mPl2m    = lndPl1m + lnPl2m;
    res->sign_dPl1mPl2m = common_sign * sign_dPl1m * sign_Pl2m;

    /* dPl1m*dPl2m */
    res->lndPl1mdPl2m    = lndPl1m + lndPl2m;
    res->sign_dPl1mdPl2m = common_sign * sign_dPl1m * sign_dPl2m;
}

/*@}*/


/**
* @name Gaunt coefficients
*/
/*@{*/

/**
 * @brief Determine qmax
 *
 * @param [in]  n  \f$n\f$
 * @param [in]  nu \f$\nu\f$
 * @param [in]  m  \f$m=\mu\f$
 * @return a0 \f$q_\mathrm{max}\f$
 */
inline int gaunt_qmax(const int n, const int nu, const int m)
{
    int xi = (n+nu-2*m)/2;
    return MIN(MIN(n,nu), xi);
}

/**
 * @brief Calculate \f$\log a_0\f$
 *
 * Cf. eq. (20).
 *
 * @param [in]  n  \f$n\f$
 * @param [in]  nu \f$\nu\f$
 * @param [in]  m  \f$m=\mu\f$
 * @return a0 \f$\log a_0\f$
 */
inline float80 gaunt_log_a0(int n, int nu, int m)
{
    return lgamma80(2*n+1)-lgamma80(n+1)+lgamma80(2*nu+1)-lgamma80(1+nu)+lgamma80(n+nu+1)-lgamma80(2*n+2*nu+1)+lgamma80(1+n+nu-2*m)-lgamma80(1+n-m)-lgamma80(1+nu-m);
}

/**
 * @brief Calculate \f$a_0\f$
 *
 * Cf. eq. (20).
 *
 * @param [in]  n  \f$n\f$
 * @param [in]  nu \f$\nu\f$
 * @param [in]  m  \f$m=\mu\f$
 * @return a0 \f$a_0\f$
 */
inline float80 gaunt_a0(int n, int nu, int m)
{
    return exp80(gaunt_log_a0(n,nu,m));
}

/* eq. (3) */
#define alpha(p, n, nu) (((float80)(pow_2(p)-pow_2(n+nu+1))*(pow_2(p)-pow_2(n-nu)))/(4*pow_2(p)-1))


/**
 * @brief Calculate Gaunt coefficients
 *
 * Determine Gaunt coefficients \f$a(m, n, mu, nu, p)\f$ for \f$m\f$, \f$n\f$,
 * \f$\mu\f$ and \f$\nu\f$ fixed.  These coefficients can be used to express
 * the product of two associated Legendre polynomials:
 *
 * \f[
 * P_n^m(x) P_{\nu}^{\mu}(x) = a_0 \sum_{q=0}^{q_\mathrm{max}} \tilde a_q P_{n+\nu-2q}^{m+mu}(x)
 * \f]
 *
 * \f$q_\mathrm{max}\f$ is the upper bound of summation, \f$a_0\f$ is the
 * prefactor and \f$\tilde a_q\f$ are normalized Gaunt coefficients.
 *
 * See [1] for more information, especially chapter 3. There is a brief
 * outline how to calculate Gaunt coefficients at the end of the chapter.
 *
 * Ref.: [1] Y.-L. Xu, J. Comp. Appl. Math. 85, 53 (1997)
 *
 * @param [in]  n  \f$n\f$
 * @param [in]  nu \f$\nu\f$
 * @param [in]  m  \f$m=\mu\f$
 * @param [out] a_tilde \f$\tilde a_q\f$ list of normalized Gaunt coefficients
 */
void gaunt(const int n_, const int nu_, const int m_, float80 a_tilde[])
{
    const float80 n  = n_;
    const float80 nu = nu_;
    const float80 m  = m_;
    const float80 n4 = n+nu-2*m;

    /* eq. (24) */
    const int qmax = gaunt_qmax(n,nu,m);

    /* eq. (28) */
    const float80 Ap = -2*m*(n-nu)*(n+nu+1);

    if(qmax < 0)
        return;

    a_tilde[0] = 1;
    if(qmax == 0)
        return;

    /* eq. (29) */
    a_tilde[1] = (n+nu-1.5)*(1-(2*n+2*nu-1)/(n4*(n4-1))*((m-n)*(m-n+1)/(2*n-1)+(m-nu)*(m-nu+1)/(2*nu-1)));
    if(qmax == 1)
        return;

    /* eq. (35) */
    a_tilde[2] = (2*n+2*nu-1)*(2*n+2*nu-7)/4*( (2*n+2*nu-3)/(n4*(n4-1)) * ( (2*n+2*nu-5)/(2*(n4-2)*(n4-3)) \
                * ( (m-n)*(m-n+1)*(m-n+2)*(m-n+3)/(2*n-1)/(2*n-3) \
                + 2*(m-n)*(m-n+1)*(m-nu)*(m-nu+1)/((2*n-1)*(2*nu-1)) \
                + (m-nu)*(m-nu+1)*(m-nu+2)*(m-nu+3)/(2*nu-1)/(2*nu-3) ) - (m-n)*(m-n+1)/(2*n-1) \
                - (m-nu)*(m-nu+1)/(2*nu-1) ) +0.5);

    for(int q = 3; q <= qmax; q++)
    {
        float80 c0,c1,c2;
        const float80 p = n+nu-2*q;
        const float80 p1 = p-2*m;
        const float80 p2 = p+2*m;

        if(Ap != 0)
        {
            /* eqs. (26), (27) */
            c0 = (p+2)*(p+3)*(p1+1)*(p1+2)*Ap*alpha(p+1,n,nu);
            c1 = Ap*(Ap*Ap \
               + (p+1)*(p+3)*(p1+2)*(p2+2)*alpha(p+2,n,nu) \
               + (p+2)*(p+4)*(p1+3)*(p2+3)*alpha(p+3,n,nu));
            c2 = -(p+2)*(p+3)*(p2+3)*(p2+4)*Ap*alpha(p+4,n,nu);

            a_tilde[q] = (c1*a_tilde[q-1] + c2*a_tilde[q-2])/c0;
        }
        else
            /* eq. (30) */
            a_tilde[q] = (p+1)*(p2+2)*alpha(p+2,n,nu)*a_tilde[q-1] / ((p+2)*(p1+1)*alpha(p+1,n,nu));
    }
}

/*@}*/
