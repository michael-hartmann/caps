/**
 * @file   sfunc.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   April, 2015
 * @brief  various functions implementing mostly special functions
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "edouble.h"
#include "sfunc.h"
#include "utils.h"

/**
 * @brief Multiply two polynomials
 *
 * Multiply the coefficients of the polynomials p1 and p2 with and store the
 * result in p. p must have length of at least len_p1+len_p2-1.
 *
 * The polynomials are stored by coefficients \f$a_0,...,a_N\f$.
 *
 * This function uses the naive algorithm to calculate the product of two polynomials and thus
 * this function has complexity \f$\mathcal{O}(\mathrm{len\_p1}\cdot\mathrm{len\_p2})\f$.
 *
 * This function is thread-safe.
 *
 * @param [in] p1 polynomial
 * @param [in] len_p1 length of array p1
 * @param [in] p2 polynomial
 * @param [in] len_p2 length of array p2
 * @param [in] p polynomial \f$p=p1 \dot p2\f$
 */
void polymult(edouble p1[], int len_p1, edouble p2[], int len_p2, edouble p[])
{
    for(int i = 0; i < len_p1+len_p2-1; i++)
        p[i] = 0;

    for(int i = 0; i < len_p1; i++)
        for(int j = 0; j < len_p2; j++)
            p[i+j] += p1[i]*p2[j];
}


/**
* @name Add numbers given as logarithms
*/
/*@{*/


/**
 * @brief Add two numbers given by their logarithms.
 *
 * Both numbers are assumed to be nonnegative.
 *
 * @param [in] log_a number
 * @param [in] log_b number
 * @return log_sum \f$\log{\left[\exp{(\mathrm{log\_a})}+\exp{(log\_b)}\right]}\f$
 */
inline edouble logadd(const edouble log_a, const edouble log_b)
{
    if(isinf(log_a) && log_a < 0)
        return log_b;
    else if(isinf(log_b) && log_b < 0)
        return log_a;

    if(log_a > log_b)
        return log_a + log1pe(expe(log_b-log_a));
    else
        return log_b + log1pe(expe(log_a-log_b));
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
edouble logadd_s(const edouble log_a, const sign_t sign_a, const edouble log_b, const sign_t sign_b, sign_t *sign)
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
        return log_a + log1pe(sign_a*sign_b*expe(log_b-log_a));
    }
    else
    {
        *sign = sign_b;
        return log_b + log1pe(sign_a*sign_b*expe(log_a-log_b));
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
inline edouble logadd_m(const edouble list[], const int len)
{
    edouble sum;
    edouble max = list[0];

    for(int i = 1; i < len; i++)
        if(list[i] > max)
            max = list[i];

    sum = expe(list[0]-max);
    for(int i = 1; i < len; i++)
        sum += expe(list[i]-max);

    return max + loge(fabse(sum));
}


inline edouble logadd_ms(const edouble list[], const sign_t signs[], const int len, sign_t *sign)
{
    edouble sum;
    edouble max = list[0];

    for(int i = 1; i < len; i++)
        if(list[i] > max)
            max = list[i];

    sum = signs[0]*expe(list[0]-max);
    for(int i = 1; i < len; i++)
        sum += signs[i]*expe(list[i]-max);

    *sign = copysigne(1, sum);
    return max + loge(fabse(sum));
}

/*@}*/

/**
 * @brief Calculate logarithm of Binomial coefficient
 *
 * @param n
 * @param k
 * @return binomial \f$\log \left(\begin{array}{c} n \\ k \end{array}\right)\f$
 */
inline edouble lbinom(int n, int k)
{
    return lgammae(1+n)-lgammae(1+k)-lgammae(1+n-k);
}


/**
* @name Bessel functions
*/
/*@{*/

void bessel_lnInuKnu(int nu, const edouble x, edouble *lnInu_p, edouble *lnKnu_p)
{
    int l;
    edouble lnKnu = 0, lnKnup = log1pe(1./x);
    edouble logx = loge(x);

    // calculate Knu, Knup
    {
        edouble prefactor = -x+0.5*(LOGPI-LOG2-logx);

        if(nu == 0)
        {
            lnKnu  = prefactor+lnKnu;
            lnKnup = prefactor+lnKnup;
        }
        else
        {
            for(l = 2; l <= nu+1; l++)
            {
                edouble lnKn_new = logadd(loge(2*l-1)+lnKnup-logx, lnKnu);
                lnKnu  = lnKnup;
                lnKnup = lnKn_new;
            }

            lnKnup = prefactor+lnKnup;
            lnKnu  = prefactor+lnKnu;
        }

        if(isnan(lnKnup) || isinf(lnKnup))
        {
            WARN(1, "Couldn't calculate Bessel functions, falling back to approximations, nu=%d, x=%Lg\n", nu, x);

            /* so, we couldn't calculate lnKnup and lnKnu. Maybe we can at
             * least use the asymptotic behaviour for small values.
             */
            if(x < sqrt(nu)*1e3)
            {
                /* small arguments */
                lnKnu  = lgammae(nu+0.5)-LOG2+(nu+0.5)*(LOG2-logx);
                lnKnup = lgammae(nu+1.5)-LOG2+(nu+1.5)*(LOG2-logx);
            }
            else
            {
                TERMINATE(1, "Couldn't calculate Bessel functions, nu=%d, x=%Lg\n", nu, x);
                lnKnu = lnKnup = -INFINITY;
            }
        }

        if(lnKnu_p != NULL)
            *lnKnu_p = lnKnu;
    }

    if(lnInu_p != NULL)
    {
        #define an(n,nu,x) (2*((nu)+0.5+(n))/(x))

        edouble num   = an(2,nu,x)+1/an(1,nu,x);
        edouble denom = an(2,nu,x);
        edouble ratio = (an(1,nu,x)*num)/denom;
        edouble ratio_last = 0;

        l = 3;
        while(1)
        {
            num   = an(l,nu,x)+1/num;
            denom = an(l,nu,x)+1/denom;
            ratio *= num/denom;

            if(ratio_last != 0 && fabse(1.L-ratio/ratio_last) < 1e-20)
                break;

            ratio_last = ratio;
            l++;
        }

        *lnInu_p = -logx-lnKnu-logadd(lnKnup-lnKnu, -loge(ratio));
        #undef an
    }
}


edouble bessel_lnKnu(const int nu, const edouble x)
{
    edouble Knu;
    bessel_lnInuKnu(nu, x, NULL, &Knu);
    return Knu;
}


edouble bessel_lnInu(const int nu, const edouble x)
{
    edouble Inu;
    bessel_lnInuKnu(nu, x, &Inu, NULL);
    return Inu;
}

/*@}*/

double linspace(double start, double stop, int N, int i)
{
    if(N == 1)
        return start;
    else
        return start+(stop-start)*i/(N-1);
}


double logspace(double start, double stop, int N, int i)
{
    if(N == 1)
        return start;
    else
        return start*pow(pow(stop/start, 1./(N-1)), i);
}

/**
 * @brief Calculate double factorial \f$n!!\f$
 *
 * @param n non-negative integer
 * @return doublefactorial \f$n!!\f$
 */
edouble ln_doublefact(int n)
{
    if(n < 0)
        return NAN;

    /* see e.g. http://en.wikipedia.org/wiki/Double_factorial */
    if(n == 0 || n == 1) /* 0!! = 1!! = 0 */
        return 0;

    if(n % 2 == 0) /* even */
    {
        int k = n/2;
        return k*LOG2 + lgammae(1+k);
    }
    else /* odd */
    {
        int k = (n+1)/2;
        return lgammae(1+2*k) - k*LOG2 - lgammae(1+k);
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
inline void plm_lnPlm_array(int lmax, int m, edouble x, edouble lnplm[], sign_t sign[])
{
    edouble logx = loge(x);

    if(m == 0)
    {
        sign[0] = +1;
        lnplm[0] = 0; // log(1)
    }
    else
    {
        sign[0]  = MPOW((int)(m/2) + m%2);
        lnplm[0] = ln_doublefact(2*m-1) + m*0.5*loge(pow_2(x)-1); // l=m,m=m
    }

    if(lmax == m)
        return;

    sign[1]  = sign[0];
    lnplm[1] = lnplm[0]+logx+loge(2*m+1); // l=m+1, m=m

    for(int l = m+2; l <= lmax; l++)
    {
        lnplm[l-m] = logadd_s(loge(2*l-1)+logx+lnplm[l-m-1], sign[l-m-1], loge(l+m-1)+lnplm[l-m-2], -sign[l-m-2], &sign[l-m]);
        lnplm[l-m]-= loge(l-m);
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
edouble plm_lnPlm(int l, int m, edouble x, sign_t *sign)
{
    edouble plm[l-m+1];
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
edouble plm_Plm(int l, int m, edouble x)
{
    sign_t sign;
    edouble value = plm_lnPlm(l, m, x, &sign);
    return sign*expe(value);
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
edouble plm_lndPlm(int l, int m, edouble x, sign_t *sign)
{
    const int lmax = l+1;
    edouble plm[lmax-m+1];
    sign_t signs[lmax-m+1];

    plm_lnPlm_array(lmax, m, x, plm, signs);

    return logadd_s(loge(l-m+1)+plm[l+1-m], signs[l+1-m], loge(l+1)+loge(x)+plm[l-m], -signs[l+1-m], sign) - loge(pow_2(x)-1);
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
edouble plm_dPlm(int l, int m, edouble x)
{
    sign_t sign;
    edouble value = plm_lndPlm(l, m, x, &sign);
    return sign*expe(value);
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
void plm_PlmPlm(int l1, int l2, int m, edouble x, plm_combination_t *res)
{
    const int lmax = MAX(l1,l2)+1;
    edouble lnPlm[lmax-m+1];
    sign_t signs[lmax-m+1];
    edouble logx = loge(x);
    edouble logx2m1 = loge(pow_2(x)-1);
    edouble lnPl1m, lnPl2m, lndPl1m, lndPl2m;
    sign_t sign_Pl1m, sign_Pl2m, sign_dPl1m, sign_dPl2m;
    sign_t common_sign = MPOW(m%2);

    plm_lnPlm_array(lmax, m, x, lnPlm, signs);

    lnPl1m    = lnPlm[l1-m];
    sign_Pl1m = signs[l1-m];
    lnPl2m    = lnPlm[l2-m];
    sign_Pl2m = signs[l2-m];

    lndPl1m = logadd_s(loge(l1-m+1)+lnPlm[l1+1-m], signs[l1+1-m], loge(l1+1)+logx+lnPlm[l1-m], -signs[l1+1-m], &sign_dPl1m) - logx2m1;
    lndPl2m = logadd_s(loge(l2-m+1)+lnPlm[l2+1-m], signs[l2+1-m], loge(l2+1)+logx+lnPlm[l2-m], -signs[l2+1-m], &sign_dPl2m) - logx2m1;

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
inline edouble gaunt_log_a0(int n, int nu, int m)
{
    return lgamma(2*n+1)-lgamma(n+1)+lgamma(2*nu+1)-lgamma(1+nu)+lgamma(n+nu+1)-lgamma(2*n+2*nu+1)+lgamma(1+n+nu-2*m)-lgamma(1+n-m)-lgamma(1+nu-m);
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
inline edouble gaunt_a0(int n, int nu, int m)
{
    return expe(gaunt_log_a0(n,nu,m));
}

/* eq. (3) */
#define alpha(p, n, nu) (((edouble)(pow_2(p)-pow_2(n+nu+1))*(pow_2(p)-pow_2(n-nu)))/(4*pow_2(p)-1))


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
void gaunt(const int n_, const int nu_, const int m_, edouble a_tilde[])
{
    const edouble n  = n_;
    const edouble nu = nu_;
    const edouble m  = m_;
    const edouble n4 = n+nu-2*m;

    /* eq. (24) */
    const int qmax = gaunt_qmax(n,nu,m);

    /* eq. (28) */
    const edouble Ap = -2*m*(n-nu)*(n+nu+1);

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
        edouble c0,c1,c2;
        const edouble p = n+nu-2*q;
        const edouble p1 = p-2*m;
        const edouble p2 = p+2*m;

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
