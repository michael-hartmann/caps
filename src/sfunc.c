/**
 * @file   sfunc.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   October, 2016
 * @brief  various functions implementing mostly special functions
 */

#include <stdlib.h>
#include <math.h>

#include "lookup.h"
#include "sfunc.h"
#include "utils.h"


/** @brief Calculate log(x) for x integer
 *
 * This function uses a lookup table to avoid calling log() for n "small".
 *
 * @param [in] n double
 * @retval log log(n)
 */
double logi(unsigned int n)
{
    if(n < lookup_logi_elems)
        return lookup_logi[n];
    else
        return log(n);
}

double lfac(unsigned int n)
{
    return lgamma(1+n);
}

/**
 * @brief Sum elements in array
 *
 * This function calculates the sum of the elements of the array input. This
 * function uses the Kahan summation algorithm to reduce numerical error.
 *
 * The algorithm is taken from Wikipedia, see
 * https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 *
 * @param [in] input array
 * @param [in] N length of array
 * @return sum sum of array elements
 */
double kahan_sum(double input[], size_t N)
{
    double sum = 0;
    double c = 0; /* running compensation for lost low-order bits */

	for(size_t i = 0; i < N; i++)
    {
        double y = input[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}

/**
 * @brief Find maximum in double array
 *
 * @param [in] input array
 * @param [in] N length of array (must be > 0)
 * @return maximum maximum of array elements
 */
double max(double input[], size_t N)
{
    double m = input[0];

    for(size_t i = 1; i < N; i++)
        m = MAX(m,input[i]);

    return m;
}

/**
 * @brief Find minimum in double array
 *
 * @param [in] input array
 * @param [in] N length of array (must be > 0)
 * @return minimum minimum of array elements
 */
double min(double input[], size_t N)
{
    double m = input[0];

    for(size_t i = 1; i < N; i++)
        m = MIN(m,input[i]);

    return m;
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
double logadd(const double log_a, const double log_b)
{
    if(isinf(log_a) && log_a < 0)
        return log_b;
    else if(isinf(log_b) && log_b < 0)
        return log_a;

    if(log_a > log_b)
        return log_a + log1p(exp(log_b-log_a));
    else
        return log_b + log1p(exp(log_a-log_b));
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
double logadd_s(const double log_a, const sign_t sign_a, const double log_b, const sign_t sign_b, sign_t *sign)
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
        return log_a + log1p(sign_a*sign_b*exp(log_b-log_a));
    }
    else
    {
        *sign = sign_b;
        return log_b + log1p(sign_a*sign_b*exp(log_a-log_b));
    }
}

/*@}*/


/**
* @name Bessel functions
*/
/*@{*/

void bessel_lnInuKnu(int nu, const double x, double *lnInu_p, double *lnKnu_p)
{
    const double logx = log(x);
    double lnKnu = 0, lnKnup = log1p(1./x);

    /* calculate Knu, Knup */
    {
        const double prefactor = -x+0.5*(M_LOGPI-M_LOG2-logx);

        if(nu == 0)
        {
            lnKnu  = prefactor+lnKnu;
            lnKnup = prefactor+lnKnup;
        }
        else
        {
            for(int l = 2; l <= nu+1; l++)
            {
                double lnKn_new = logadd(log(2*l-1)+lnKnup-logx, lnKnu);
                lnKnu  = lnKnup;
                lnKnup = lnKn_new;
            }

            lnKnup = prefactor+lnKnup;
            lnKnu  = prefactor+lnKnu;
        }

        TERMINATE(!isfinite(lnKnup), "Couldn't calculate Bessel functions, nu=%d, x=%g\n", nu, x);

        if(lnKnu_p != NULL)
            *lnKnu_p = lnKnu;
    }

    if(lnInu_p != NULL)
    {
        #define an(n,nu,x) (2*((nu)+0.5+(n))/(x))

        double num   = an(2,nu,x)+1/an(1,nu,x);
        double denom = an(2,nu,x);
        double ratio = (an(1,nu,x)*num)/denom;
        double ratio_last = 0;

        for(int l = 3; 1; l++)
        {
            num   = an(l,nu,x)+1/num;
            denom = an(l,nu,x)+1/denom;
            ratio *= num/denom;

            if(ratio_last != 0 && fabs(1.-ratio/ratio_last) < 1e-20)
                break;

            ratio_last = ratio;
        }

        *lnInu_p = -logx-lnKnu-logadd(lnKnup-lnKnu, -log(ratio));
        #undef an
    }
}


double bessel_lnKnu(const int nu, const double x)
{
    double Knu;
    bessel_lnInuKnu(nu, x, NULL, &Knu);
    return Knu;
}


double bessel_lnInu(const int nu, const double x)
{
    double Inu;
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
double ln_factorial2(unsigned int n)
{
    /* see e.g. http://en.wikipedia.org/wiki/Double_factorial */
    if(n == 0 || n == 1) /* 0!! = 1!! = 0 */
        return 0;

    if(n % 2 == 0) /* even */
    {
        int k = n/2;
        return k*M_LOG2 + lgamma(1+k);
    }
    else /* odd */
    {
        int k = (n+1)/2;
        return lgamma(1+2*k) - k*M_LOG2 - lgamma(1+k);
    }
}

double factorial2(unsigned int n)
{
    if(n < lookup_factorial2_elems)
        return lookup_factorial2[n];
    else
        return exp(ln_factorial2(n));
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
 * See https://en.wikipedia.org/wiki/Associated_Legendre_polynomials .
 *
 * @param [in] lmax maximum value of \f$\ell\f$
 * @param [in] m order
 * @param [in] x position to evaluate associated Legendre polynomial
 * @param [out] lnplm array of logarithms of values
 * @param [out] sign corressponding signs
 */
void plm_lnPlm_array(int lmax, int m, double x, double lnplm[], sign_t sign[])
{
    double logx = log(x);

    if(m == 0)
    {
        sign[0] = +1;
        lnplm[0] = 0; // log(1)
    }
    else
    {
        sign[0]  = MPOW((int)(m/2) + m%2);
        lnplm[0] = ln_factorial2(2*m-1) + m*0.5*log(pow_2(x)-1); // l=m,m=m
    }

    if(lmax == m)
        return;

    sign[1]  = sign[0];
    lnplm[1] = lnplm[0]+logx+logi(2*m+1); // l=m+1, m=m

    /* (l-m)*P_l^m(x) = x*(2l-1)*P_(l-1)^m(x) - (l+m-1)*P_(l-2)^m(x) */
    for(int l = m+2; l <= lmax; l++)
        lnplm[l-m] = logadd_s(logi(2*l-1)+logx+lnplm[l-m-1], sign[l-m-1], logi(l+m-1)+lnplm[l-m-2], -sign[l-m-2], &sign[l-m]) - logi(l-m);
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
double plm_lnPlm(int l, int m, double x, sign_t *sign)
{
    double plm[l-m+1];
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
double plm_Plm(int l, int m, double x)
{
    sign_t sign;
    double value = plm_lnPlm(l, m, x, &sign);
    return sign*exp(value);
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
double plm_lndPlm(int l, int m, double x, sign_t *sign)
{
    if(m == 0 && l == 1)
        return 0;
    else
    {
        const int lmax = l+1;
        double plm[lmax-m+1];
        sign_t signs[lmax-m+1];

        plm_lnPlm_array(lmax, m, x, plm, signs);

        return logadd_s(logi(l-m+1)+plm[l+1-m], signs[l+1-m], logi(l+1)+log(x)+plm[l-m], -signs[l+1-m], sign) - log(pow_2(x)-1);
    }
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
double plm_dPlm(int l, int m, double x)
{
    sign_t sign = 0;
    double value = plm_lndPlm(l, m, x, &sign);
    return sign*exp(value);
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
void plm_PlmPlm(int l1, int l2, int m, double x, plm_combination_t *res)
{
    const int lmax = MAX(l1,l2)+1;
    double lnPlm[lmax-m+1];
    sign_t signs[lmax-m+1];

    plm_lnPlm_array(lmax, m, x, lnPlm, signs);

    plm_PlmPlm_from_array(l1, l2, m, x, lnPlm, signs, res);
}

void plm_PlmPlm_from_array(int l1, int l2, int m, double x, double lnPlm[], sign_t signs[], plm_combination_t *res)
{
    const double logx = log(x);
    const double logx2m1 = log(pow_2(x)-1);
    double lnPl1m, lnPl2m, lndPl1m, lndPl2m;
    sign_t sign_Pl1m, sign_Pl2m, sign_dPl1m, sign_dPl2m;
    sign_t common_sign = MPOW(m%2);

    lnPl1m    = lnPlm[l1-m];
    sign_Pl1m = signs[l1-m];
    lnPl2m    = lnPlm[l2-m];
    sign_Pl2m = signs[l2-m];

    /* Plm(l=1,m=0,x) = x */
    if(m == 0 && l1 == 1)
    {
        sign_dPl1m = +1;
        lndPl1m = 0;
    }
    else
        lndPl1m = logadd_s(logi(l1-m+1)+lnPlm[l1+1-m], signs[l1+1-m], logi(l1+1)+logx+lnPlm[l1-m], -signs[l1+1-m], &sign_dPl1m) - logx2m1;
    if(m == 0 && l2 == 1)
    {
        sign_dPl2m = +1;
        lndPl2m = 0;
    }
    else
        lndPl2m = logadd_s(logi(l2-m+1)+lnPlm[l2+1-m], signs[l2+1-m], logi(l2+1)+logx+lnPlm[l2-m], -signs[l2+1-m], &sign_dPl2m) - logx2m1;

    /* Pl1m*Pl2m */
    res->lnPl1mPl2m    = lnPl1m + lnPl2m;
    res->sign_Pl1mPl2m = common_sign * sign_Pl1m * sign_Pl2m;

    /* Pl1m*dPl2m */
    res->lnPl1mdPl2m    = lnPl1m + lndPl2m;
    res->sign_Pl1mdPl2m = common_sign * sign_Pl1m * sign_dPl2m;

    /* dPl1m*Pl2m */
    res->lndPl1mPl2m    = lndPl1m + lnPl2m;
    res->sign_dPl1mPl2m = common_sign * sign_dPl1m * sign_Pl2m;

    /* dPl1m*dPl2m */
    res->lndPl1mdPl2m    = lndPl1m + lndPl2m;
    res->sign_dPl1mdPl2m = common_sign * sign_dPl1m * sign_dPl2m;
}

/*@}*/
