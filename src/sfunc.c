/**
 * @file   sfunc.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   September, 2016
 * @brief  various functions implementing mostly special functions
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "sfunc.h"
#include "utils.h"


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
inline double logadd(const double log_a, const double log_b)
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


/**
 * @brief Add numbers given by their logarithms.
 *
 * len numbers in list will be added. The numbers are assumed to be positive.
 *
 * @param [in]  list array of numbers given by logarithm
 * @param [in]  len length of list
 * @return log_sum \f$\log{\sum_{i=1}^\mathrm{len} \mathrm{sign\_i}\cdot\exp{(\mathrm{log\_i})}}\f$
 */
inline double logadd_m(const double list[], const int len)
{
    double max = list[0];

    for(int i = 1; i < len; i++)
        if(list[i] > max)
            max = list[i];

    double sum = exp(list[0]-max);
    for(int i = 1; i < len; i++)
        sum += exp(list[i]-max);

    return max + log(sum);
}


inline double logadd_ms(const double list[], const sign_t signs[], const int len, sign_t *sign)
{
    double sum;
    double max = list[0];

    for(int i = 1; i < len; i++)
        if(list[i] > max)
            max = list[i];

    /* return +0 */
    if(max == -INFINITY)
    {
        *sign = +1;
        return -INFINITY;
    }

    sum = signs[0]*exp(list[0]-max);
    for(int i = 1; i < len; i++)
        sum += signs[i]*exp(list[i]-max);

    *sign = copysign(1, sum);
    return max + log(fabs(sum));
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
double ln_doublefact(int n)
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
        lnplm[0] = ln_doublefact(2*m-1) + m*0.5*log(pow_2(x)-1); // l=m,m=m
    }

    if(lmax == m)
        return;

    sign[1]  = sign[0];
    lnplm[1] = lnplm[0]+logx+log(2*m+1); // l=m+1, m=m

    /* (l-m)*P_l^m(x) = x*(2l-1)*P_(l-1)^m(x) - (l+m-1)*P_(l-2)^m(x) */
    for(int l = m+2; l <= lmax; l++)
        lnplm[l-m] = logadd_s(log(2*l-1)+logx+lnplm[l-m-1], sign[l-m-1], log(l+m-1)+lnplm[l-m-2], -sign[l-m-2], &sign[l-m]) - log(l-m);
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

        return logadd_s(log(l-m+1)+plm[l+1-m], signs[l+1-m], log(l+1)+log(x)+plm[l-m], -signs[l+1-m], sign) - log(pow_2(x)-1);
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
    double logx = log(x);
    double logx2m1 = log(pow_2(x)-1);
    double lnPl1m, lnPl2m, lndPl1m, lndPl2m;
    sign_t sign_Pl1m, sign_Pl2m, sign_dPl1m, sign_dPl2m;
    sign_t common_sign = MPOW(m%2);

    plm_lnPlm_array(lmax, m, x, lnPlm, signs);

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
        lndPl1m = logadd_s(log(l1-m+1)+lnPlm[l1+1-m], signs[l1+1-m], log(l1+1)+logx+lnPlm[l1-m], -signs[l1+1-m], &sign_dPl1m) - logx2m1;
    if(m == 0 && l2 == 1)
    {
        sign_dPl2m = +1;
        lndPl2m = 0;
    }
    else
        lndPl2m = logadd_s(log(l2-m+1)+lnPlm[l2+1-m], signs[l2+1-m], log(l2+1)+logx+lnPlm[l2-m], -signs[l2+1-m], &sign_dPl2m) - logx2m1;

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
inline double gaunt_log_a0(int n, int nu, int m)
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
inline double gaunt_a0(int n, int nu, int m)
{
    return exp(gaunt_log_a0(n,nu,m));
}

/* eq. (3) */
#define alpha(p, n, nu) (((pow_2(p)-pow_2(n+nu+1))*(pow_2(p)-pow_2(n-nu)))/(4*pow_2(p)-1))


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
void gaunt(const int n_, const int nu_, const int m_, double a_tilde[])
{
    const double n  = n_;
    const double nu = nu_;
    const double m  = m_;
    const double n4 = n+nu-2*m;

    /* eq. (24) */
    const int qmax = gaunt_qmax(n,nu,m);

    /* eq. (28) */
    const double Ap = -2*m*(n-nu)*(n+nu+1);

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
        double c0,c1,c2;
        const double p = n+nu-2*q;
        const double p1 = p-2*m;
        const double p2 = p+2*m;

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
