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
 * @param [in] n integer
 * @retval log log(n)
 */
double logi(unsigned int n)
{
    if(n < lookup_logi_elems)
        return lookup_logi[n];
    else
        return log(n);
}


/** @brief Calculate log(n!) = log(Γ(n+1))
 *
 * @param [in] n integer
 * @retval lfac log(n!) = log(Γ(n+1))
 */
double lfac(unsigned int n)
{
    if(n < lookup_lfac_elems)
        return lookup_lfac[n];
    else
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

double logadd_ms(log_t list[], const int len, sign_t *sign)
{
    double max = list[0].v;

    for(int i = 1; i < len; i++)
        if(list[i].v > max)
            max = list[i].v;

    double sum = list[0].s*exp(list[0].v-max);
    for(int i = 1; i < len; i++)
        sum += list[i].s*exp(list[i].v-max);

    *sign = SGN(sum);
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
                double lnKn_new = logadd(logi(2*l-1)+lnKnup-logx, lnKnu);
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
        return k*M_LOG2 + lfac(k);
    }
    else /* odd */
    {
        int k = (n+1)/2;
        return lfac(2*k) - k*M_LOG2 - lfac(k);
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
 * This function calculates associated Legendre functions for \f$m \ge 0\f$ and
 * \f$x \ge 1\f$.
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
 *     (+i)^2 = (-i)^2 = -1.
 * \f]
 *
 * So, we actually calculate
 * \f[
 *     P_\ell^m(x) = (x^2-1)^{m/2} \frac{\mathrm{d}^m}{\mathrm{d}x^m} P_\ell(x) .
 * \f]
 * Note that we don't include the factor i^m in our calculuation.
 *
 * The values of Plm are calculated from Plm(l=m,m=m,x) to Plm(l=lmax,m=m,x).
 * The associated Legendre polynomials are calculated using a recurrence
 * relation. Each Plm is normalized as Plm(l,m,x)/factor^l.
 *
 * This function calculates the associated Legendre polynomials for
 * \f$\ell=m,...,\ell_\mathrm{max}\f$.
 *
 * See https://en.wikipedia.org/wiki/Associated_Legendre_polynomials .
 *
 * @param [in] lmax maximum value of \f$\ell\f$
 * @param [in] m order
 * @param [in] x argument
 * @param [in] factor scaling factor
 * @param [out] array array of values
 */
void Plm_array(int lmax, int m, double x, double factor, double array[])
{
    array[0] = 1;

    if(lmax == m)
    {
        array[0] = exp(ln_factorial2(2*m-1) + m/2.*log((x+1)/factor*(x-1)/factor));
        for(int ll = 1; ll < lmax-m+1; ll++)
            array[ll] = 0;
        return;
    }
    else if(m > 0)
        //array[0] = exp(ln_factorial2(2*m-1) + m/2.*log((x+1)/factor*(x-1)/factor));
        factor /= exp((ln_factorial2(2*m-1) + m/2.*log((x+1)/factor*(x-1)/factor))/(lmax-m));

    array[1] = x*(2*m+1)*array[0]/factor;

    if(array[1] == 0)
    {
        for(int ll = 2; ll < lmax-m+1; ll++)
            array[ll] = 0;

        return;
    }

    const double y = factor*x;
    const double factor2 = pow_2(factor);

    if(isinf(factor2))
    {
        for(int ll = m+2; ll < lmax+1; ll++)
            array[ll-m] = 0;

        return;
    }

    for(int ll = m+2; ll < lmax+1; ll++)
        array[ll-m] = ((2*ll-1)*y*array[ll-m-1] - (ll+m-1)*array[ll-m-2])/((ll-m)*factor2);
}


double Plm(int l, int m, double x, double factor)
{
    double array[l-m+1];
    Plm_array(l, m, x, factor, array);
    return array[l-m];
}


/** @brief Estimate value of Plm(x) for x >> 1
 *
 * This function computes the value of Plm(x) using an approximation for large
 * arguments x >> 1 and returns the logarithm of the estimate.
 *
 * @param [in] l l
 * @param [in] m m
 * @param [in] x argument
 * @retval estimate, ≈log(Plm(x))
 */
double Plm_estimate(int l, int m, double x)
{
    return lfac(2*l)-l*M_LOG2-lfac(l)-lfac(l-m)+l*log(x);
}

/*@}*/
