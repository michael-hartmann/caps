/**
 * @file   sfunc.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   July, 2017
 * @brief  various functions implementing mostly special functions
 */

#include <stdlib.h>
#include <math.h>

#include "besselI.h"
#include "logfac.h"
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

/** @brief Compute sqrt(1+x)-1
 *
 * If x is small, sqrt(1+x)≈1 and a loss of significance occures when
 * calculating sqrt(1+x)-1.
 *
 * For this reason we use the Taylor series
 *     sqrt(1+x) ≈ 1 + x/2 - x²/8 + x³/16 - 5*x^4/128 + ...
 * to avoid a loss of significance if x is small.
 *
 * @param [in] x
 * @retval sqrt(1+x)-1
 */
double sqrtpm1(double x)
{
    if(fabs(x) < 1e-5)
        /* sqrt(1+x)-1 = x/2 - x²/8 + x³/16 - 5*x^4/128 + O(x^5) */
        return x/2 - pow_2(x)/8 + pow_3(x)/16 - 5*pow_4(x)/128;
    else
        return sqrt(1+x)-1;
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

/* @brief Calculate I_{nu+1/2}/I_{nu+3/2}
 *
 * @param nu order
 * @param x argument
 * @retval ratio
 */
double bessel_continued_fraction(int nu, double x)
{
    /* it's faster to calculate the inverse of x only once */
    const double invx = 1/x;

    const double a1 = (2*(nu+1)+1)*invx;
    const double a2 = (2*(nu+2)+1)*invx;

    double num   = a2+1/a1;
    double denom = a2;
    double ratio = a1*num/denom;
    double ratio_last = 0;

    for(int l = 3; 1; l++)
    {
        const double an = (2*nu+1+2*l)*invx;
        num   = an+1/num;
        denom = an+1/denom;
        ratio *= num/denom;

        if(ratio == ratio_last)
            return ratio;

        ratio_last = ratio;
    }
}

double bessel_lnInu(int nu, double x)
{
    double lnInu;
    bessel_lnInuKnu(nu, x, &lnInu, NULL);
    return lnInu;
}

double bessel_lnKnu(int nu, double x)
{
    double lnKnu;
    bessel_lnInuKnu(nu, x, NULL, &lnKnu);
    return lnKnu;
}

void bessel_lnInuKnu(int nu, const double x, double *lnInu_p, double *lnKnu_p)
{
    const double logx = log(x);
    const double invx = 1/x;

    double lnKnu, lnKnup;
    double Knu = 1, Knup = 1+invx;
    double prefactor = -x+0.5*(M_LOGPI-M_LOG2-logx);

    /* calculate Knu, Knup */
    if(nu == 0)
    {
        lnKnu  = prefactor+log(Knu);
        lnKnup = prefactor+log(Knup);
    }
    else
    {
        for(int l = 2; l <= nu+1; l++)
        {
            double Kn_new = (2*l-1)*Knup*invx + Knu;
            Knu  = Knup;
            Knup = Kn_new;

            if(Knu > 1e100)
            {
                Knu  *= 1e-100;
                Knup *= 1e-100;
                prefactor += log(1e100);
            }
        }

        lnKnup = prefactor+log(Knup);
        lnKnu  = prefactor+log(Knu);
    }

    TERMINATE(!isfinite(lnKnup), "Couldn't calculate Bessel functions, nu=%d, x=%g\n", nu, x);

    if(lnKnu_p != NULL)
        *lnKnu_p = lnKnu;

    if(lnInu_p != NULL)
    {
        double ratio = bessel_continued_fraction(nu,x);
        *lnInu_p = -logx-logadd(lnKnup, lnKnu-log(ratio));
    }
}

/*@}*/


/**
* @name Associated Legendre polynomials
*/
/*@{*/

/**
 * @brief Associated Legendre polynomials for argument x > 1
 *
 * This function calculates associated Legendre functions for m >= 0 and x > 0.
 *
 * Associated Legendre polynomials are defined as follows:
 *     Plm(x) = (-1)^m (1-x^2)^(m/2) D^m/Dx^m  Pl(x)
 * where Pl(x) denotes a Legendre polynomial.
 *
 * As Pl(x) are ordinary polynomials, the only problem is the term
 * (1-x^2)^(m/2) when extending the domain to values of x > 1. We will use the
 * convention sqrt(-x) = +i sqrt(x).
 *
 * Note: Products of associated legendre polynomials with common m are
 * unambiguous, because (+i)^2 = (-i)^2 = -1.
 *
 * What this functions actually computes:
 *     Plm(x) = (x^2-1)^(m/2) D^m/Dx^m  Pl(x) .
 * Note that we don't include the factor i^m in our calculuation.
 *
 * See also https://en.wikipedia.org/wiki/Associated_Legendre_polynomials .
 *
 * @param [in] l degree
 * @param [in] m order
 * @param [in] x argument
 * @retval log(Plm(x))
 */
double Plm(int l, int m, double x)
{
    if((l-m) <= 200)
        return Plm_upwards(l, m, x);
    else
        return Plm_downwards(l, m, x);
}

/* @brief Associated Legendre polynomials using upwards recurrence relation
 *
 * The values of Plm are calculated from Plm(l=m,m=m,x) to Plm(l=lmax,m=m,x).
 * The associated Legendre polynomials are calculated using a recurrence
 * relation http://dlmf.nist.gov/14.10.E3 .
 *
 * @param [in] l degree
 * @param [in] m order
 * @param [in] x argument
 */
double Plm_upwards(int l, int m, double x)
{
    double array[l-m+1];
    double log_prefactor = ln_factorial2(2*m-1) + m/2.*log((x+1)*(x-1));

    if(l == m)
        return log_prefactor;

    array[0] = 1;
    array[1] = x*(2*m+1)*array[0];

    if(array[1] == 0)
        return -INFINITY;

    for(int ll = 2; ll < l+1-m; ll++)
    {
        const double k = (2.*m-1.)/ll;

        array[ll] = (2+k)*x*array[ll-1] - (1+k)*array[ll-2];

        const double elem = fabs(array[ll]);
        if(elem < 1e-100)
        {
            log_prefactor -= log(1e100);
            array[ll]   *= 1e100;
            array[ll-1] *= 1e100;
        }
        else if(elem > 1e+100)
        {
            log_prefactor += log(1e100);
            array[ll]   *= 1e-100;
            array[ll-1] *= 1e-100;
        }
    }

    if(isnan(array[l-m]))
        return NAN;

    return log_prefactor+log(fabs(array[l-m]));
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

/* Legendre polynomial Pl
 *
 * Evaluation of Pl(x) for x>=1 using an asymptotic expansion provided that
 *      (l+1)*sqrt((x+1)*(x-1)) >= 25.
 *
 * O(1) computation of Legendre polynomials and Gauss-Legendre nodes and
 * weights for parallel computing, section 3.2.
 */
static double _Pl1(int l, double x, double sinhxi)
{
    const double xi = acosh(x);
    /* exp(xi) = exp(acosh(x)) = x+ sqrt((x+1)*(x-1)) */
    const double expxi = x+sinhxi;

    /* Cl0 */
    double Clm = exp( lfac(l)-lfac(2*l+2)-M_LOGPI/2+(l+1)*log(4)+lfac(l+1) );

    double sum = 0;
    double sinhxi_m = 1; /* sinh(xi)**m */
    double expxi_m  = 1; /* exp(xi)**m */
    for(int m = 0; m < 17; m++)
    {
        sum += Clm*(expxi_m+exp(-(m+2*l+1)*xi))/sinhxi_m;

        sinhxi_m *= sinhxi;
        expxi_m  *= expxi;
        Clm *= pow_2(m+0.5)/((2*m+2)*(l+1.5+m));
    }
    sum += Clm*(expxi_m+exp(-(2*l+18)*xi))/sinhxi_m; /* m=17 */

    return -(M_LOG2+M_LOGPI)/2 - log(sinhxi)/2 + log(sum) + (l+0.5)*xi;
}

/* equations (3.27)-(3.31) */
static double _fn(int n, double hn[13])
{
    switch(n)
    {
        case 0:
            return hn[0];
        case 2:
            return 1./8*hn[1] - 1./12*hn[2];
        case 4:
            return 11./384*hn[2] - 7./160*hn[3] + 1./160*hn[4];
        case 6:
            return 173./15360*hn[3] - 101./3584*hn[4] + 671./80640*hn[5] - 61./120960*hn[6];
        case 8:
            return 22931./3440640*hn[4] - 90497./3870720*hn[5] + 217./20480*hn[6] - 1261./967680*hn[7] + 1261./29030400*hn[8];
        case 10:
            return 1319183./247726080*hn[5] - 10918993./454164480*hn[6] + 1676287./113541120*hn[7] - 7034857./2554675200*hn[8] + 1501./8110080*hn[9] - 79./20275200*hn[10];
        case 12:
            return 233526463./43599790080*hn[6] - 1396004969./47233105920*hn[7] + 2323237523./101213798400*hn[8] - 72836747./12651724800*hn[9] + 3135577./5367398400*hn[10] - 1532789./61993451520*hn[11] + 66643./185980354560*hn[12];
        default:
            return 0;
    }
}


/* Legendre polynomial Pl
 *
 * Evaluation of Pl(x) for x>=1 using an asymptotic expansion provided that
 *      (l+1)*sqrt((x+1)*(x-1)) < 25.
 *
 * O(1) computation of Legendre polynomials and Gauss-Legendre nodes and
 * weights for parallel computing, section 3.3.
 */
static double _Pl2(int l, double x)
{
    const int N = 14;
    const double v = l+0.5;
    const double xi = acosh(x);

    const double y = xi*v;
    double hn[13];

    double yn = 1; /* y^n */
    for(int n = 0; n < 13; n++)
    {
        hn[n] = yn*besselI(n,-y);
        yn *= y;
    }

    double s = 0;
    double vn = 1; /* v**n */
    for(int n = 0; n < N; n += 2)
    {
        s += _fn(n,hn)/vn;
        vn *= v*v;
    }

    return log(s);
}


/* Legendre polynomial Pl
 *
 * Evaluation of Pl(x) for x>=1 using the recurrence relation
 *      (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x).
 */
static double _Pl3(int l, double x)
{
    double log_prefactor = 0;

    if(l == 0)
        return 0;
    if(l == 1)
        return log(x);

    double v0 = 1, v1 = x;
    for(int n = 1; n < l; n++)
    {
        double v = ((2*n+1)*x*v1-n*v0)/(n+1);
        v0 = v1;
        v1 = v;

        if(v > 1e100)
        {
            log_prefactor += log(1e100);
            v0 *= 1e-100;
            v1 *= 1e-100;
        }
    }

    return log_prefactor+log(v1);
}


/* @brief Legendre polynomial Pl
 *
 * Evaluation of Pl(x) for x>=1.
 *
 * For l < 100 a recurrence relation is used (see _Pl3), otherwise asymptotic
 * expansions are used (see _Pl1 and _Pl2).
 *
 * Function returns log(Pl(x)).
 *
 * @param [in] l degree
 * @param [in] x argument
 * @retval log(Pl(x))
 */
double Pl(int l, double x)
{
    if(l < 100)
        return _Pl3(l,x);

    const double sinhxi = sqrt((x+1)*(x-1));
    if((l+1)*sinhxi > 25)
        /* (l+1)*sinh(xi) >= 25. */
        return _Pl1(l,x,sinhxi);
    else
        /* (l+1)*sinh(xi) < 25. */
        return _Pl2(l,x);
}

/* Calculate fraction P_l^m/P_l^{m-1}
 *
 * The fraction is computed using a continued fraction, see http://dlmf.nist.gov/14.14.E1 .
 *
 * To evaluate the continued fraction, we use http://dlmf.nist.gov/1.12#E5 and
 * http://dlmf.nist.gov/1.12#E6 .
 * See also Numerical Recipes in C, ch. 5.2, Evaluation of Continued Fractions
 *
 * @param [in] l degree
 * @param [in] m order
 * @param [in] x argument
 * @retval Plm(l,m,x)/Plm(l,m-1,x)
 */
static double _cf(const int l, const int m, const double x)
{
    const double alpha = (1-1/(x*x))/4;

    /* initial values for recrrence */
    double A,B, Amm = 0, Bmm = 1;
    double Am = x*(l-m+1)*(l+m)*alpha; /* a0 */
    double Bm = m;                     /* b0 */

    int j = 1;
    for(int i = 1; i < 2048; i++)
    {
        /* do 8 iterations */
        for(int k = 0; k < 8; k++)
        {
            const double aj = (l-m+1-j)*alpha*(l+m+j);
            const double bj = m+j;

            A = bj*Am + aj*Amm;
            B = bj*Bm + aj*Bmm;

            Amm = Am;
            Am  = A;
            Bmm = Bm;
            Bm  = B;

            j++;
        }

        /* compute last and current continued fraction */
        double f1 = Amm/Bmm, f2 = Am/Bm;

        /* if they are identical, we're done */
        if(f1 == f2)
            return 2*f1/sqrt((x+1)*(x-1));

        /* rescale */
        Amm /= Bm;
        Am  /= Bm;
        Bmm /= Bm;
        Bm   = 1;
    }

    TERMINATE(true, "l=%d, m=%d, x=%.15g", l,m,x);

    return NAN;
}


/* @brief Associated Legendre polynomials using downwards recurrence relation
 *
 * First, the fraction Plm(l,m,x)/Plm(l,m-1,x) is computed using a continued
 * fraction, see _cf. Then the downwards recurrence relation
 * http://dlmf.nist.gov/14.10.E1 is used from Plm(l,m,x) to Plm(l,0,x).
 * Together with Pl(x) one has the solution.
 *
 * @param [in] l degree
 * @param [in] m order
 * @param [in] x argument
 */
double Plm_downwards(int l, int m, double x)
{
    double vp,vm,c;
    double prefactor = Pl(l,x); /* value of Pl(x) */

    if(m == 0)
        return prefactor;

    c = 2*x/sqrt((x-1)*(x+1));

    vp = +1;
    vm = -1/_cf(l,m,x);

    for(int mm = m-2; mm >= 0; mm--)
    {
        /* prevent integer overflows in denominator */
        double v = (vp-(mm+1)*vm*c)/((double)(l+1+mm)*(l-mm));
        vp = vm;
        vm = v;

        if(fabs(vm) < 1e-100)
        {
            vm *= 1e100;
            vp *= 1e100;
            prefactor += log(1e100);
        }
    }

    return prefactor-log(fabs(vm));
}

/*@}*/
