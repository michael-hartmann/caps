/**
 * @file   bessel.c
 * @author Stephen L. Moshier, Cephes Math Library Release 2.8, June 2000
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   July, 2017
 * @brief  Computation of Bessel functions
 */

#include <stdlib.h>
#include <math.h>

#include "sfunc.h"
#include "constants.h"
#include "bessel.h"

/* Chebyshev coefficients for exp(-x) I0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I0(x) } = 1.
 */
static double A0[] =
{
    -4.41534164647933937950E-18,
     3.33079451882223809783E-17,
    -2.43127984654795469359E-16,
     1.71539128555513303061E-15,
    -1.16853328779934516808E-14,
     7.67618549860493561688E-14,
    -4.85644678311192946090E-13,
     2.95505266312963983461E-12,
    -1.72682629144155570723E-11,
     9.67580903537323691224E-11,
    -5.18979560163526290666E-10,
     2.65982372468238665035E-9,
    -1.30002500998624804212E-8,
     6.04699502254191894932E-8,
    -2.67079385394061173391E-7,
     1.11738753912010371815E-6,
    -4.41673835845875056359E-6,
     1.64484480707288970893E-5,
    -5.75419501008210370398E-5,
     1.88502885095841655729E-4,
    -5.76375574538582365885E-4,
     1.63947561694133579842E-3,
    -4.32430999505057594430E-3,
     1.05464603945949983183E-2,
    -2.37374148058994688156E-2,
     4.93052842396707084878E-2,
    -9.49010970480476444210E-2,
     1.71620901522208775349E-1,
    -3.04682672343198398683E-1,
     6.76795274409476084995E-1
};

/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */
static double B0[] =
{
    -7.23318048787475395456E-18,
    -4.83050448594418207126E-18,
     4.46562142029675999901E-17,
     3.46122286769746109310E-17,
    -2.82762398051658348494E-16,
    -3.42548561967721913462E-16,
     1.77256013305652638360E-15,
     3.81168066935262242075E-15,
    -9.55484669882830764870E-15,
    -4.15056934728722208663E-14,
     1.54008621752140982691E-14,
     3.85277838274214270114E-13,
     7.18012445138366623367E-13,
    -1.79417853150680611778E-12,
    -1.32158118404477131188E-11,
    -3.14991652796324136454E-11,
     1.18891471078464383424E-11,
     4.94060238822496958910E-10,
     3.39623202570838634515E-9,
     2.26666899049817806459E-8,
     2.04891858946906374183E-7,
     2.89137052083475648297E-6,
     6.88975834691682398426E-5,
     3.36911647825569408990E-3,
     8.04490411014108831608E-1
};

/* Chebyshev coefficients for exp(-x) I1(x) / x
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
 */
static double A1[] =
{
     2.77791411276104639959E-18,
    -2.11142121435816608115E-17,
     1.55363195773620046921E-16,
    -1.10559694773538630805E-15,
     7.60068429473540693410E-15,
    -5.04218550472791168711E-14,
     3.22379336594557470981E-13,
    -1.98397439776494371520E-12,
     1.17361862988909016308E-11,
    -6.66348972350202774223E-11,
     3.62559028155211703701E-10,
    -1.88724975172282928790E-9,
     9.38153738649577178388E-9,
    -4.44505912879632808065E-8,
     2.00329475355213526229E-7,
    -8.56872026469545474066E-7,
     3.47025130813767847674E-6,
    -1.32731636560394358279E-5,
     4.78156510755005422638E-5,
    -1.61760815825896745588E-4,
     5.12285956168575772895E-4,
    -1.51357245063125314899E-3,
     4.15642294431288815669E-3,
    -1.05640848946261981558E-2,
     2.47264490306265168283E-2,
    -5.29459812080949914269E-2,
     1.02643658689847095384E-1,
    -1.76416518357834055153E-1,
     2.52587186443633654823E-1
};

/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
 */
static double B1[] =
{
     7.51729631084210481353E-18,
     4.41434832307170791151E-18,
    -4.65030536848935832153E-17,
    -3.20952592199342395980E-17,
     2.96262899764595013876E-16,
     3.30820231092092828324E-16,
    -1.88035477551078244854E-15,
    -3.81440307243700780478E-15,
     1.04202769841288027642E-14,
     4.27244001671195135429E-14,
    -2.10154184277266431302E-14,
    -4.08355111109219731823E-13,
    -7.19855177624590851209E-13,
     2.03562854414708950722E-12,
     1.41258074366137813316E-11,
     3.25260358301548823856E-11,
    -1.89749581235054123450E-11,
    -5.58974346219658380687E-10,
    -3.83538038596423702205E-9,
    -2.63146884688951950684E-8,
    -2.51223623787020892529E-7,
    -3.88256480887769039346E-6,
    -1.10588938762623716291E-4,
    -9.76109749136146840777E-3,
     7.78576235018280120474E-1
};

/** @brief Evaluate Chebyshev series
 *
 * Evaluates the series
 *      y = Sum( coef[i] * T_i(x/2), from i=0 to N-1)
 * of Chebyshev polynomials Ti at argument x/2.
 *
 * Coefficients are stored in reverse order, i.e. the zero order term is last
 * in the array.
 * Note: n is the number of coefficients, not the order.
 *
 * If coefficients are for the interval a to b, x must have been transformed to
 * x->2(2x-b-a)/(b-a) before entering the routine. This maps x from (a, b) to
 * (-1, 1), over which the Chebyshev polynomials are defined.
 *
 * If the coefficients are for the inverted interval, in which (a, b) is mapped
 * to (1/b, 1/a), the transformation required is x->2(2ab/x-b-a)/(b-a). If b is
 * infinity, this becomes x->4a/x-1.
 *
 * SPEED:
 * Taking advantage of the recurrence properties of the
 * Chebyshev polynomials, the routine requires one more
 * addition per loop than evaluating a nested polynomial of
 * the same degree.
 *
 * @param [in] x Chebyshev series is evaluated at this point
 * @param [in] array Chebyshev coefficients
 * @param [in] n number of Chebyshev coefficients, number of elements of array
 * @retval Chebychev series evaluated at x
 */
static double chbevl(double x, double array[], int n)
{
    double *p = array;
    double b0 = *p++, b1 = 0.0, b2;
    int i = n-1;

    do
    {
        b2 = b1;
        b1 = b0;
        b0 = x*b1-b2+*p++;
    }
    while(--i);

    return 0.5*(b0-b2);
}

/** @brief Modified Bessel function of order zero
 *
 * Returns modified Bessel function of order zero of the argument.
 *
 * The function is defined as i0(x) = j0(ix).
 *
 * The range is partitioned into the two intervals [0,8] and (8, infinity).
 * Chebyshev polynomial expansions are employed in each interval.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,30         6000       8.2e-17     1.9e-17
 *    IEEE      0,30        30000       5.8e-16     1.4e-16
 *
 * @param [in] x argument
 * @retval I0(x)
 */
double besselI0(double x)
{
    if(x < 0)
        x = -x;

    if(x <= 8.0)
    {
        double y = (x/2.0) - 2.0;
        return(exp(x)*chbevl(y,A0,30));
    }

    return exp(x)*chbevl(32.0/x-2.0,B0,25)/sqrt(x);
}


/** @brief Modified Bessel function of order zero, exponentially scaled
 *
 * Returns exponentially scaled modified Bessel function of order zero of the
 * argument.
 *
 * The function is defined as i0e(x) = exp(-|x|) j0( ix ).
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30        30000       5.4e-16     1.2e-16
 *
 * See besselI0.
 *
 * @param [in] x argument
 * @reval exp(-|x|)*I0(x)
 */
double besselI0e(double x)
{
    if(x < 0)
        x = -x;

    if(x <= 8.0)
    {
        double y = (x/2.0) - 2.0;
        return chbevl(y,A0,30);
    }

    return chbevl(32.0/x-2.0,B0,25)/sqrt(x);
}


/** @brief Modified Bessel function of order one
 *
 * Returns modified Bessel function of order one of the argument.
 *
 * The function is defined as i1(x) = -i j1( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and (8, infinity).
 * Chebyshev polynomial expansions are employed in each interval.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        3400       1.2e-16     2.3e-17
 *    IEEE      0, 30       30000       1.9e-15     2.1e-16
 *
 * @param [in] x argument
 * @reval I1(x)
 */
double besselI1(double x)
{ 
    double z = fabs(x);

    if(z <= 8.0)
    {
        double y = (z/2.0)-2.0;
        z = chbevl(y,A1,29)*z*exp(z);
    }
    else
    {
        z = exp(z)*chbevl(32.0/z-2.0,B1,25)/sqrt(z);
    }
    if( x < 0.0 )
        z = -z;
    return z;
}

/** @brief Modified Bessel function of order one, exponentially scaled
 *
 * Returns exponentially scaled modified Bessel function * of order one of the
 * argument.
 *
 * The function is defined as i1(x) = -i exp(-|x|) j1( ix ).
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       2.0e-15     2.0e-16
 *
 * See besselI1.
 *
 * @param [in] x argument
 * @reval exp(-|x|)*I1(x)
 */
double besselI1e(double x)
{ 
    double z = fabs(x);

    if(z <= 8.0)
    {
        double y = (z/2.0)-2.0;
        z = chbevl(y,A1,29)*z;
    }
    else
        z = chbevl(32.0/z-2.0,B1,25)/sqrt(z);

    if(x < 0.0)
        z = -z;

    return z;
}

#define ACC 40.0
#define BIGNO 1e10
#define BIGNI 1e-10

/** @brief Modified Bessel function of integer order
 *
 * Returns modified Bessel function of order n for the argument x.
 *
 * The function is defined as in(x) = jn( ix ).
 *
 * The algorithm is taken from Numerical Recipes in C.
 *
 * @param [in] n order
 * @param [in] x argument
 * @retval In(x)
 */
double besselI(int n, double x)
{
    if(n < 0)
        return NAN;
    if(n == 0)
        return besselI0(x);
    if(n == 1)
        return besselI1(x);

    if(x == 0)
        return 0;

    double tox = 2/fabs(x);
    double bip = 0, ans = 0, bi = 1;

    for(int j = 2*(n+(int)sqrt(ACC*n)); j > 0; j--)
    {
        double bim = bip+j*tox*bi;
        bip = bi;
        bi  = bim;

        if(fabs(bi) > BIGNO)
        {
            ans *= BIGNI;
            bi  *= BIGNI;
            bip *= BIGNI;
        }

        if(j == n)
            ans = bip;
    }

    ans *= besselI0(x)/bi;
    return x < 0.0 && (n & 1) ? -ans : ans;
}

/* @brief Calculate I_{nu+1/2}/I_{nu+3/2}
 *
 * Compute the ratio of the modified Bessel functions of the first kind
 * I_{nu+1/2}(x)/I_{nu+3/2}(x) using a continued fraction.
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

/** @brief Compute log I_{nu+1/2}(x)
 *
 * Compute logarithm of modified Bessel function of the first kind
 * I_{nu+1/2}(x).
 *
 * @param [in] nu order
 * @param [in] x argument
 * @retval log I_{nu+1/2}(x)
 */
double bessel_lnInu(int nu, double x)
{
    double lnInu;
    bessel_lnInuKnu(nu, x, &lnInu, NULL);
    return lnInu;
}

/** @brief Compute log K_{nu+1/2}(x)
 *
 * Compute logarithm of modified Bessel function of the second kind
 * K_{nu+1/2}(x).
 *
 * @param [in] nu order
 * @param [in] x argument
 * @retval log K_{nu+1/2}(x)
 */
double bessel_lnKnu(int nu, double x)
{
    double lnKnu;
    bessel_lnInuKnu(nu, x, NULL, &lnKnu);
    return lnKnu;
}

/** @brief Compute modified Bessel functions of first and second kind
 *
 * This function computes the logarithm of the modified Bessel functions
 * I_{nu+1/2}(x) and K_{nu+1/2}(x). The values are saved in lnInu_p and
 * lnKnu_p.
 *
 * If lnInu_p or lnKnu_p is NULL, the pointer is not accessed.
 *
 * @param [in] nu order
 * @param [in] x argument
 * @param [out] lnInu_p pointer for log I_{nu+1/2}(x)
 * @param [out] lnKnu_p pointer for log K_{nu+1/2}(x)
 */
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

    if(lnKnu_p != NULL)
        *lnKnu_p = lnKnu;

    if(lnInu_p != NULL)
    {
        double ratio = bessel_continued_fraction(nu,x);
        *lnInu_p = -logx-logadd(lnKnup, lnKnu-log(ratio));
    }
}
