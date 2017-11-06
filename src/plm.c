/**
 * @file   plm.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   July, 2017
 * @brief  computation of Legendre and associated Legendre polynomials
 */

#include <math.h>
#include <stdbool.h>

#include "constants.h"
#include "plm.h"
#include "logfac.h"
#include "misc.h"
#include "bessel.h"
#include "utils.h"

/**
 * @brief Associated Legendre polynomials for argument x > 1
 *
 * This function calculates associated Legendre functions for m >= 0 and x > 0.
 *
 * Associated Legendre polynomials are defined as follows:
 * \f[
 *  P_l^m(x) = (-1)^m (1-x^2)^{m/2} \frac{d}{dx^m} P_l(x)
 * \f]
 * where \f$P_l(x)\f$ denotes a Legendre polynomial.
 *
 * As \f$P_l(x)\f$ are ordinary polynomials, the only problem is the term
 * \f$(1-x^2)^{m/2}\f$ when extending the domain to values of x > 1. We will use the
 * convention \f$\sqrt{-x} = +i \sqrt{x}\f$.
 *
 * Note: Products of associated legendre polynomials with common m are
 * unambiguous, because \f$(+i)^2 = (-i)^2 = -1\f$.
 *
 * What this functions actually computes:
 * \f[
 *     P_l^m(x) = (x^2-1)^{m/2} \frac{d}{dx^m} P_l(x) .
 * \f]
 * Note that we don't include the factor \f$i^m\f$ in our calculuation.
 *
 * For (l-m) <= 200 we use an upwards recurrence relation, see \ref lnPlm_upwards, otherwise we use a
 * downwards recurrence relation, see Plm_downwards .
 *
 * See also https://en.wikipedia.org/wiki/Associated_Legendre_polynomials .
 *
 * @param [in] l degree
 * @param [in] m order
 * @param [in] x argument
 * @retval log(Plm(x))
 */
double lnPlm(int l, int m, double x)
{
    if(m == 0)
        return lnPl(l, x);
    else if((l-m) <= 200)
        return lnPlm_upwards(l, m, x);
    else
        return lnPlm_downwards(l, m, x);
}

/**
 * @brief Associated Legendre polynomials using upwards recurrence relation
 *
 * The values of Plm are calculated from \f$P_m^m(x)\f$ to \f$P_l^m(x)\f$.
 * The associated Legendre polynomials are calculated using the recurrence
 * relation http://dlmf.nist.gov/14.10.E3 with 
 * \f[ 
 *      P_m^m(x) = \frac{(2m)!}{2^m m!} (x^2-1)^{m/2}
 * \f]
 * (http://dlmf.nist.gov/14.7.E15).
 *
 * @param [in] l degree
 * @param [in] m order
 * @param [in] x argument
 * @retval logPlm \f$\log P_l^m(x)\f$
 */
double lnPlm_upwards(int l, int m, double x)
{
    double array[l-m+1];
    /* P_m^m = (2m)!/(2^m*m!) (x²-1)^(m/2), http://dlmf.nist.gov/14.7.E15 */
    double log_prefactor = lfac(2*m)-m*log(2)-lfac(m) + m/2.*log((x+1)*(x-1));

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

/**
 * @brief Estimate value of \f$P_l^m(x)\f$ for \f$x \gg 1\f$
 *
 * This function computes the value of \f$P_l^m(x)\f$ using an approximation
 * for large arguments \f$x \gg 1\f$ and returns the logarithm of the estimate.
 *
 * @param [in] l l
 * @param [in] m m
 * @param [in] x argument
 * @retval estimate \f$ \approx \log(P_l^m(x))\f$
 */
double lnPlm_estimate(int l, int m, double x)
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


/**
 * @brief Compute Legendre polynomial \f$P_l(x)\f$
 *
 * Evaluation of \f$P_l(x)\f$ for x>=1.
 *
 * For l < 100 a recurrence relation is used (see _Pl3), otherwise asymptotic
 * expansions are used (see _Pl1 and _Pl2).
 *
 * The function returns \f$\log P_l(x)\f$.
 *
 * @param [in] l degree
 * @param [in] x argument
 * @retval logPl \f$\log P_l(x)\f$
 */
double lnPl(int l, double x)
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

/**
 * @brief Calculate fraction \f$P_l^{m-1}(x)/P_l^m(x)\f$
 *
 * The fraction is computed using a continued fraction, see http://dlmf.nist.gov/14.14.E1 .
 *
 * To evaluate the continued fraction, we use http://dlmf.nist.gov/1.12#E5 and
 * http://dlmf.nist.gov/1.12#E6 .
 *
 * See also Numerical Recipes in C, ch. 5.2, Evaluation of Continued Fractions
 *
 * @param [in] l degree
 * @param [in] m order
 * @param [in] x argument
 * @retval ratio \f$P_l^{m-1}(x)/P_l^m(x)\f$
 */
double Plm_continued_fraction(const long l, const long m, const double x)
{
    const double alpha = (1-1/(x*x))/4;

    /* initial values for recurrence */
    double Am = 0, Bm = 1;
    double A = x*(l-m+1)*(l+m)*alpha; /* a0 */
    double B = m;                     /* b0 */

    int j = 1+m;
    for(int i = 1; i < 2048; i++)
    {
        /* do 16 iterations */
        int k = j+16;
        for(; j < k; j++)
        {
            const double aj = (l+1-j)*(l+j)*alpha;

            double A_ = j*A + aj*Am;
            Am = A;
            A  = A_;

            double B_ = j*B + aj*Bm;
            Bm = B;
            B  = B_;
        }

        /* check if current and last approximation of continued fraction are
         * identical:       A/B ?= Am/Bm     <=>     A*Bm ?= Am*B
         * if they are identical, we're done
         */
        if(A*Bm == Am*B)
            return (B*x*sqrt(1-1/(x*x)))/(2*A);

        /* rescale */
        const double invB = 1/B;
        Am *= invB;
        A  *= invB;
        Bm *= invB;
        B   = 1;
    }

    TERMINATE(true, "l=%ld, m=%ld, x=%.15g", l,m,x);

    return NAN;
}


/**
 * @brief Compute associated Legendre polynomials using downwards recurrence relation
 *
 * First, the fraction \f$P_l^m(x)/P_l^{m-1}(x)\f$ is computed using \ref Plm_continued_fraction.
 * Then the downwards recurrence relation http://dlmf.nist.gov/14.10.E6 is used
 * from \f$P_l^m(x)\f$ to \f$P_l^0(x)\f$. Together with \f$P_l(x)\f$ (see \ref Pl) one has the solution.
 *
 * This routine is efficient if \f$l \gg m\f$.
 *
 * @param [in] l degree
 * @param [in] m order
 * @param [in] x argument
 * @retval logPlm \f$\log P_l^m(x)\f$
 */
double lnPlm_downwards(int l, int m, double x)
{
    double vp,vm,c;
    double prefactor = lnPl(l,x); /* value of Pl(x) */

    if(m == 0)
        return prefactor;

    /* herbie XXX */
    c = 2*x/sqrt((x-1)*(x+1));

    vp = +1;
    vm = -Plm_continued_fraction(l,m,x);

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


/** @brief Compute 1st and 2nd logarithmic derivative of associated Legendre polynomial
 *
 * Compute d/dx ln(Plm(x)) and d²/dx² ln(Plm(x)).
 *
 * If d2lnPlm is NULL, the 2nd logarithmic derivative will not be computed.
 *
 * @param [in] l degree
 * @param [in] m order
 * @param [in] x position
 * @param [out] d2lnPlm 2nd logarithmic derivative of Plm(x)
 * @retval dlnPlm first logarithmic derivative of Plm(x)
 */
double dlnPlm(int l, int m, double x, double *d2lnPlm)
{
    double p = 0, q = 0, df;
    const double c2 = (x+1)*(x-1);
    const double c  = sqrt(c2);

    if(m+1 <= l)
        p = 1/Plm_continued_fraction(l,m+1,x);

    df = p/c + m*x/c2;

    if(d2lnPlm != NULL)
    {
        if(m+2 <= l)
            q = p/Plm_continued_fraction(l,m+2,x);
        *d2lnPlm = (q - p*p - m*(x*x+1)/c2)/c2;
    }

    return df;
}
