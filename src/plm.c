/**
 * @file   plm.c
 * @author Michael Hartmann <caps@speicherleck.de>
 * @date   January, 2019
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
 * @brief Associated Legendre polynomials for argument \f$x>1\f$
 *
 * This function calculates associated Legendre functions for \f$m \ge 0\f$ and
 * \f$x>1\f$.
 *
 * The associated Legendre polynomials for \f$x>1\f$ are defined as follows
 * (see references)
 * \f[
 *     P_l^m(x) = (x^2-1)^{m/2} \frac{\mathrm{d}}{\mathrm{d}x^m} P_l(x) .
 * \f]
 * Note that in contrast to the common choice in physics, we omit the
 * Condon-Shortly phase \f$(-1)^m\f$, and interchange the factors \f$x^2\f$ and
 * \f$1\f$ in the first bracket after the equal sign. With this definition the
 * associated Legendre polynomials are real and positive functions.
 *
 * For \f$l-m \le 200\f$ we use an upwards recurrence relation in \f$m\f$, see
 * \ref lnPlm_upwards, otherwise we use a downwards recurrence relation in
 * \f$m\f$, see \ref lnPlm_downwards .
 *
 * References:
 * - DLMF, §14.7.11, http://dlmf.nist.gov/14.7#E11
 * - Zhang, Jin, Computation of Special Functions, 1996
 *
 * @param [in] l degree
 * @param [in] m order
 * @param [in] x argument
 * @retval logPlm \f$\log P_l^m(x)\f$
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
 * The values of \f$P_l^m(x)\f$ is computed using the recurrence relation
 * http://dlmf.nist.gov/14.10.E3 in upwards direction starting from
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
    /* P_m^m = (2m)!/(2^m*m!) (x²-1)^(m/2), http://dlmf.nist.gov/14.7.E15 */
    double log_prefactor = lfac(2*m)-m*log(2)-lfac(m) + m/2.*log((x+1)*(x-1));

    if(l == m)
        return log_prefactor;

    double *array = xmalloc((l-m+1)*sizeof(double));
    array[0] = 1;
    array[1] = x*(2*m+1)*array[0];

    if(array[1] == 0)
    {
        xfree(array);
        return -INFINITY;
    }

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

    const double v = array[l-m];
    xfree(array);

    if(isnan(v))
        return NAN;

    return log_prefactor+log(fabs(v));
}


/** @brief Compute Legendre polynomial \f$\log P_l(x)\f$ for large \f$x\f$
 *
 * Evaluation of \f$\log P_l(x)\f$ for \f$x\ge1\f$ using an asymptotic expansion
 * provided that
 * \f[
 * (l+1) \sqrt{(x+1)(x-1)} \ge 25.
 * \f]
 *
 * \f$\mathcal{O}(1)\f$ computation of Legendre polynomials and Gauss-Legendre
 * nodes and weights for parallel computing, section 3.2.
 *
 * See \ref lnPl.
 *
 * @param [in] l degree
 * @param [in] x argument
 * @param [in] sinhxi \f$\sinh\xi = \sqrt{(x+1)(x-1)}\f$
 * @retval logPl \f$\log P_l(x)\f$
 */
static double _Pl1(int l, double x, double sinhxi)
{
    const int M = 17;
    const double xi = acosh(x);
    /* exp(xi) = exp(acosh(x)) = x+ sqrt((x+1)*(x-1)) */
    const double expxi = x+sinhxi;
    double Clm;

    /* Clm */
    {
        const double k = 1./l;
        const double k2 = k*k;

        /* asymptotic expansion of Gamma(l+1)/Gamma(l+3/2) with machine precision for l>=100 */
        /* Clm = (1 - 3./8*k + 25./128*k2 - 105./1024*k*k2 + 1659./32768*k4 - 6237./262144*k*k4 + 50765./4194304*k2*k4)/sqrt(l); */
        Clm = (1 + k*(-3./8 + 25./128*k + k2*(-105./1024 + 1659./32768*k + k2*(-6237./262144 + 50765./4194304*k))))/sqrt(l);
    }

    double sum = 0;
    double sinhxi_m = 1; /* sinh(xi)**m */
    double expxi_m  = 1; /* exp(xi)**m */
    for(int m = 0; m < M-1; m++)
    {
        double v = Clm*(expxi_m+exp(-(m+2*l+1)*xi))/sinhxi_m;
        sum += v;

        sinhxi_m *= sinhxi;
        expxi_m  *= expxi;
        Clm *= pow_2(m+0.5)/((2*m+2)*(l+1.5+m));
    }
    sum += Clm*(expxi_m+exp(-(2*l+M-1+1)*xi))/sinhxi_m; /* m=M-1 */

    return -(M_LOG2+M_LOGPI)/2 - log(sinhxi)/2 + log(sum) + (l+0.5)*xi;
}

/** see equations (3.27)-(3.31) */
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


/** @brief Compute Legendre polynomial \f$\log P_l(x)\f$ for small \f$x\f$
 *
 * Evaluation of \f$\log P_l(x)\f$ for \f$\ge1\f$ using an asymptotic expansion
 * provided that
 * \f[
 *      (l+1)\sqrt{(x+1)(x-1)} < 25.
 * \f]
 *
 * \f$\mathcal{O}(1)\f$ computation of Legendre polynomials and Gauss-Legendre
 * nodes and weights for parallel computing, section 3.3.
 *
 * See \ref lnPl.
 *
 * @param [in] l
 * @param [in] x
 * @retval logPl \f$\log P_l(x)\f$
 */
static double _Pl2(int l, double x)
{
    const double v = l+0.5;
    const double xi = acosh(x);

    const double y = xi*v;
    double hn[13];

    double yn = 1; /* y^n */
    for(int n = 0; n < 13; n++)
    {
        /* yn = (-y)**n */
        hn[n] = yn*bessel_In(n,y);
        yn *= -y;
    }

    double s = 0;
    double vn = 1; /* v**n */
    for(int n = 0; n < 13; n += 2)
    {
        s += _fn(n,hn)/vn;
        vn *= v*v;
    }

    return log(s);
}


/** @brief Compute Legendre polynomial \f$\log P_l(x)\f$ using recurrence relation
 *
 * Evaluation of \f$\log P_l(x)\f$ for \f$x\ge1\f$ using the recurrence
 * relation
 * \f[
 *      (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x).
 * \f]
 *
 * See \ref lnPl.
 *
 * @param [in] l order
 * @param [in] x argument
 * @retval logPl \f$\log P_l(x)\f$
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
 * @brief Compute Legendre polynomial \f$\log P_l(x)\f$
 *
 * Evaluation of \f$\log P_l(x)\f$ for \f$x\ge1\f$.
 *
 * For \f$l < 100\f$ a recurrence relation is used (see \ref _Pl3), otherwise
 * asymptotic expansions are used (see \ref _Pl1 and \ref _Pl2).
 *
 * The function returns \f$\log P_l(x)\f$.
 *
 * Reference:
 * - Bogaert, Michiels, Fostier, O(1) Computation of Legendre Polynomials and
 *   Gauss--Legendre Nodes and Weights for Parallel Computing, SIAM J. Sci.
 *   Comput. 3, 34 (2012)
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
 * http://dlmf.nist.gov/1.12#E6.
 *
 * See also Numerical Recipes in C, chapter 5.2, Evaluation of Continued Fractions.
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
 * First, the fraction \f$P_l^m(x)/P_l^{m-1}(x)\f$ is computed using \ref
 * Plm_continued_fraction.  Then the downwards recurrence relation
 * http://dlmf.nist.gov/14.10.E6 is used from \f$P_l^m(x)\f$ to \f$P_l^0(x)\f$.
 * Together with \f$P_l(x)\f$ (see \ref lnPl) one can compute \f$P_l^m(x)\f$.
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
    double prefactor = lnPl(l,x); /* value of Pl(x) */

    if(m == 0)
        return prefactor;

    /* 2x/sqrt(x²-1) */
    const double c = 2*x/sqrt((x-1)*(x+1));

    double vp = +1;
    double vm = Plm_continued_fraction(l,m,x);

    for(int mm = m-2; mm >= 0; mm--)
    {
        /* prevent integer overflows in denominator */
        double v = (vp+(mm+1)*vm*c)/((l+1.+mm)*(l-mm));
        vp = vm;
        vm = v;

        if(vm < 1e-100)
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
 * Compute \f$\frac{\mathrm{d}}{\mathrm{d}x} \log P_l^m(x)\f$ and
 * \f$\frac{\mathrm{d}^2}{\mathrm{d}x^2} \log P_l^m(x)\f$.
 *
 * If d2lnPlm is NULL, the 2nd logarithmic derivative will not be computed.
 *
 * @param [in] l degree
 * @param [in] m order
 * @param [in] x argument
 * @param [out] d2lnPlm 2nd logarithmic derivative of \f$P_l^m(x)\f$
 * @retval dlnPlm first logarithmic derivative of \f$P_l^m(x)\f$
 */
double dlnPlm(int l, int m, double x, double *d2lnPlm)
{
    double p = 0, q = 0, df;
    const double c2 = (x+1)*(x-1); /* x²-1 */
    const double c  = sqrt(c2); /* sqrt(x²-1) */

    if(m+1 <= l)
        p = 1/Plm_continued_fraction(l,m+1,x); /* p = P_l^{m+1}/P_l^m */

    df = p/c + m*x/c2; /* d/dx log(Plm(x)) */

    if(d2lnPlm != NULL)
    {
        if(m+2 <= l)
            q = p/Plm_continued_fraction(l,m+2,x); /* P_l^{m+1}/P_l^m */
        *d2lnPlm = (q - p*p - m*(x*x+1)/c2)/c2; /* d²/dx² log(Plm(x)) */
    }

    return df;
}
