#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "besselI.h"
#include "plm.h"
#include "sfunc.h"


/* Legendre polynomial Pl
 *
 * Evaluation of Pl(x) for x>=1 using an asymptotic expansion provided that
 *      (l+1)*sqrt((x+1)*(x-1)) >= 25.
 *
 * O(1) computation of Legendre polynomials and Gauss-Legendre nodes and
 * weights for parallel computing, section 3.2.
 */
static double _Pl1(int l, double x)
{
    const double xi = acosh(x);
    const double expxi = exp(xi);
    const double sinhxi = sqrt((x+1)*(x-1));

    double sum = 0;
    double sinhxi_m = 1; /* sinh(xi)**m */
    double expxi_m  = 1; /* exp(xi)**m */
    for(int m = 0; m < 18; m++)
    {
        //return exp(2*lgamma(m+0.5)+lgamma(l+1) - log(M_PI) - m*log(2) -lgamma(l+m+1.5) -lgamma(m+1));
        double Clm = exp( 2*lfac(2*m) - 3*lfac(m) + lfac(l) - lfac(2*(l+m+1)) - M_LOGPI/2 + (2*(l+1)-3*m)*log(2) + lfac(l+m+1)  );

        sum += Clm * (expxi_m+exp(-(m+2*l+1)*xi)) /sinhxi_m;
        sinhxi_m *= sinhxi;
        expxi_m  *= expxi;
    }

    return log(2/(M_PI*sinhxi))/2 + log(sum) + (l+0.5)*xi - log(2);
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
 *      (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - x P_{n-1}(x).
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


/* Legendre polynomial Pl
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
        return _Pl1(l,x);
    else
        /* (l+1)*sinh(xi) < 25. */
        return _Pl2(l,x);
}


/* derivative of Legendre polynomial Pl
 *
 * Evaluation of dPl(x) for x>=1.
 *
 * Function returns log(dPl(x)).
 *
 * @param [in] l degree
 * @param [in] x argument
 * @retval log(dPl(x))
 */
double dPl(int l, double x)
{
    double lnPl  = Pl(l,x);   /* P_l */
    double lnPlm = Pl(l-1,x); /* P_{l-1} */

    return log((l*x)/((x+1)*(x-1))) + lnPl + log1p( -exp(lnPlm-lnPl)/x );
}

/* Calculate fraction P_l^m/P_l^{m-1}
 *
 * The fraction is computed using a continued fraction, see http://dlmf.nist.gov/14.14.E1 .
 *
 * To evaluate the continued fraction, we use http://dlmf.nist.gov/1.12#E5 and
 * http://dlmf.nist.gov/1.12#E6 .
 */
double _cf(int l, int m, double x)
{
    const double c = (x+1)*(x-1)/4;
    double Amm = 1, Am = 0, Bmm = 0;
    double last = 0;

    for(int n = 0; n < 10000; n++)
    {
        double f, an,bn, A, B;

        an = (l-m-n+1)*(l+m+n)*c;
        bn = (m+n)*x;

        A = Am*bn + Amm*an;
        B = 1/(bn + Bmm*an);

        f = A*B;

        Amm = Am*B;
        Am = A*B;

        Bmm = B;

        if(fabs(1-last/f) < 1e-16)
            return 2*f/sqrt((x+1)*(x-1));

        last = f;
    }

    return NAN;
}


/* Evaluate associated Legendre polynomial */
double Plm2(int l, int m, double x)
{
    if(l < 200 || ((double)m/l) > 0.01)
        return Plm(l,m,x,1,1);

    if(m == 0)
        return Pl(l,x);

    double root = 1/sqrt((x+1)*(x-1));
    double prefactor = Pl(l,x);
    double v0 = 1;
    double v1 = -exp(dPl(l,x)-prefactor)/root;

    for(int mm = 1; mm < m; mm++)
    {
        double v = (l+mm)*(l-mm+1)*v0 + 2*mm*x*v1*root;
        v0 = v1;
        v1 = v;

        if(fabs(v) > 1e150)
        {
            prefactor += log(1e150);
            v0 /= 1e150;
            v1 /= 1e150;
        }
    }

    return prefactor+log(fabs(v1));
}
