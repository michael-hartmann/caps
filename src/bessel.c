/**
 * @file   bessel.c
 * @author Stephen L. Moshier, Cephes Math Library Release 2.8, June 2000
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   October, 2019
 * @brief  Computation of Bessel functions
 */

#include <stdio.h>

#include <stdlib.h>
#include <math.h>

#include "bessel.h"
#include "constants.h"
#include "misc.h"
#include "utils.h"

/**
 * @name modified Bessel functions for orders n=0,1
 * @{
 */

/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */
static double I0_coeffs[] =
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

/* Chebyshev coefficients for K0(x) + log(x/2) I0(x)
 * in the interval [0,2].  The odd order coefficients are all
 * zero; only the even order coefficients are listed.
 *
 * lim(x->0){ K0(x) + log(x/2) I0(x) } = -EUL.
 */
static double K0_coeffsA[] =
{
    1.37446543561352307156E-16,
    4.25981614279661018399E-14,
    1.03496952576338420167E-11,
    1.90451637722020886025E-9,
    2.53479107902614945675E-7,
    2.28621210311945178607E-5,
    1.26461541144692592338E-3,
    3.59799365153615016266E-2,
    3.44289899924628486886E-1,
    -5.35327393233902768720E-1
};

/* Chebyshev coefficients for exp(x) sqrt(x) K0(x)
 * in the inverted interval [2,infinity].
 *
 * lim(x->inf){ exp(x) sqrt(x) K0(x) } = sqrt(pi/2).
 */
static double K0_coeffsB[] = {
     5.30043377268626276149E-18,
    -1.64758043015242134646E-17,
     5.21039150503902756861E-17,
    -1.67823109680541210385E-16,
     5.51205597852431940784E-16,
    -1.84859337734377901440E-15,
     6.34007647740507060557E-15,
    -2.22751332699166985548E-14,
     8.03289077536357521100E-14,
    -2.98009692317273043925E-13,
     1.14034058820847496303E-12,
    -4.51459788337394416547E-12,
     1.85594911495471785253E-11,
    -7.95748924447710747776E-11,
     3.57739728140030116597E-10,
    -1.69753450938905987466E-9,
     8.57403401741422608519E-9,
    -4.66048989768794782956E-8,
     2.76681363944501510342E-7,
    -1.83175552271911948767E-6,
     1.39498137188764993662E-5,
    -1.28495495816278026384E-4,
     1.56988388573005337491E-3,
    -3.14481013119645005427E-2,
     2.44030308206595545468E0
};

/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
 */
static double I1_coeffs[] =
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

/* Chebyshev coefficients for x(K1(x) - log(x/2) I1(x))
 * in the interval [0,2].
 *
 * lim(x->0){ x(K1(x) - log(x/2) I1(x)) } = 1.
 */
static double K1_coeffsA[] =
{
    -7.02386347938628759343E-18,
    -2.42744985051936593393E-15,
    -6.66690169419932900609E-13,
    -1.41148839263352776110E-10,
    -2.21338763073472585583E-8,
    -2.43340614156596823496E-6,
    -1.73028895751305206302E-4,
    -6.97572385963986435018E-3,
    -1.22611180822657148235E-1,
    -3.53155960776544875667E-1,
     1.52530022733894777053E0
};

/* Chebyshev coefficients for exp(x) sqrt(x) K1(x)
 * in the interval [2,infinity].
 *
 * lim(x->inf){ exp(x) sqrt(x) K1(x) } = sqrt(pi/2).
 */
static double K1_coeffsB[] =
{
    -5.75674448366501715755E-18,
     1.79405087314755922667E-17,
    -5.68946255844285935196E-17,
     1.83809354436663880070E-16,
    -6.05704724837331885336E-16,
     2.03870316562433424052E-15,
    -7.01983709041831346144E-15,
     2.47715442448130437068E-14,
    -8.97670518232499435011E-14,
     3.34841966607842919884E-13,
    -1.28917396095102890680E-12,
     5.13963967348173025100E-12,
    -2.12996783842756842877E-11,
     9.21831518760500529508E-11,
    -4.19035475934189648750E-10,
     2.01504975519703286596E-9,
    -1.03457624656780970260E-8,
     5.74108412545004946722E-8,
    -3.50196060308781257119E-7,
     2.40648494783721712015E-6,
    -1.93619797416608296024E-5,
     1.95215518471351631108E-4,
    -2.85781685962277938680E-3,
     1.03923736576817238437E-1,
     2.72062619048444266945E0
};

/** @brief Evaluate Chebyshev series
 *
 * Evaluates the series
 *      y = Sum'( coef[i] * T_i(x/2), from i=0 to N-1)
 * of Chebyshev polynomials Ti at argument x/2. The prime indicates that the
 * term for i=0 has to be weighted by a factor 1/2.
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

/** @brief Modified Bessel function \f$I_0(x)\f$
 *
 * See \ref bessel_logI0.
 *
 * @param [in] x argument
 * @retval I0 \f$I_0(x)\f$
 */
double bessel_I0(double x)
{
    return exp(bessel_logI0(x));
}


/** @brief Logarithm of modified Bessel function \f$I_0(x)\f$
 *
 * For \f$x=0\f$ the value \f$\log I_0(0) = \log(1)=0\f$ is returned.
 *
 * For \f$0<x<8\f$ a series expansion is used, see \ref bessel_logInu_series.
 *
 * For \f$8\le x<800\f$ a Chebychev expansion is used.
 *
 * For \f$x\ge800\f$ the Hankel expansion
 * \f[
 * I_0(x) \approx \frac{e^x}{\sqrt{2\pi x}} \left( 1 + k + \frac{9}{2} k^2 + \frac{225}{6} k^3 \right), \quad k=\frac{1}{8x}
 * \f]
 * is used.
 *
 * The range is partitioned into the two intervals \f$[0,8]\f$ and \f$(8, \infty)\f$.
 * Chebyshev polynomial expansions are employed in each interval.
 *
 * @param [in] x argument
 * @retval logI0 \f$\log I_0(x)\f$
 */
double bessel_logI0(double x)
{
    if(x < 0)
        return NAN;
    if(x == 0)
        return 0; /* log(I0(0)) = log(1) = 0 */
    if(x < 8)
        return bessel_logInu_series(0, x); /* series expansion */
    if(x < 800)
        /* Chebychev expansion */
        return x+log(chbevl(32/x-2,I0_coeffs,25))-0.5*log(x);
    else
    {
        /* Hankel expansion
         * I_0(x) ≈ exp(x)/sqrt(2*pi*x) * ( 1 + k + 9/2*k² + 225/6*k³ ), k=1/8x
         */
         const double k = 1./8/x;
         return x-0.5*log(2*M_PI*x) + log1p( k*(1+9./2*k*(1+25./3*k)) );
    }
}


/** @brief Modified Bessel function \f$K_0(x)\f$
 *
 * See \ref bessel_logK0.
 *
 * @param [in] x argument
 * @retval K0 \f$K_0(x)\f$
 */
double bessel_K0(double x)
{
    return exp(bessel_logK0(x));
}


/** @brief Logarithm of modified Bessel function \f$K_0(x)\f$
 *
 * For small arguments \f$0<x<10^{-8}\f$, the limiting form
 * \f$K_0(x) \approx -\log(x/2)-\gamma\f$ for \f$x\to0\f$ where \f$\gamma\f$
 * denotes the Euler-Mascheroni constant is used.
 *
 * For large arguments \f$x>800\f$, the Hankel expansion
 * \f[
 * K_0(x) \approx \sqrt{\frac{\pi}{2x}} e^{-x} \left(1-k+\frac{9}{2} k^2 - \frac{225}{6} k^3\right), \quad k=\frac{1}{8x}
 * \f]
 * is used.
 *
 * For intermediate values, the range is partitioned into the two intervals
 * \f$[10^{-8},2)\f$ and \f$(2,800)\f$ and Chebyshev polynomial expansions are
 * employed in each interval.
 *
 * @param [in] x argument
 * @retval logK0 \f$\log K_0(x)\f$
 */
double bessel_logK0(double x)
{
    if(x <= 0)
        return NAN;
    if(x < 1e-8)
    {
        /* K_0(x) ≈ -log(x/2)-gamma */
        const double gamma = 0.5772156649015328606;
        return log(-log(x/2)-gamma);
    }
    if(x <= 2)
        return log(chbevl(x*x-2, K0_coeffsA, 10)-log(0.5*x)*bessel_I0(x));
    if(x < 800)
        return -x-0.5*log(x)+log(chbevl(8/x-2, K0_coeffsB, 25));
    else /* x >= 800 */
    {
        /* Hankel expansion
         * K_0(x) ≈ sqrt(pi/2x)*exp(-x) * ( 1 - k + 9/2*k² - 225/6*k³ ), k=1/8x
         */
         const double k = 1./8/x;
         return 0.5*log(M_PI/2/x) - x + log1p( k*(-1+9./2*k*(1-25./3*k)) );
    }

}


/** @brief Modified Bessel function \f$I_1(x)\f$
 *
 * See \ref bessel_logI1.
 *
 * @param [in] x argument
 * @retval I1 \f$I_1(x)\f$
 */
double bessel_I1(double x)
{
    return exp(bessel_logI1(x));
}


/** @brief Logarithm of modified Bessel function \f$I_1(x)\f$
 *
 * For \f$0<x<8\f$ a series expansion is used, see \ref bessel_logInu_series.
 *
 * For \f$8\le x<800\f$ a Chebychev expansion is used.
 *
 * For \f$x\ge800\f$ the Hankel expansion
 * \f[
 * I_0(x) \approx \frac{e^x}{\sqrt{2\pi x}} \left( 1 - 3k - \frac{15}{2} k^2 - \frac{105}{2} k^3 \right), \quad k=\frac{1}{8x}
 * \f]
 * is used.
 *
 * @param [in] x argument
 * @retval logI1 \f$\log I_1(x)\f$
 */
double bessel_logI1(double x)
{
    if(x <= 0)
        return NAN;
    if(x < 8)
        return bessel_logInu_series(1, x); /* series expansion */
    if(x < 800)
        /* Chebychev expansion */
        return x+log(chbevl(32/x-2,I1_coeffs,25))-0.5*log(x);
    else
    {
        /* Hankel expansion
         * I_1(x) ≈ exp(x)/sqrt(2*pi*x) * ( 1 - 3k - 15/2*k² - 105/2*k³ ), k=1/8x
         */
         const double k = 1./8/x;
         return x-0.5*log(2*M_PI*x) + log1p( -k*(3+k*(15./2+105./2*k)) );
    }
}

/** @brief Modified Bessel function \f$K_1(x)\f$
 *
 * See \ref bessel_logK1.
 *
 * @param [in] x argument
 * @retval K1 \f$K_1(x)\f$
 */
double bessel_K1(double x)
{
    return exp(bessel_logK1(x));
}

/** @brief Logarithm of modified Bessel function \f$K_1(x)\f$
 *
 * For small arguments \f$x<10^{-8}\f$, the limiting form
 * \f$K_1(x)\approx 1/x\f$ for \f$x\to0\f$ is used.
 *
 * For large arguments \f$x\le800\f$, the Hankel expansion
 * \f[
 * K_1(x) \approx \sqrt{\frac{\pi}{2x}} e^{-x} \left(1 + 3k - \frac{15}{2} k^2 + \frac{315}{6} k^3\right), \quad k = \frac{1}{8x}
 * \f]
 * is used.
 *
 * For intermediate values, the range is partitioned into the two intervals
 * \f$[10^{-8},8)\f$ and \f$[8,800)\f$ and Chebyshev polynomial expansions are
 * employed in each interval.
 *
 * @param [in] x argument
 * @retval logK1 \f$\log K_1(x)\f$
 */
double bessel_logK1(double x)
{
    if(x <= 0)
        return NAN;
    if(x < 1e-8)
        /* K_1(x) ≈ 1/x */
        return -log(x);
    if(x < 2)
        /* Chebychev expansion */
        return log(log(0.5*x)*bessel_I1(x)+chbevl(x*x-2, K1_coeffsA, 11)/x);
    if(x < 800)
        /* Chebychev expansion */
        return -x+log(chbevl(8/x-2, K1_coeffsB, 25))-0.5*log(x);
    else /* x >= 800 */
    {
        /* Hankel expansion
         * K_1(x) ≈ sqrt(pi/2x)*exp(-x) * ( 1 + 3k - 15/2 k² + 315/6 k³ )
         */
        double k = 1/(8*x);
        return 0.5*log(M_PI/2/x) -x + log1p( k*(3+k*(-15./2+315./6*k)));
    }
}

/*@}*/

/**
 * @name modified Bessel functions for integer orders
 * @{
 */

/** @brief Modified Bessel function \f$I_n(x)\f$ for integer order \f$n\f$
 *
 * See \ref bessel_logIn.
 *
 * @param [in] n order
 * @param [in] x argument
 * @retval In \f$I_n(x)\f$
 */
double bessel_In(int n, double x)
{
    return exp(bessel_logIn(n,x));
}

/** @brief Modified Bessel function \f$K_n(x)\f$ for integer order \f$n\f$
 *
 * See \ref bessel_logKn.
 *
 * @param [in] n order
 * @param [in] x argument
 * @retval Kn \f$K_n(x)\f$
 */
double bessel_Kn(int n, double x)
{
    return exp(bessel_logKn(n,x));
}

/** @brief Logarithm of modified Bessel functions \f$K_n(x)\f$ for integer orders \f$n=0,1,\dots,n_\mathrm{max}\f$
 *
 * The Bessel functions \f$K_n(x)\f$ are computed using the recurrence relation
 * \f[
 * K_{n+1}(x) = K_{n-1}(x) + \frac{2n}{x} K_n(x)
 * \f]
 * in upwards direction. The Bessel functions \f$K_0(x)\f$ and \f$K_1(x)\f$
 * are computed using \ref bessel_logK0 and \ref bessel_logK1.
 *
 * @param [in]  nmax \f$n_\mathrm{max}\f$ maximum order
 * @param [in]  x argument
 * @param [out] logKn array of n+1 elements with the values of \f$K_0(x), K_1(x),\dots, K_{n_\mathrm{max}}(x)\f$
 */
void bessel_logKn_recursive(int nmax, double x, double *logKn)
{
    if(nmax < 0)
        return;

    /* K_0(x) */
    logKn[0] = bessel_logK0(x);

    if(nmax < 1)
        return;

    /* K_1(x) */
    logKn[1] = bessel_logK1(x);

    for(int n = 1; n < nmax; n++)
    {
        /* K_{n+1} = K_{n-1} + 2n/x K_n = 2n/x * K_n * (1+x/2n*K_{n-1}/K_n) */
        double k = 0.5*x/n;
        logKn[n+1] = -log(k)+logKn[n]+log1p(exp(logKn[n-1]-logKn[n])*k);
    }
}

/** @brief Logarithm of modified Bessel function \f$K_n(x)\f$ for integer order \f$n\f$
 *
 * For \f$n=0\f$ and \f$n=1\f$ the function calls \ref bessel_logK0 or \ref
 * bessel_logK1.
 *
 * For \f$n\ge100\f$ an asymptotic expansion is used, see \ref
 * bessel_logKnu_asymp.
 *
 * Otherwise, the function is computed using a recursion relation, see \ref
 * bessel_logKn_recursive.
 *
 * @param [in]  n order
 * @param [in]  x argument
 * @retval Kn \f$\log K_n(x)\f$
 */
double bessel_logKn(int n, double x)
{
    double logKn[101];

    if(n < 0)
        return NAN;
    if(n == 0)
        return bessel_logK0(x);
    if(n == 1)
        return bessel_logK1(x);

    /* for n>=100 use asymptotic expansion */
    if(n >= 100)
        return bessel_logKnu_asymp(n, x);

    bessel_logKn_recursive(n, x, logKn);
    return logKn[n];
}

/** @brief Logarithm of modified Bessel function \f$I_n(x)\f$ for integer order \f$n\f$
 *
 * For \f$n=0\f$ and \f$n=1\f$ the function calls \ref bessel_logI0 or \ref
 * bessel_logI1.
 *
 * For \f$n>=100\f$ an asymptotic expansion is used, see \ref
 * bessel_logInu_asymp.
 *
 * For \f$n<100\f$ and \f$x<5\sqrt{n}\f$ a series expansion is used, see \ref
 * bessel_logInu_series.
 *
 * Otherwise, the functions \f$I_n(x)\f$ is computed using the recurrence relation
 * \f[
 * I_{n-1}(x) = I_{n+1}(x) + \frac{2n}{x} I_n(x)
 * \f]
 * in downwards direction using Miller's algorithm. See also \ref bessel_logI0
 * and \ref bessel_ratioI.
 *
 * @param [in]  n order
 * @param [in]  x argument
 * @retval In \f$\log I_n(x)\f$
 */
double bessel_logIn(int n, double x)
{
    if(n < 0)
        return NAN;
    if(n == 0)
        return bessel_logI0(x);
    if(n == 1)
        return bessel_logI1(x);

    if(x == 0)
        return -INFINITY;

    /* for n>=100 use asymptotic expansion */
    if(n >= 100)
        return bessel_logInu_asymp(n, x);

    /* for small values use series expansion */
    if(x < 5*sqrt(n))
        return bessel_logInu_series(n,x);

    /* miller algorithm */
    double Ip = 1; /* I_n */
    double I = bessel_ratioI(n-1,x); /* I_{n-1}/I_n */

    for(int k = n-1; k > 0; k--)
    {
        double Im = Ip + 2*k/x*I;
        Ip = I;
        I = Im;
    }

    return bessel_logI0(x)-log(I);
}

/*@}*/

/**
 * @name modified Bessel functions for arbitrary orders
 * @{
 */

/**
 * @brief Calculate \f$I_\nu(x)/I_{\nu+1}(x)\f$
 *
 * Compute the ratio of the modified Bessel functions of the first kind
 * \f$I_\nu(x)/I_{\nu+1}(x)\f$ using a continued fraction, see
 * https://dlmf.nist.gov/10.33.
 *
 * @param nu order
 * @param x argument
 * @retval ratio \f$I_\nu(x)/I_{\nu+1}(x)\f$
 */
double bessel_ratioI(double nu, double x)
{
    /* it's faster to calculate the inverse of x only once */
    const double invx2 = 2/x;

    const double a1 = (nu+1)*invx2;
    const double a2 = (nu+2)*invx2;

    double num   = a2+1/a1;
    double denom = a2;
    double ratio = a1*num/denom;
    double ratio_last = 0;

    for(int l = 3; 1; l++)
    {
        const double an = invx2*(nu+l);
        num   = an+1/num;
        denom = an+1/denom;
        ratio *= num/denom;

        if(ratio == ratio_last)
            return ratio;

        ratio_last = ratio;
    }
}


/** @brief Compute modified Bessel function \f$I_\nu(x)\f$ using asymptotic expansion
 *
 * For \f$n\ge100\f$ the asymptotic expansion is accurate.
 *
 * See also https://dlmf.nist.gov/10.41#ii.
 *
 * @param [in] nu order
 * @param [in] x argument
 * @retval logI \f$\log I_\nu(x)\f$
 */
double bessel_logInu_asymp(double nu, double x)
{
    const double z = x/nu;
    const double p2 = 1/(1+z*z); /* p² */

    const double U5 = 59535/262144.+p2*(-67608983/9175040.+p2*(250881631/5898240.+p2*(-108313205/1179648.+p2*(5391411025/63700992.-5391411025/191102976.*p2))));
    const double U4 = 3675/32768. + p2*(-96833/40960.+p2*(144001/16384.+p2*(-7436429/663552.+37182145/7962624.*p2)));
    const double U3 = 75/1024.+p2*(-4563/5120.+p2*(17017/9216.-85085/82944.*p2));
    const double U2 = 9/128.+p2*(-77/192.+385/1152.*p2);
    const double U1 = 1/8.-5/24.*p2;
    /* U0 = 1 */

    const double a = sqrt(1+z*z);
    const double y = 1/(nu*a);
    const double sum = log1p(y*(U1 + y*(U2 + y*(U3 + y*(U4 + U5*y)))));

    const double eta = a+log(z/(1+a));
    const double prefactor = nu*eta-log(2*M_PI*nu*a)/2;

    return prefactor+sum;
}


/** @brief Compute modified Bessel function \f$K_\nu(x)\f$ using asymptotic expansion
 *
 * For \f$n\ge100\f$ the asymptotic expansion is accurate.
 *
 * See also https://dlmf.nist.gov/10.41#ii.
 *
 * @param [in] nu order
 * @param [in] x argument
 * @retval logI \f$\log I_\nu(x)\f$
 */
double bessel_logKnu_asymp(double nu, double x)
{
    const double z = x/nu;
    const double p2 = 1/(1+z*z); /* p² */

    const double U5 = 59535/262144.+p2*(-67608983/9175040.+p2*(250881631/5898240.+p2*(-108313205/1179648.+p2*(5391411025/63700992.-5391411025/191102976.*p2))));
    const double U4 = 3675/32768. + p2*(-96833/40960.+p2*(144001/16384.+p2*(-7436429/663552.+37182145/7962624.*p2)));
    const double U3 = 75/1024.+p2*(-4563/5120.+p2*(17017/9216.-85085/82944.*p2));
    const double U2 = 9/128.+p2*(-77/192.+385/1152.*p2);
    const double U1 = 1/8.-5/24.*p2;
    /* U0 = 1 */

    const double a = sqrt(1+z*z);
    const double y = -1/(nu*a);
    const double sum = log1p(y*(U1 + y*(U2 + y*(U3 + y*(U4 + U5*y)))));

    const double eta = a+log(z/(1+a));
    const double prefactor = -nu*eta-0.5*log(a)+0.5*log(0.5*M_PI/nu);

    return prefactor+sum;
}


/** @brief Compute modified Bessel functions \f$I_\nu(x)\f$ using series expansion
 *
 * The modified Bessel function is computed using the series expansion
 * \f[
 * I_\nu(x) = \sum_{m=0}^\infty \frac{1}{m! \Gamma(1+m+\nu)}\left(\frac{x}{2}\right)^{2m+\nu} .
 * \f]
 *
 * The functions succeeds for orders up to \f$\nu\le100000\f$ when
 * \f$x\le10\sqrt{\nu}\f$. For larger values of \f$x\f$ the function might
 * return NAN.
 *
 * @param [in] nu order
 * @param [in] x argument
 * @retval Inu \f$I_\nu(x)\f$ if successful or NAN otherwise
 */
double bessel_logInu_series(double nu, double x)
{
    const double y = 0.25*x*x; /* x²/4 */
    double v = 1, sum = 0;

    for(int m = 1; m < 2048; m++)
    {
        v *= y/(m*(m+nu)); /* y**m */
        sum += v;

        if(v < 1e-16*sum)
            return 0.5*nu*log(y)-lgamma(1+nu)+log1p(sum);
    }

    return NAN;
}

/*@}*/

/**
 * @name modified Bessel functions for half-integer orders
 * @{
 */

/** @brief Compute modified Bessel functions of first and second kind
 *
 * This function computes the logarithm of the modified Bessel functions
 * \f$I_{\nu+1/2}(x)\f$ and \f$K_{\nu+1/2}(x)\f$. The results are saved in
 * logInu_p and logKnu_p.
 *
 * If logInu_p or logKnu_p is NULL, the variable is not referenced.
 *
 * @param [in] nu order
 * @param [in] x argument
 * @param [out] logInu_p pointer for \f$\log I_{\nu+1/2}(x)\f$
 * @param [out] logKnu_p pointer for \f$\log K_{\nu+1/2}(x)\f$
 */
void bessel_logInuKnu_half(int nu, const double x, double *logInu_p, double *logKnu_p)
{
    const double logx = log(x);
    const double invx = 1/x;

    double logKnu, logKnup;
    double Knu = 1, Knup = 1+invx;
    double prefactor = -x+0.5*(M_LOGPI-M_LOG2-logx);

    /* calculate Knu, Knup */
    if(nu == 0)
    {
        logKnu  = prefactor+log(Knu);
        logKnup = prefactor+log(Knup);
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

        logKnup = prefactor+log(Knup);
        logKnu  = prefactor+log(Knu);
    }

    if(logKnu_p != NULL)
        *logKnu_p = logKnu;

    if(logInu_p != NULL)
    {
        if(nu >= 100)
        {
            *logInu_p = bessel_logInu_asymp(nu+0.5, x);
            return;
        }

        double ratio = bessel_ratioI(nu+0.5,x);
        *logInu_p = -logx-logadd(logKnup, logKnu-log(ratio));
    }
}

/** @brief Compute \f$\log I_{\nu+1/2}(x)\f$
 *
 * Compute logarithm of modified Bessel function of the first kind
 * \f$I_{\nu+1/2}(x)\f$.
 *
 * @param [in] nu order
 * @param [in] x argument
 * @retval logI \f$\log I_{\nu+1/2}(x)\f$
 */
double bessel_logInu_half(int nu, double x)
{
    double logInu;
    bessel_logInuKnu_half(nu, x, &logInu, NULL);
    return logInu;
}

/** @brief Compute \f$\log K_{\nu+1/2}(x)\f$
 *
 * Compute logarithm of modified Bessel function of the second kind
 * \f$K_{\nu+1/2}(x)\f$.
 *
 * @param [in] nu order
 * @param [in] x argument
 * @retval logK \f$K_{\nu+1/2}(x)\f$
 */
double bessel_logKnu_half(int nu, double x)
{
    double logKnu;
    bessel_logInuKnu_half(nu, x, NULL, &logKnu);
    return logKnu;
}

/*@}*/
