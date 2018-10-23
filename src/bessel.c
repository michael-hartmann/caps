/**
 * @file   bessel.c
 * @author Stephen L. Moshier, Cephes Math Library Release 2.8, June 2000
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   October, 2019
 * @brief  Computation of Bessel functions
 */

#include <stdlib.h>
#include <math.h>

#include "misc.h"
#include "constants.h"
#include "bessel.h"

/**
 * @name modified Bessel functions for orders n=0,1
 * @{
 */

/* Chebyshev coefficients for exp(-x) I0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I0(x) } = 1.
 */
static double I0_A[] =
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
static double I0_B[] =
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
static double K0_A[] =
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
static double K0_B[] = {
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

/* Chebyshev coefficients for exp(-x) I1(x) / x
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
 */
static double I1_A[] =
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
static double I1_B[] =
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
static double K1_A[] =
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
static double K1_B[] =
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
 * The range is partitioned into the two intervals \f$[0,8]\f$ and \f$(8, \infty)\f$.
 * Chebyshev polynomial expansions are employed in each interval.
 *
 * @param [in] x argument
 * @retval I0 \f$I_0(x)\f$
 */
double bessel_I0(double x)
{
    if(x < 0)
        x = -x;

    if(x <= 8.0)
        return exp(x)*chbevl(x/2-2,I0_A,30);

    return exp(x)*chbevl(32.0/x-2.0,I0_B,25)/sqrt(x);
}


/** @brief Modified Bessel function \f$I_0(x)\f$, exponentially scaled
 *
 * The function returns \f$\exp(-|x|) I_0(x)\f$.
 *
 * See \ref bessel_I0.
 *
 * @param [in] x argument
 * @retval I0exp \f$\exp(-|x|) I_0(x)\f$
 */
double bessel_I0e(double x)
{
    if(x < 0)
        x = -x;

    if(x <= 8.0)
        return chbevl(x/2-2,I0_A,30);

    return chbevl(32.0/x-2.0,I0_B,25)/sqrt(x);
}

/** @brief Modified Bessel function \f$K_0(x)\f$
 *
 * The range is partitioned into the two intervals \f$[0,8]\f$ and
 * \f$(8,\infty)\f$. Chebyshev polynomial expansions are employed in each
 * interval.
 *
 * @param [in] x argument
 * @retval K0 \f$K_0(x)\f$
 */
double bessel_K0(double x)
{
    double y,z;

    if(x <= 0)
        return NAN;

    if(x <= 2.0)
    {
        y = x*x-2.0;
        return chbevl(y, K0_A, 10)-log(0.5*x)*bessel_I0(x);
    }

    z = 8.0/x-2.0;
    y = exp(-x)*chbevl(z, K0_B, 25)/sqrt(x);
    return y;
}


/** @brief Modified Bessel function \f$K_0(x)\f$, exponentially scaled
 *
 * The function returns \f$\exp(x) K_0(x)\f$.
 *
 * See \ref bessel_K0.
 *
 * @param [in] x argument
 * @retval K0exp \f$\exp(x) K_0(x)\f$
 */
double bessel_K0e(double x)
{
    double y;

    if(x <= 0)
        return NAN;

    if(x <= 2)
    {
        y = x*x-2.0;
        y = chbevl(y, K0_A, 10)-log(0.5*x)*bessel_I0(x);
        return(y*exp(x));
    }

    y = chbevl(8.0/x-2.0, K0_B, 25)/sqrt(x);

    return y;
}

/** @brief Modified Bessel function \f$I_1(x)\f$
 *
 * The range is partitioned into the two intervals \f$[0,8]\f$ and \f$(8,
 * \infty)\f$. Chebyshev polynomial expansions are employed in each interval.
 *
 * @param [in] x argument
 * @retval I1 \f$I_1(x)\f$
 */
double bessel_I1(double x)
{
    double z = fabs(x);

    if(z <= 8.0)
    {
        double y = (z/2.0)-2.0;
        z = chbevl(y,I1_A,29)*z*exp(z);
    }
    else
        z = exp(z)*chbevl(32.0/z-2.0,I1_B,25)/sqrt(z);
    if(x < 0)
        z = -z;
    return z;
}

/** @brief Modified Bessel function \f$I_1(x)\f$, exponentially scaled
 *
 * The function returns \f$\exp(-|x|) I_1(x)\f$.
 *
 * See \ref bessel_I1.
 *
 * @param [in] x argument
 * @retval I1e \f$\exp(-|x|) I_1(x)\f$
 */
double bessel_I1e(double x)
{
    double z = fabs(x);

    if(z <= 8.0)
    {
        double y = (z/2.0)-2.0;
        z = chbevl(y,I1_A,29)*z;
    }
    else
        z = chbevl(32.0/z-2.0,I1_B,25)/sqrt(z);

    if(x < 0.0)
        z = -z;

    return z;
}

/** @brief Modified Bessel function \f$K_1(x)\f$
 *
 * The range is partitioned into the two intervals \f$[0,2]\f$ and \f$(2,
 * \infty)\f$. Chebyshev polynomial expansions are employed in each interval.
 *
 * @param [in] x argument
 * @retval K1 \f$K_1(x)\f$
 */
double bessel_K1(double x)
{
    double y, z;

    z = 0.5*x;
    if(z <= 0)
        return NAN;

    if(x <= 2.0)
    {
        y = x*x-2.0;
        return log(z)*bessel_I1(x)+chbevl(y, K1_A, 11)/x;
    }

    return exp(-x)*chbevl(8.0/x-2.0, K1_B, 25)/sqrt(x);
}

/** @brief Modified Bessel function \f$K_1(x)\f$, exponentially scaled
 *
 * The function returns \f$\exp(x) K_1(x)\f$.
 *
 * See \ref bessel_K1.
 *
 * @param [in] x argument
 * @retval K1exp \f$\exp(x) K_1(x)\f$
 */
double bessel_K1e(double x)
{
    double y;

    if(x <= 0)
        return NAN;

    if(x <= 2)
    {
        y = x*x-2.0;
        y = log(0.5*x)*bessel_I1(x)+chbevl(y, K1_A, 11)/x;
        return y*exp(x);
    }

    return chbevl(8.0/x-2.0, K1_B, 25)/sqrt(x);
}

/*@}*/

/**
 * @name modified Bessel functions for integer orders
 * @{
 */

#define ACC 40.0
#define BIGNO 1e10
#define BIGNI 1e-10

/** @brief Modified Bessel function \f$I_n(x)\f$ for integer order \f$n\f$
 *
 * For \f$n<0\f$ NAN is returned.
 *
 * The algorithm is taken from Numerical Recipes in C.
 *
 * @param [in] n order
 * @param [in] x argument
 * @retval In \f$I_n(x)\f$
 */
double bessel_In(int n, double x)
{
    if(n < 0)
        return NAN;
    if(n == 0)
        return bessel_I0(x);
    if(n == 1)
        return bessel_I1(x);

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

    ans *= bessel_I0(x)/bi;
    return x < 0 && (n & 1) ? -ans : ans;
}

/** @brief Logarithm of modified Bessel function \f$K_n(x)\f$ for integer orders
 *
 * The Bessel functions are computed using a recurrence relation.
 *
 * @param [in]  n order
 * @param [in]  x argument
 * @param [out] K_n(x) array of n+1 elements with the values of \f$K_0(x), K_1(x),\dots, K_n(x)\f$
 */
void log_besselKn_array(int n, double x, double out[])
{
    /* K_0(x) */
    {
        if(x < 1e-100)
        {
            /* K_0(x) ≈ -log(x/2)-gamma */
            double gamma = 0.5772156649015328606;
            out[0] = log(-log(x/2)-gamma);
        }
        else if(x > 650)
        {
            /* Hankel expansion
             * K_0(x) ≈ sqrt(pi/2x)*exp(-x) * ( 1 - k + 9/2*k² - 225/6*k³ ), k=1/8x
             */
             double k = 1./8/x;
             out[0] = 0.5*log(M_PI/2/x) - x + log1p( k*(-1+9./2*k*(1-25./3*k)) );
        }
        else
            out[0] = log(bessel_K0(x));
    }

    if(n == 0)
        return;

    /* K_1(x) */
    {
        if(x < 1e-8)
            /* K_1(x) ≈ 1/x */
            out[1] = -log(x);
        else if(x > 600)
        {
            /* Hankel expansion
             * K_1(x) ≈ sqrt(pi/2x)*exp(-x) * ( 1 + 3/8x - 15/128x² + 315/3072x³ )
             */
             double k = 1/(8*x);
             out[1] = 0.5*log(M_PI/2/x) -x + log1p(3*k*(1+5./2*k *(-1+21./3*k)));
        }
        else
            out[1] = log(bessel_K1(x));
    }

    for(int l = 1; l < n; l++)
    {
        /* K_{n+1} = K_{n-1} + 2n/x K_n = 2n/x * K_n * (1+x/2n*K_{n-1}/K_n) */
        double k = 0.5*x/l;
        out[l+1] = -log(k)+out[l]+log1p(exp(out[l-1]-out[l])*k);
    }
}

/** @brief Logarithm of modified Bessel function \f$K_n(x)\f$ for integer order \f$n\f$
 *
 * The Bessel functions are computed using a recurrence relation.
 *
 * @param [in]  n order
 * @param [in]  x argument
 * @retval Kn \f$\log K_n(x)\f$
 */
double log_besselKn(int n, double x)
{
    double *out = malloc((n+1)*sizeof(double));
    log_besselKn_array(n, x, out);
    double v = out[n];
    free(out);
    return v;
}

/*@}*/

/**
 * @name modified Bessel functions for non-integer orders
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
double bessel_continued_fraction(double nu, double x)
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


/** @brief Compute \f$\log I_{\nu+1/2}(x)\f$
 *
 * Compute logarithm of modified Bessel function of the first kind
 * \f$I_{\nu+1/2}(x)\f$.
 *
 * @param [in] nu order
 * @param [in] x argument
 * @retval logI \f$\log I_{\nu+1/2}(x)\f$
 */
double bessel_lnInu(int nu, double x)
{
    double lnInu;
    bessel_lnInuKnu(nu, x, &lnInu, NULL);
    return lnInu;
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
double bessel_lnKnu(int nu, double x)
{
    double lnKnu;
    bessel_lnInuKnu(nu, x, NULL, &lnKnu);
    return lnKnu;
}

/** @brief Compute modified Bessel function \f$I_\nu(x)\f$ using asymptotic expansion
 *
 * See https://dlmf.nist.gov/10.41#ii
 *
 * @param [in] order
 * @param [in] argument
 * @param [out] relerror estimated relative error
 * @retval logI \f$\log I_\nu(x)\f$
 */
static double __lnbesselI_asymp(double nu, double x, double *relerror)
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

	*relerror = U5*y*y*y*y*y;

    return prefactor+sum;
}

/** @brief Compute modified Bessel functions of first and second kind
 *
 * This function computes the logarithm of the modified Bessel functions
 * \f$I_{\nu+1/2}(x)\f$ and \f$K_{\nu+1/2}(x)\f$. The results are saved in
 * lnInu_p and lnKnu_p.
 *
 * If lnInu_p or lnKnu_p is NULL, the variable is not referenced.
 *
 * @param [in] nu order
 * @param [in] x argument
 * @param [out] lnInu_p pointer for \f$\log I_{\nu+1/2}(x)\f$
 * @param [out] lnKnu_p pointer for \f$\log K_{\nu+1/2}(x)\f$
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
		if(nu > 100)
        {
            double relerr;
            *lnInu_p = __lnbesselI_asymp(nu+0.5, x, &relerr);

            if(relerr < 1e-12)
                return;
        }

        double ratio = bessel_continued_fraction(nu+0.5,x);
        *lnInu_p = -logx-logadd(lnKnup, lnKnu-log(ratio));
    }
}

/*@}*/
