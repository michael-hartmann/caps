/**
 * @file   fcqs.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   December, 2018
 * @brief  exponentially convergent Fourier-Chebshev quadrature scheme (experimental)
 */

#include <math.h>
#include <stdio.h>

#include "constants.h"
#include "fcqs.h"

/** MMIN and MMAX must be chosen in a way that there exists a positive
 * integer k such that MMAX = MMIN * 2**k.
 */
#define MMIN 5
#define MMAX 2560

static double cot2(double x) __attribute__ ((pure));
static double wi_semiinf(double ti, double L, double N) __attribute__ ((pure));
static double wi_finite(double ti, double N) __attribute__ ((pure));

/** @brief Squared cotangent
 *
 * Compute square of cotangent of \f$x\f$, i.e. \f$(\cos x/\sin x)^2\f$
 *
 * @param [in] x argument
 * @retval cot2 \f$\mathrm{cot}^2(x)\f$
 */
static double cot2(double x)
{
    const double cot = cos(x)/sin(x);
    return cot*cot;
}


/** @brief Weights for quadrature scheme (semiinfinite interval)
 *
 * The weights correspond to (3.2e) of [1]. Here we have used that
 * \f$\cos(j\pi)=(-1)^j\f$.
 *
 * References:
 * - [1] Boyd, Exponentially Convergent Fourier-Chebychev Quadrature Schemes on
 * Bounded and Infinite Intervals, Journal of Scientific Computing, Vol. 2, No.
 * 2 (1987)
 *
 * @param [in] ti node
 * @param [in] L  boosting parameter
 * @param [in] N  order / number of points
 * @retval wi weight
 */
static double wi_semiinf(double ti, double L, double N)
{
    double sum = 0;
    for(int j = 1; j <= N; j+=2)
        sum += sin(j*ti)/j;

    return 8*L*sin(ti)/pow_2(1-cos(ti))/(N+1)*sum;
}

/** @brief Weights for quadrature scheme (infinite interval)
 *
 * The weights correspond to (3.1e) of [1]. Here we have used that
 * \f$\cos(j\pi)=(-1)^j\f$.
 *
 * References:
 * - [1] Boyd, Exponentially Convergent Fourier-Chebychev Quadrature Schemes on
 * Bounded and Infinite Intervals, Journal of Scientific Computing, Vol. 2, No.
 * 2 (1987)
 *
 * @param [in] ti node
 * @param [in] N  order / number of points
 * @retval wi weight
 */
static double wi_finite(double ti, double N)
{
    double sum = 0;
    for(int j = 1; j <= N; j+=2)
        sum += sin(j*ti)/j;

    return 4*sin(ti)/(N+1)*sum;
}

/** @brief Integrate function \f$f(x)\f$ over interval \f$[0,\infty)\f$
 *
 * This method uses an adaptive exponentially convergent Fourier-Chebshev
 * quadrature to compute the integral over the interval \f$[0,\infty)\f$. The
 * method approximately doubles the number of nodes until the desired accuracy
 * is achieved.
 *
 * Values of ier after integration:
 * * ier=0: evaluation successful
 * * ier=1: relative accuracy epsrel must be positive
 * * ier=2: integrand returned NAN
 * * ier=3: integrand returned +inf or -inf
 * * ier=4: could not achieve desired accuracy
 *
 * @param [in]     f integrand
 * @param [in]     args pointer given to f when called
 * @param [in,out] epsrel on begin desired accuracy, afterwards achieved accuracy
 * @param [in]     neval number of evaluations of integrand (may be set to NULL)
 * @param [in]     L boosting parameter
 * @param [out]    ier exit code
 * @retval integral numerical value of integral
 */
double fcqs_semiinf(double f(double, void *), void *args, double *epsrel, int *neval, double L, int *ier)
{
    /* initialize cache */
    double f_cache[MMAX];
    for(size_t i = 0; i < sizeof(f_cache)/sizeof(f_cache[0]); i++)
        f_cache[i] = NAN;

    if(neval)
        *neval = 0;

    if(*epsrel <= 0)
    {
        *ier = 1;
        return NAN;
    }

    int M = MMIN;
    double Ilast = NAN;
    while(M <= MMAX)
    {
        const int ratio = MMAX/M;
        const int N = M-1;
        double I = 0;

        for(int i = 1; i <= N; i++)
        {
            const int index = ratio*i;
            double fti = f_cache[index];

            const double ti = M_PI*i/(N+1);

            /* if necessary compute fti=f(ti) */
            if(isnan(fti))
            {
                const double xi = L*cot2(ti/2);

                if(neval)
                    *neval = *neval+1;
                fti = f(xi, args);

                f_cache[index] = fti;
                if(isnan(fti))
                {
                    *ier = 2;
                    return Ilast;
                }
                if(isinf(fti))
                {
                    *ier = 3;
                    return Ilast;
                }
            }

            I += wi_semiinf(ti,L,N)*fti;
        }

        /* check estimated accuracy */
        if(!isnan(Ilast))
        {
            const double eps = fabs(1-Ilast/I);
            if(eps < *epsrel)
            {
                *ier = 0;
                *epsrel = eps;
                return I;
            }
        }

        Ilast = I;
        M *= 2;
    }

    /* could not achieve accuracy */
    *ier = 4;

    return Ilast;
}

/** @brief Integrate function \f$f(x)\f$ over interval \f$[a,b]\f$
 *
 * This method uses an adaptive exponentially convergent Fourier-Chebshev
 * quadrature to compute the integral over the interval \f$[a,b]\f$. The method
 * approximately doubles the number of nodes until the desired accuracy is
 * achieved.
 *
 * Values of ier after integration:
 * * ier=0: evaluation successful
 * * ier=1: relative accuracy epsrel must be positive
 * * ier=2: integrand returned NAN
 * * ier=3: integrand returned +inf or -inf
 * * ier=4: could not achieve desired accuracy
 *
 * @param [in]     f integrand
 * @param [in]     args pointer given to f when called
 * @param [in]     a left border of integration
 * @param [in]     b right border of integration
 * @param [in,out] epsrel on begin desired accuracy, afterwards achieved accuracy
 * @param [out]    neval number of evaluations of integrand (may be set to NULL)
 * @param [out]    ier exit code
 * @retval integral numerical value of integral
 */
double fcqs_finite(double f(double, void *), void *args, double a, double b, double *epsrel, int *neval, int *ier)
{
    const double dx = (b-a)/2;

    /* initialize cache */
    double f_cache[MMAX];
    for(size_t i = 0; i < sizeof(f_cache)/sizeof(f_cache[0]); i++)
        f_cache[i] = NAN;

    if(neval)
        *neval = 0;

    if(*epsrel <= 0)
    {
        *ier = 1;
        return NAN;
    }

    int M = MMIN;
    double Ilast = NAN;
    while(M <= MMAX)
    {
        const int ratio = MMAX/M;
        const int N = M-1;
        double I = 0;

        for(int i = 1; i <= N; i++)
        {
            const int index = ratio*i;
            double fxi = f_cache[index];

            const double ti = M_PI*i/(N+1);

            /* if necessary compute fti=f(ti) */
            if(isnan(fxi))
            {
                const double xi = (cos(ti)*(b-a)+a+b)/2;

                if(neval)
                    *neval = *neval+1;
                fxi = f(xi, args);

                f_cache[index] = fxi;
                if(isnan(fxi))
                {
                    *ier = 2;
                    return Ilast;
                }
                if(isinf(fxi))
                {
                    *ier = 3;
                    return Ilast;
                }
            }

            I += wi_finite(ti,N)*fxi*dx;
        }

        /* check estimated accuracy */
        if(!isnan(Ilast))
        {
            const double eps = fabs(1-Ilast/I);
            if(eps < *epsrel)
            {
                *ier = 0;
                *epsrel = eps;
                return I;
            }
        }

        Ilast = I;
        M *= 2;
    }

    /* could not achieve accuracy */
    *ier = 4;

    return Ilast;
}
