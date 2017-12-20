/* exponentially convergent Fourier-Chebshev quadrature scheme for a
 * semi-infinite interval [0,inf]
 */

#include <math.h>
#include <stdio.h>

#include "constants.h"

/* @brief Cotangent
 *
 * Compute cotangent of x, cos(x)/sin(x)
 *
 * @param [in] x argument
 * @retval cot(x)
 */
static double cot(double x)
{
    return cos(x)/sin(x);
}


/* @brief Weights for quadrature scheme
 *
 * The weights correspond to (3.2e) of [1]. Here we have used that
 * cos(j*pi)=(-1)**j.
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


/* @brief Integrate function over interval [0,inf]
 *
 * This method uses an exponentially convergent Fourier-Chebshev quadrature to
 * compute the integral over the interval [0,inf]. The method approximately
 * doubles the number of nodes until the desired accuracy is achieved.
 *
 * Values of ier after integration:
 *   ier=0: evaluation successful
 *   ier=1: relative accuracy epsrel must be positive
 *   ier=2: integrand returned NAN
 *   ier=3: integrand returned +inf or -inf
 *   ier=4: could not acchieve desired accuracy
 *
 * @param [in]     f integrand
 * @param [in]     args NULL pointer given to f when called
 * @param [in,out] epsrel on begin desired accuracy, afterwards achieved accuracy
 * @param [in]     neval number of evaluations of integrand (may be set to NULL)
 * @param [in]     L boosting parameter
 * @param [out]    ier exit code
 * @retval I numerical value of integral
 */ 
double fcgs_semiinf(double f(double, void *), void *args, double *epsrel, int *neval, double L, int *ier)
{
    /* MMIN and MMAX must be chosen in a way that there exists a positive
     * integer k such that MMAX = MMIN * 2**k.
     */
    #define MMIN 5
    #define MMAX 2560

    /* initialize cache */
    double f_cache[MMAX];
    for(size_t i = 0; i < sizeof(f_cache)/sizeof(f_cache[0]); i++)
        f_cache[i] = NAN;

    if(neval)
        *neval = 0;

    if(epsrel < 0)
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
                const double xi = L*pow_2(cot(ti/2));

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
