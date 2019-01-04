/**
 * @file   quadpack.h
 * @date   January, 2019
 * @brief  library for numerical integration of one-dimensional functions
*/

#ifndef QUADPACK_H
#define QUADPACK_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define GK_7_15   1 /**< Gauss-Kronrod 7-15 rule */
#define GK_10_21  2 /**< Gauss-Kronrod 10-21 rule */
#define GK_15_31  3 /**< Gauss-Kronrod 15-31 rule */
#define GK_20_41  4 /**< Gauss-Kronrod 20-41 rule */
#define GK_25_51  5 /**< Gauss-Kronrod 25-51 rule */
#define GK_30_61  6 /**< Gauss-Kronrod 30-61 rule */


/** @brief Integration over (semi-) infinite intervals
 *
 * Adaptive integration routine which handles functions to be integrated
 * between -infinity to +infinity, or between either of those limits and some
 * finite, real boundary.
 *
 * The adaptive strategy compares results of integration over the interval with
 * the sum of results obtained from integration of bisected interval. Since
 * error estimates are available from each regional integration, the interval
 * with the largest error is bisected and new results are computed. This
 * bisection process is continued until the error is less than the prescribed
 * limit or convergence failure is determined.
 *
 * Note that bisection, in the sense used above, refers to bisection of the
 * transformed interval.
 *
 * error messages:
 * - ier=0: Normal and reliable termination of the routine. It is assumed that
 *   the requested accuracy has been achieved.
 * - ier=1: Maximum number of subdivisions allowed has been achieved. It is
 *   advised to analyze the integrand in order to determine the integration
 *   difficulties.
 * - ier=2: The occurrence of roundoff error is detected, which prevents the
 *   requested tolerance from being achieved. The error may be under-estimated.
 * - ier=3: Extremely bad integrand behaviour occurs at some points of the
 *   integration interval.
 * - ier=4: The algorithm does not converge. Roundoff error is detected in the
 *   extrapolation table. It is assumed that the requested tolerance cannot be
 *   achieved, and that the returned result is the best which can be obtained.
 * - ier=5: The integral is probably divergent, or slowly convergent. It must
 *   be noted that divergence can occur with any other value of ier.
 * - ier=6: The input is invalid.
 *
 * @param [in]  f       double precision function to be integrated
 * @param [in]  bound   optional finite bound on integral
 * @param [in]  inf     specifies range of integration as follows:
 *          * inf=-1: range is from -infinity to bound,
 *          * inf=+1: range is from bound to +infinity,
 *          * inf=+2: range is from -infinity to +infinity, (bound is ignored in this case)
 * @param [in]  epsabs  absolute accuracy requested
 * @param [in]  epsrel  relative accuracy requested
 * @param [out] abserr  estimate of the modulus of the absolute error, which should equal or exceed abs(I-result)
 * @param [out] neval   number of integrand evaluations
 * @param [out] ier     error message; ier=0 for normal and reliable termination, otherwise ier>0
 * @param [out] user_data   pointer that will be passed as second argument to integrand function f
 * @retval result       approximation to the integral
 */
double dqagi(double f(double, void *), double bound, int inf, double epsabs, double epsrel, double *abserr, int *neval, int *ier, void *user_data);


/** @brief Integration over finite intervals
 *
 * Adaptive integration routine which handles functions to be integrated
 * between two finite bounds.
 *
 * The adaptive strategy compares results of integration over the given
 * interval with the sum of results obtained from integration over a bisected
 * interval. Since error estimates are available from each regional
 * integration, the region with the largest error is bisected and new results
 * are computed. This bisection process is continued until the error is less
 * than the prescribed limit or convergence failure is determined.
 *
 * error messages:
 * - ier=1 Maximum number of subdivisions allowed has been achieved. It is
 *   advised to analyze the integrand in order to determine the integration
 *   difficulties.
 * - ier=2: The occurrence of roundoff error is detected, which prevents the
 *   requested tolerance from being achieved. The error may be
 *   under-estimated.
 * - ier=3: Extremely bad integrand behaviour occurs at some points of the
 *   integration interval.
 * - ier=4: The algorithm does not converge. Roundoff error is detected in the
 *   extrapolation table. It is presumed that the requested tolerance cannot be
 *   achieved, and that the returned result is the best which can be obtained.
 * - ier=5: The integral is probably divergent, or slowly convergent. It must
 *   be noted that divergence can occur with any other value of ier.
 * - ier=6: The input is invalid.
 *
 * @param [in]  f       double precision function to be integrated
 * @param [in]  a       lower limit of integration
 * @param [in]  b       upper limit of integration
 * @param [in]  epsabs  absolute accuracy requested
 * @param [in]  epsrel  relative accuracy requested
 * @param [out] abserr  estimate of the modulus of the absolute error, which should equal or exceed abs(I-result)
 * @param [out] neval   number of integrand evaluations
 * @param [out] ier     error message; ier=0 for normal and reliable termination, otherwise ier>0
 * @param [out] user_data pointer that will be passed as second argument to integrand function f
 * @retval result       approximation to the integral
 */
double dqags(double f(double, void *), double a, double b, double epsabs, double epsrel, double *abserr, int *neval, int *ier, void *user_data);


/** @bierf Approximation to definite integral
 *
 * Allows user's choice of Gauss-Kronrod integration rule.
 *
 * error messages:
 * - ier=1: Maximum number of subdivisions allowed has been achieved. It is
 *   advised to analyze the integrand in order to determine the integration
 *   difficulties.
 * - ier=2: The occurrence of roundoff error is detected, which prevents the
 *   requested tolerance from being achieved.
 * - ier=3: Extremely bad integrand behaviour occurs at some points of the
 *   integration interval.
 * - ier=6: The input is invalid.

 *
 * @param [in]  f       double precision function to be integrated
 * @param [in]  a       lower limit of integration
 * @param [in]  b       upper limit of integration
 * @param [in]  epsabs  absolute accuracy requested
 * @param [in]  epsrel  relative accuracy requested
 * @param [in]  epsrel  relative accuracy requested
 * @param [in]  irule   integration rule to be used (GK_7_15, GK_7_15, GK_10_21, GK_15_31, GK_20_41, GK_25_51, or GK_30_61)
 * @param [out] abserr  estimate of the modulus of the absolute error, which should equal or exceed abs(I-result)
 * @param [out] neval   number of integrand evaluations
 * @param [out] ier     error message; ier=0 for normal and reliable termination, otherwise ier>0
 * @param [out] last    number of subintervals actually produced in the subdivision process
 * @param [out] user_data pointer that will be passed as second argument to integrand function f
 * @retval result       approximation to the integral
 */
double dqage(double f(double, void *), double a, double b, double epsabs, double epsrel, int irule, double *abserr, int *neval, int *ier, int *last, void *user_data);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
