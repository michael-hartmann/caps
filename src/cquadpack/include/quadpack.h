#ifndef QUADPACK_H
#define QUADPACK_H

#define GK_7_15   1
#define GK_10_21  2
#define GK_15_31  3
#define GK_20_41  4
#define GK_25_51  5
#define GK_30_61  6


/* DQAGI - Integration over (semi-) infinite intervals. (From QUADPACK)
 *
 *    Adaptive integration routine which handles functions
 *    to be integrated between -infinity to +infinity, or
 *    between either of those limits and some finite,
 *    real boundary.
 *
 *    The adaptive strategy compares results of integration
 *    over the interval with the sum of results obtained from
 *    integration of bisected interval. Since error estimates
 *    are available from each regional integration, the interval
 *    with the largest error is bisected and new results are
 *    computed. This bisection process is continued until the
 *    error is less than the prescribed limit or convergence
 *    failure is determined.
 *
 *    Note that bisection, in the sense used above, refers to
 *    bisection of the transformed interval.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    bound - optional finite bound on integral.
 *
 *    inf - specifies range of integration as follows:
 *        inf = -1 -- range is from -infinity to bound,
 *        inf =  1 -- range is from bound to +infinity,
 *        inf =  2 -- range is from -infinity to +infinity,
 *                (bound is immaterial in this case).
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 */
double dqagi(double f(double, void *), double bound, int inf, double epsabs, double epsrel, double *abserr, int *neval, int *ier, void *user_data);


/* DQAGS - Integration over finite intervals. (From QUADPACK)
 *
 *    Adaptive integration routine which handles functions
 *    to be integrated between two finite bounds.
 *
 *    The adaptive strategy compares results of integration
 *    over the given interval with the sum of results obtained
 *    from integration over a bisected interval. Since error
 *    estimates are available from each regional integration, the
 *    region with the largest error is bisected and new results
 *    are computed. This bisection process is continued until the
 *    error is less than the prescribed limit or convergence
 *    failure is determined.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 */
double dqags(double f(double, void *), double a, double b, double epsabs, double epsrel, double *abserr, int *neval, int *ier, void *user_data);


/* DQAGE - Approximation to definite integral. (From QUADPACK)
 *
 *    Allows user's choice of Gauss-Kronrod integration rule.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 *
 *    irule - integration rule to be used as follows:
 *        irule = 1 -- G_K 7-15
 *        irule = 2 -- G_K 10-21
 *        irule = 3 -- G_K 15-31
 *        irule = 4 -- G_K 20-41
 *        irule = 5 -- G_K 25-51
 *        irule = 6 -- G_K 30-61
 *
 *    limit - maximum number of subintervals.
 */
double dqage(double f(double, void *), double a, double b, double epsabs, double epsrel, int irule, double *abserr, int *neval, int *ier, int *last, void *user_data);

#endif
