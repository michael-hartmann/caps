/**
 * @file   misc.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   July, 2017
 * @brief  various mathematical functions
 */

#include <stdlib.h>
#include <math.h>

#include "misc.h"

/**
 * @brief Compute sum of array elements
 *
 * This function calculates the sum of the elements of the array input. This
 * function uses the Kahan summation algorithm to reduce numerical error.
 *
 * The algorithm is taken from Wikipedia, see
 * https://en.wikipedia.org/wiki/Kahan_summation_algorithm.
 *
 * @param [in] input array
 * @param [in] N length of array
 * @return sum sum of array elements
 */
double kahan_sum(double input[], size_t N)
{
    double sum = 0;
    double c = 0; /* running compensation for lost low-order bits */

	for(size_t i = 0; i < N; i++)
    {
        double y = input[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}

/** @brief Compute \f$\sqrt{1+x}-1\f$
 *
 * If \f$x\f$ is small, \f$\sqrt{1+x}\approx1\f$ and a loss of significance occurs
 * when calculating \f$\sqrt{1+x}-1\f$.
 *
 * For this reason we compute
 * \f[
 *     \sqrt{1+x}-1 = \frac{x}{\sqrt{1+x}+1}
 * \f]
 * to avoid a loss of significance if x is small.
 *
 * @param [in] x
 * @retval sqrt(1+x)-1
 */
double sqrtpm1(double x)
{
    return x/(sqrt(1+x)+1);
}

/**
 * @brief Add two numbers given by their logarithms.
 *
 * Both numbers are assumed to be nonnegative.
 *
 * @param [in] log_a number
 * @param [in] log_b number
 * @return log_sum \f$\log{\left[\exp({\mathrm{log\_a}})+\exp{(\mathrm{log\_b})}\right]}\f$
 */
double logadd(const double log_a, const double log_b)
{
    if(isinf(log_a) && log_a < 0)
        return log_b;
    else if(isinf(log_b) && log_b < 0)
        return log_a;

    if(log_a > log_b)
        return log_a + log1p(exp(log_b-log_a));
    else
        return log_b + log1p(exp(log_a-log_b));
}

/**
 * @brief Add N numbers given by their logarithms.
 *
 * The logarithm and the sign of the N numbers are given by list. The numbers
 * of elements of list must be N, the sign of the result will be stored in
 * sign.
 *
 * @param [in] list list of numbers given by logarithm and sign
 * @param [in] N number of elements of list
 * @return logsum log(sum_i list_i)
 */
double logadd_ms(log_t list[], const int N, sign_t *sign)
{
    double max = list[0].v;

    for(int i = 1; i < N; i++)
        if(list[i].v > max)
            max = list[i].v;

    double sum = list[0].s*exp(list[0].v-max);
    for(int i = 1; i < N; i++)
        sum += list[i].s*exp(list[i].v-max);

    *sign = SGN(sum);
    return max + log(fabs(sum));
}
