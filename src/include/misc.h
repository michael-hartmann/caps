#ifndef __MISC_H
#define __MISC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <math.h>

#include "constants.h"

/** represent number \f$v\f$ by its sign and \f$\log|v|\f$ */
typedef struct
{
    sign_t s; /**< sign of number*/
    double v; /**< logarithm of absolute value of number */
} log_t;

double sqrtpm1(double x) __attribute__ ((pure));

double kahan_sum(double input[], size_t len) __attribute__ ((pure));

double logadd(const double a, const double b) __attribute__ ((pure));
double logadd_ms(log_t list[], const int len, sign_t *sign);

#ifdef __cplusplus
}
#endif

#endif
