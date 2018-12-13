#ifndef __MISC_H
#define __MISC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <math.h>

#include "constants.h"

typedef struct
{
    sign_t s;
    double v;
} log_t;

double sqrtpm1(double x) __attribute__ ((pure));

double kahan_sum(double input[], size_t len) __attribute__ ((pure));

double logadd(const double a, const double b) __attribute__ ((pure));
double logadd_ms(log_t list[], const int len, sign_t *sign);

#ifdef __cplusplus
}
#endif

#endif
