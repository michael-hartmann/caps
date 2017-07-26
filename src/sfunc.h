#ifndef __SFUNC_H
#define __SFUNC_H

#include <stdlib.h>
#include <math.h>

#include "constants.h"

typedef struct
{
    sign_t s;
    double v;
} log_t;

double sqrtpm1(double x) __attribute__ ((pure));

double kahan_sum(double input[], size_t len);

double logadd(const double a, const double b);
double logadd_ms(log_t list[], const int len, sign_t *sign);

#endif
