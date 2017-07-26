#ifndef __SFUNC_H
#define __SFUNC_H

#include <stdlib.h>
#include <math.h>

#include "libcasimir.h"

#define pow_2(x) ((x)*(x))
#define pow_3(x) ((x)*(x)*(x))
#define pow_4(x) ((x)*(x)*(x)*(x))

/* calculate pow(-1,a) = -1**a */
#define MPOW(a) (1-2*((signed char)(a) & 1))

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) ((((a))>((b)))?((a)):((b)))
#endif

#define SGN(val) ((0 < (val)) - ((val) < 0))

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
