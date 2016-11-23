#ifndef __SFUNC_H
#define __SFUNC_H

#include <stdlib.h>
#include <math.h>

#include "libcasimir.h"

#define pow_2(x) ((x)*(x))
#define pow_3(x) ((x)*(x)*(x))
#define pow_4(x) ((x)*(x)*(x)*(x))
#define pow_5(x) ((x)*(x)*(x)*(x)*(x))
#define pow_6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define pow_7(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x))

/* calculate pow(-1,a) = -1**a */
#define MPOW(a) (1-2*((signed char)(a) & 1))

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) ((((a))>((b)))?((a)):((b)))
#endif

typedef struct {
    double lnPl1mPl2m;
    int sign_Pl1mPl2m;

    double lndPl1mPl2m;
    int sign_dPl1mPl2m;

    double lnPl1mdPl2m;
    int sign_Pl1mdPl2m;

    double lndPl1mdPl2m;
    int sign_dPl1mdPl2m;
} plm_combination_t;

double lfac(unsigned int n);
double logi(unsigned int x);

double kahan_sum(double input[], size_t len);

double max(double input[], size_t len);
double min(double input[], size_t len);

double logadd(const double a, const double b);
double logadd_s(const double a, const sign_t sign_a, const double b, const sign_t sign_b, sign_t *sign);

double bessel_lnInu(const int n, const double x);
double bessel_lnKnu(const int n, const double x);
void bessel_lnInuKnu(int nu, const double x, double *lnInu_p, double *lnKnu_p);

double ln_factorial2(unsigned int n);
double factorial2(unsigned int n);

void Plm_array(int lmax, int m, double x, double factor, double array[]);
double Plm_estimate(int l, int m, double x);


#endif
