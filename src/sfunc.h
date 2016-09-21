#ifndef __SFUNC_H
#define __SFUNC_H

#include <stdlib.h>
#include <math.h>

#include "floattypes.h"
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

double kahan_sum(double input[], size_t len);

inline double logadd(const double a, const double b);
inline double logadd_s(const double a, const sign_t sign_a, const double b, const sign_t sign_b, sign_t *sign);
inline double logadd_m(const double list[], const int len);
inline double logadd_ms(const double list[], const sign_t signs[], const int len, sign_t *sign);

double bessel_lnInu(const int n, const double x);
double bessel_lnKnu(const int n, const double x);
void bessel_lnInuKnu(int nu, const double x, double *lnInu_p, double *lnKnu_p);

double ln_doublefact(int n);

void plm_lnPlm_array(int lmax, int m, double x, double lnplm[], sign_t sign[]);
double plm_lnPlm (int l, int m, double x, sign_t *sign);
double plm_Plm   (int l, int m, double x);
double plm_lndPlm(int l, int m, double x, sign_t *sign);
double plm_dPlm  (int l, int m, double x);

void plm_PlmPlm(int l1, int l2, int m, double x, plm_combination_t *res);

inline double gaunt_log_a0(int n, int nu, int m);
inline double gaunt_a0(int n,int nu,int m);
void gaunt(const int n, const int nu, const int m, double a_tilde[]);

inline int gaunt_qmax(const int n, const int nu, const int m);

#endif
