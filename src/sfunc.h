#ifndef __SFUNC_H
#define __SFUNC_H

#include <stdlib.h>
#include <math.h>

#include "edouble.h"
#include "libcasimir.h"

#define PI     3.141592653589793238462643383279502884197169L
#define LOG2   0.6931471805599453094172321214581765680755007L
#define LOGPI  1.14472988584940017414342735135305871164729L
#define LOG4   1.386294361119890618834464242916353136151001L

/* abbrevations for functions */
#define lngamma(x) (lgammae(x))
#define pow_2(x) ((x)*(x))
#define pow_3(x) ((x)*(x)*(x))
#define pow_4(x) ((x)*(x)*(x)*(x))
#define pow_5(x) ((x)*(x)*(x)*(x)*(x))
#define pow_6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define pow_7(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x))
#define lnfac(x)  (lgammae(1+(x)))

/* calculate pow(-1,a) = -1**a */
#define MPOW(a) (1-2*((a) & 1))

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) ((((a))>((b)))?((a)):((b)))
#endif

typedef struct {
    edouble lnPl1mPl2m;
    int sign_Pl1mPl2m;

    edouble lndPl1mPl2m;
    int sign_dPl1mPl2m;

    edouble lnPl1mdPl2m;
    int sign_Pl1mdPl2m;

    edouble lndPl1mdPl2m;
    int sign_dPl1mdPl2m;
} plm_combination_t;

void polymult(edouble p1[], int len_p1, edouble p2[], int len_p2, edouble p[]);

edouble inline logadd(const edouble a, const edouble b);
edouble inline logadd_s(const edouble a, const sign_t sign_a, const edouble b, const sign_t sign_b, sign_t *sign);
edouble inline logadd_m(const edouble list[], const size_t len);
edouble inline logadd_ms(const edouble list[], const sign_t signs[], const size_t len, sign_t *sign);

edouble inline lbinom(int n, int k);

edouble bessel_lnInu(const int n, const edouble x);
edouble bessel_lnKnu(const int n, const edouble x);
void bessel_lnInuKnu(int nu, const edouble x, edouble *lnInu_p, edouble *lnKnu_p);

double linspace(double start, double stop, int N, int i);
double logspace(double start, double stop, int N, int i);

edouble ln_doublefact(int n);

inline void plm_lnPlm_array(int lmax, int m, edouble x, edouble lnplm[], sign_t sign[]);
edouble plm_lnPlm (int l, int m, edouble x, sign_t *sign);
edouble plm_Plm   (int l, int m, edouble x);
edouble plm_lndPlm(int l, int m, edouble x, sign_t *sign);
edouble plm_dPlm  (int l, int m, edouble x);

void plm_PlmPlm(int l1, int l2, int m, edouble x, plm_combination_t *res);

edouble inline gaunt_log_a0(int n, int nu, int m);
edouble inline gaunt_a0(int n,int nu,int m);
void gaunt(const int n, const int nu, const int m, edouble a_tilde[]);

int inline gaunt_qmax(const int n, const int nu, const int m);

#endif
