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
    float64 lnPl1mPl2m;
    int sign_Pl1mPl2m;

    float64 lndPl1mPl2m;
    int sign_dPl1mPl2m;

    float64 lnPl1mdPl2m;
    int sign_Pl1mdPl2m;

    float64 lndPl1mdPl2m;
    int sign_dPl1mdPl2m;
} plm_combination_t;

float64 kahan_sum(float64 input[], size_t len);

inline float64 logadd(const float64 a, const float64 b);
inline float64 logadd_s(const float64 a, const sign_t sign_a, const float64 b, const sign_t sign_b, sign_t *sign);
inline float64 logadd_m(const float64 list[], const int len);
inline float64 logadd_ms(const float64 list[], const sign_t signs[], const int len, sign_t *sign);

float64 bessel_lnInu(const int n, const float64 x);
float64 bessel_lnKnu(const int n, const float64 x);
void bessel_lnInuKnu(int nu, const float64 x, float64 *lnInu_p, float64 *lnKnu_p);

float64 ln_doublefact(int n);

void plm_lnPlm_array(int lmax, int m, float64 x, float64 lnplm[], sign_t sign[]);
float64 plm_lnPlm (int l, int m, float64 x, sign_t *sign);
float64 plm_Plm   (int l, int m, float64 x);
float64 plm_lndPlm(int l, int m, float64 x, sign_t *sign);
float64 plm_dPlm  (int l, int m, float64 x);

void plm_PlmPlm(int l1, int l2, int m, float64 x, plm_combination_t *res);

inline float64 gaunt_log_a0(int n, int nu, int m);
inline float64 gaunt_a0(int n,int nu,int m);
void gaunt(const int n, const int nu, const int m, float64 a_tilde[]);

inline int gaunt_qmax(const int n, const int nu, const int m);

#endif
