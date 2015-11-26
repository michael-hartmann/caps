#ifndef INTEGRATION_PERF_H
#define INTEGRATION_PERF_H

#include "edouble.h"
#include "libcasimir.h"

#define TE 1
#define TM 0

#define log_A0(m,nT) (2*log(m)-2*nT)
#define sign_A0(l2,m,rp) MPOW(-1, l2+m+rp)

#define log_B0(m,nT) (-2*nT)
#define sign_B0(l2,m,rp) MPOW(-1, l2+m+1+rp)

#define log_C0(m,nT) (log(m)-2*nT)
#define sign_C0(l2,m,rp) MPOW(-1, l2+m+rp)

#define log_D0(m,nT) (log(m)-2*nT)
#define sign_D0(l2,m,rp) MPOW(-1, l2+m+rp+1)


typedef struct {
    edouble tau;
    edouble *cache_I;
    edouble *cache_K
    int dim_K;
    int m,lmax;
    double nT;
} integration_perf_t;

void casimir_integrate_perf_init(integration_perf_t *self, double nT, int m, int lmax);
void casimir_integrate_perf_free(integration_perf_t *self);

edouble casimir_integrate_I(integration_perf_t *self, int nu);
edouble casimir_integrate_K(integration_perf_t *self, const int l1, const int l2, sign_t *sign);

#endif
