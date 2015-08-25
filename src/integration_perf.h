#ifndef INTEGRATION_PERF_H
#define INTEGRATION_PERF_H

#include "edouble.h"
#include "libcasimir.h"

typedef struct {
    int nu_max;
    int m2_max;
    edouble tau;
    edouble *cache_I;
    
    edouble **cache;
    sign_t **signs;
    int N,elems,m,lmax;
} integration_perf_t;

void casimir_integrate_perf_init(integration_perf_t *self, double nT, int lmax);
void casimir_integrate_perf_free(integration_perf_t *self);
void casimir_integrate_perf(integration_perf_t *self, int l1, int l2, int m, casimir_integrals_t *cint);

edouble *cache_gaunt_get(integration_perf_t *cache, int n, int nu, int m, sign_t **signs);

edouble I(integration_perf_t *self, int nu, int m2);

#endif
