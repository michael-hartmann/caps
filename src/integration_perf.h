#ifndef INTEGRATION_PERF_H
#define INTEGRATION_PERF_H

#include "edouble.h"
#include "libcasimir.h"

typedef struct {
    edouble **cache;
    sign_t **signs;
    int N,elems;
} gaunt_cache_t;

typedef struct {
    int nu_max;
    int m2_max;
    edouble tau;
    edouble *cache_I;
    gaunt_cache_t *cache_gaunt;
} integration_perf_t;

void casimir_integrate_perf_init(integration_perf_t *self, double nT, int lmax);
void casimir_integrate_perf_free(integration_perf_t *self);
void casimir_integrate_perf(integration_perf_t *self, int l1, int l2, int m, casimir_integrals_t *cint);

void cache_gaunt_init(gaunt_cache_t *cache, int lmax);
void cache_gaunt_free(gaunt_cache_t *cache);
edouble *cache_gaunt_get(gaunt_cache_t *cache, int n, int nu, int m, sign_t **signs);

#endif
