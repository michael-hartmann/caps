#ifndef INTEGRATION_PERF_H
#define INTEGRATION_PERF_H

#include "edouble.h"

typedef struct {
    int nu_max;
    int m2_max;
    edouble tau;
    edouble *cache_I;
} integration_perf_t;

void casimir_integrate_perf_init(integration_perf_t *self, double nT, int lmax);
void casimir_integrate_perf_free(integration_perf_t *self);
void casimir_integrate_perf(integration_perf_t *self, int l1, int l2, int m, casimir_integrals_t *cint);

#endif
