#ifndef INTEGRATION_PERF_H
#define INTEGRATION_PERF_H

#include "floattypes.h"
#include "libcasimir.h"

#define TE 1
#define TM 0

/* here are all signs considered, except the sign of the evaluation of the
 * integral. But also the -1 of Lambda(l1,l2,m) is considered here.
 */
#define log_A0(m,nT) (2*log64(m)-2*(nT))
#define sign_A0(l2,m,rp) MPOW(1+(l2)+(rp))

#define log_B0(m,nT) (-2*(nT))
#define sign_B0(l2,m,rp) MPOW((l2)+(rp))

#define log_C0(m,nT) (log64(m)-2*(nT))
#define sign_C0(l2,m,rp) MPOW((l2)+(rp))

#define log_D0(m,nT) (log64(m)-2*(nT))
#define sign_D0(l2,m,rp) MPOW((l2)+(rp)+1)


typedef struct {
    float64 *cache_I;
    float64 *cache_K;
    sign_t *cache_K_signs;
    int dim_K;
    int m,lmax;
    double nT;
} integration_perf_t;

void casimir_integrate_perf_init(integration_perf_t *self, double nT, int m, int lmax);
void casimir_integrate_perf_free(integration_perf_t *self);
void casimir_integrate_perf(integration_perf_t *self, int l1, int l2, casimir_integrals_t *cint);

float64 casimir_integrate_perf_I(integration_perf_t *self, int nu);
float64 casimir_integrate_K(integration_perf_t *self, const int l1, const int l2, sign_t *sign);

#endif
