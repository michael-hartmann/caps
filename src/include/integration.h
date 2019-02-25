#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#include "libcaps.h"
#include "integration.h"

integration_t *caps_integrate_init(caps_t *self, double xi_, int m, double epsrel);
void caps_integrate_free(integration_t *integration);

double caps_integrate_I(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double caps_integrate_K(integration_t *self, int nu, polarization_t p, sign_t *sign);

double caps_integrate_A(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double caps_integrate_B(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double caps_integrate_C(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double caps_integrate_D(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);

integration_plasma_t *caps_integrate_plasma_init(caps_t *caps, double omegap, double epsrel);
double caps_integrate_plasma(integration_plasma_t *self, int l1, int l2, int m, double *ratio1, double *ratio2);
void caps_integrate_plasma_free(integration_plasma_t *self);

double K_estimate(int nu, int m, double alpha, double eps, double *a, double *b, double *approx);

#ifdef __cplusplus
}
#endif

#endif
