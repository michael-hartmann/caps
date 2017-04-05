#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#include <stdbool.h>

#include "libcasimir.h"
#include "integration.h"

typedef struct {
    double LbyR, omegap, epsrel;
    HashTable *cache;
} integration_plasma_t;

integration_t *casimir_integrate_init(casimir_t *self, double xi, int m, double epsrel);
void casimir_integrate_free(integration_t *integration);

double casimir_integrate_I(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double casimir_integrate_K(integration_t *self, int nu, polarization_t p, sign_t *sign);

double casimir_integrate_A(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double casimir_integrate_B(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double casimir_integrate_C(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double casimir_integrate_D(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);

integration_plasma_t *casimir_integrate_plasma_init(double RbyL, double omegap, double epsrel);
double casimir_integrate_plasma(integration_plasma_t *self, int l1, int l2, int m);
void casimir_integrate_plasma_free(integration_plasma_t *self);

#endif
