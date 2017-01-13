#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#include <stdbool.h>

#include "libcasimir.h"
#include "integration.h"
#include "hash-table.h"


typedef struct {
    casimir_t *casimir;
    int m;
    double tau,epsrel;
    HashTable *hash_table_I;
    HashTable *hash_table_K;
    bool is_pc;
} integration_t;

integration_t *casimir_integrate_init(casimir_t *self, double xi, int m, double epsrel);
void casimir_integrate_free(integration_t *integration);

double casimir_integrate_I(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double casimir_integrate_K(integration_t *self, int nu, polarization_t p, sign_t *sign);

double casimir_integrate_A(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double casimir_integrate_B(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double casimir_integrate_C(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);
double casimir_integrate_D(integration_t *self, int l1, int l2, polarization_t p, sign_t *sign);

#endif
