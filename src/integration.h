#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#include "libcasimir.h"
#include "integration.h"
#include "hash-table.h"


#define A_TE 0
#define A_TM 1
#define B_TE 2
#define B_TM 3
#define C_TE 4
#define C_TM 5
#define D_TE 6
#define D_TM 7

typedef struct {
    casimir_t *casimir;
    int n,m;
    double tau,epsrel;
    HashTable *hash_table;
} integration_t;

integration_t *casimir_integrate_init(casimir_t *self, int n, int m, double epsrel);
void casimir_integrate_free(integration_t *integration);
int casimir_integrate(integration_t *int_obj, int l1, int l2, double v[8]);

double casimir_integrate_I(integration_t *self, int l1, int l2, polarization_t p, double *prefactor);
double casimir_integrate_A(integration_t *self, int l1, int l2, polarization_t p, double *prefactor);
double casimir_integrate_B(integration_t *self, int l1, int l2, polarization_t p, double *prefactor);
double casimir_integrate_C(integration_t *self, int l1, int l2, polarization_t p, double *prefactor);
double casimir_integrate_D(integration_t *self, int l1, int l2, polarization_t p, double *prefactor);

#endif
