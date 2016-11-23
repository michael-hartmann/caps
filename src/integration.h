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

/* maximum number of intervals */
#define INTEGRATE_INTERVALS_MAX 200

/* XXX check signs! XXX */
/*  sign of Lambda --\
 *                   |
 *                   V */
#define A0(l1,l2,m) (-MPOW((l2)+(m)))
#define B0(l1,l2,m) (-MPOW((l2)+(m)+1))
#define C0(l1,l2,m) (-MPOW((l2)+(m)))
#define D0(l1,l2,m) (-MPOW((l2)+(m)+1))

typedef struct {
    casimir_t *casimir;
    int n,m;
    double nT,tau,log_tau;
    HashTable *hash_table;
} integration_t;

integration_t *casimir_integrate_init(casimir_t *self, int n, int m);
void casimir_integrate_free(integration_t *integration);
int casimir_integrate(integration_t *int_obj, int l1, int l2, double v[8]);

double casimir_integrate_I(integration_t *self, int l1, int l2, polarization_t p, double *prefactor);
double casimir_integrate_A(integration_t *self, int l1, int l2, polarization_t p, double *prefactor);
double casimir_integrate_B(integration_t *self, int l1, int l2, polarization_t p, double *prefactor);
double casimir_integrate_C(integration_t *self, int l1, int l2, polarization_t p, double *prefactor);
double casimir_integrate_D(integration_t *self, int l1, int l2, polarization_t p, double *prefactor);

#endif
