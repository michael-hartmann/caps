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

typedef struct {
    /* a = k/N, b = (k+1)/N */
    int k,N;
    double K15[8];
    double err[8];
} interval_t;

void casimir_integrate_integrands(integration_t *int_obj, double t, int l1, int l2, double log_prefactor, double v[8]);

int casimir_integrate_init(casimir_t *self, integration_t *int_obj, int n, int m);
int casimir_integrate_free(integration_t *int_obj);
int casimir_integrate(integration_t *int_obj, int l1, int l2, double v[8]);

#endif
