#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#include "libcasimir.h"

#define A0(l1,l2,m) (MPOW((l2)+(m))*pow_2(m))
#define B0(l1,l2,m) (MPOW((l2)+(m)+1))
#define C0(l1,l2,m) (MPOW((l2)+(m)))

typedef struct {
    casimir_t *casimir;

    int m;
    double nT;
    int lmax;
} integration_t;

typedef struct {
    sign_t signs[8];
    double a,b;
    double K15[8];
    double err[8];
    double maxerr;
} interval_t;

int casimir_integrate_init(casimir_t *self, integration_t *int_obj, double nT, int m);
int casimir_integrate_free(integration_t* int_obj);
int casimir_integrate(integration_t *int_obj, int l1, int l2, casimir_integrals_t *cint);

#endif
