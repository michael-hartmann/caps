#ifndef __INTEGRATION_DRUDE_H
#define __INTEGRATION_DRUDE_H

#include "libcasimir.h"
#include "floattypes.h"

typedef struct {
    sign_t sign_A;
    float80 lnA_TE;
    float80 lnA_TM;

    sign_t sign_B;
    float80 lnB_TE;
    float80 lnB_TM;
    
    sign_t sign_C;
    float80 lnC_TE;
    float80 lnC_TM;
    
    sign_t sign_D;
    float80 lnD_TE;
    float80 lnD_TM;
} integrands_drude_t;

void casimir_integrate_drude(casimir_t *self, casimir_integrals_t *cint, int l1, int l2, int m, int n, double T);

void integrands_drude(float80 x, integrands_drude_t *integrands, casimir_t *self, double nT, int l1, int l2, int m);

double log_polyintegrate(float80 p[], size_t len, int l1, int l2, int m, double tau, sign_t *sign);

#endif
