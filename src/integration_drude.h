#ifndef __INTEGRATION_DRUDE_H
#define __INTEGRATION_DRUDE_H

#include "libcasimir.h"
#include "edouble.h"

typedef struct {
    sign_t sign_A;
    edouble lnA_TE;
    edouble lnA_TM;

    sign_t sign_B;
    edouble lnB_TE;
    edouble lnB_TM;
    
    sign_t sign_C;
    edouble lnC_TE;
    edouble lnC_TM;
    
    sign_t sign_D;
    edouble lnD_TE;
    edouble lnD_TM;
} integrands_drude_t;

void casimir_integrate_drude(casimir_t *self, casimir_integrals_t *cint, int l1, int l2, int m, int n, double T);

void integrands_drude(edouble x, integrands_drude_t *integrands, casimir_t *self, double nT, int l1, int l2, int m);

double log_polyintegrate(edouble p[], size_t len, int l1, int l2, int m, double tau, sign_t *sign);

#endif
