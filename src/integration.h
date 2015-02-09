#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#include "libcasimir.h"
#include "edouble.h"

typedef struct {
    int sign_A;
    edouble lnA_TE;
    edouble lnA_TM;

    int sign_B;
    edouble lnB_TE;
    edouble lnB_TM;
    
    int sign_C;
    edouble lnC_TE;
    edouble lnC_TM;
    
    int sign_D;
    edouble lnD_TE;
    edouble lnD_TM;
} integrands_drude_t;

void casimir_integrate_perf(casimir_integrals_t *cint, int l1, int l2, int m, int n, double T);
void casimir_integrate_drude(casimir_t *self, casimir_integrals_t *cint, int l1, int l2, int m, int n, double T);

void integrands_drude(edouble x, integrands_drude_t *integrands, casimir_t *self, double nT, int l1, int l2, int m);

double log_polyintegrate(edouble p[], size_t len, int l1, int l2, int m, double tau, int *sign);
void polym(edouble p[], int m);
void polyplm(edouble pl1[], edouble pl2[], int l1, int l2, int m);
void polydplm(edouble pl1[], edouble pl2[], int l1, int l2, int m);

void casimir_integrate_coefficients(int l1, int l2, int m, edouble pmppl1mppl2m[], edouble pmpdpl1mpdpl2m[], edouble pmpdpl1mppl2m[], edouble pmppl1mpdpl2m[]);

#endif
