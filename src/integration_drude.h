#ifndef __INTEGRATION_DRUDE_H
#define __INTEGRATION_DRUDE_H

#include "libcasimir.h"
#include "floattypes.h"


/*
 * The accuracy of the Drude-Integration
 */
#define DRUDE_INTEG_ACCURACY 1e-9

typedef struct {
    sign_t sign_A;
    float64 lnA_TE;
    float64 lnA_TM;

    sign_t sign_B;
    float64 lnB_TE;
    float64 lnB_TM;
    
    sign_t sign_C;
    float64 lnC_TE;
    float64 lnC_TM;
    
    sign_t sign_D;
    float64 lnD_TE;
    float64 lnD_TM;
} integrands_drude_t;


// Forward Declaration
struct plm_cache;

typedef struct {
    struct plm_cache*    plm_cache;
    casimir_t*           casimir;

    int m;
    double nT;
    int lmax;
} integration_drude_t;


struct integ_context {
    casimir_t*     casimir;
    integration_drude_t* int_drude;
    
    double         nT;
    int            l1, l2, m;
    float64        c0, c_max;
};


void casimir_integrate_drude_init(casimir_t* self, integration_drude_t* int_drude, double nT, int m, int lmax);

void casimir_integrate_drude_free(integration_drude_t* int_drude);


void casimir_integrate_drude(integration_drude_t* int_drude, int l1, int l2, casimir_integrals_t* cint);

//void integrands_drude(float64 x, integrands_drude_t *integrands, casimir_t *self, double nT, int l1, int l2, int m);

double log_polyintegrate(float64 p[], size_t len, int l1, int l2, int m, double tau, sign_t *sign);

#endif
