#ifndef __CASIMIR_CYLINDER
#define __CASIMIR_CYLINDER

#include "matrix.h"

typedef struct {
    double R; /* radius of the cylinder */
    double d; /* smallest separation beteween cylinder and plate */
    double H; /* a+R*/
    int lmax; /* truncation of vector space */
    int verbose;
    detalg_t detalg;
} casimir_cp_t;

typedef struct {
    double *cacheI, *cacheK1, *cacheK2;
    int lmax;
    char DN;
} kernel_args_t;

casimir_cp_t *casimir_cp_init(double R, double d);
double casimir_cp_logdetD(casimir_cp_t *self, double q, char DN);
void casimir_cp_free(casimir_cp_t *self);

int casimir_cp_get_lmax(casimir_cp_t *self);
int casimir_cp_set_lmax(casimir_cp_t *self, int lmax);

int casimir_cp_get_verbose(casimir_cp_t *self);
int casimir_cp_set_verbose(casimir_cp_t *self, int verbose);

detalg_t casimir_cp_get_detalg(casimir_cp_t *self);
int casimir_cp_set_detalg(casimir_cp_t *self, detalg_t detalg);

kernel_args_t *kernel_init(casimir_cp_t *self, double q, char DN);
void kernel_free(kernel_args_t *args);

#endif
