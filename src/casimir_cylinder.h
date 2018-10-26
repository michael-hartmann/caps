#ifndef __CASIMIR_CYLINDER
#define __CASIMIR_CYLINDER

typedef struct {
    double R; /* radius of the cylinder */
    double d; /* smallest separation beteween cylinder and plate */
    double H; /* a+R*/
    int lmax; /* truncation of vector space */
} casimir_cp_t;

casimir_cp_t *casimir_cp_init(double R, double d);
double casimir_cp_dirichlet(casimir_cp_t *self, double q);
double casimir_cp_neumann(casimir_cp_t *self, double q);
void casimir_cp_free(casimir_cp_t *self);

int casimir_cp_get_lmax(casimir_cp_t *self);
int casimir_cp_set_lmax(casimir_cp_t *self, int lmax);

#endif
