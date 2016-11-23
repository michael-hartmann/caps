/**
 * @file   integration.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   November, 2016
 * @brief  Perform integration for arbitrary materials
 */

#include <math.h>

#include "utils.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"
#include "hash-table.h"


integration_t *casimir_integrate_init(casimir_t *self, int n, int m)
{
    return NULL;
}

void casimir_integrate_free(integration_t *integration)
{

}

int casimir_integrate(integration_t *int_obj, int l1, int l2, double v[8])
{
    return 0;
}

double casimir_integrate_A(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    return 0;
}

double casimir_integrate_B(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    return 0;
}

double casimir_integrate_C(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    return 0;
}

double casimir_integrate_D(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    return 0;
}

double casimir_integrate_I(integration_t *self, int l1, int l2, polarization_t p, double *prefactor)
{
    /*
    int m = self->m;
    double tau = self->tau;
    */

    return 0;
}
