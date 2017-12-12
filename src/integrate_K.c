#include <stdio.h>

#include "libcasimir.h"
#include "integration.h"

int main(int argc, char *argv[])
{
    double epsrel = 1e-8;

    if(argc < 4)
    {
        fprintf(stderr, "%s nu m alpha\n", argv[0]);
        return 1;
    }

    int nu       = atoi(argv[1]);
    int m        = atoi(argv[2]);
    double alpha = atof(argv[3]);
    polarization_t p = TE;
    sign_t sign;

    printf("m     = %d\n", m);
    printf("nu    = %d\n", nu);
    printf("alpha = %g\n", alpha);

    casimir_t *casimir = casimir_init(1);
    integration_t *integration = casimir_integrate_init(casimir, alpha/2, m, epsrel);
    double v = casimir_integrate_K(integration, nu, p, &sign);
    casimir_free(casimir);

    printf("%.16g, %d\n", v, sign);

    return 0;
}
