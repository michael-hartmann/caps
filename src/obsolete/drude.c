#include <stdio.h>

#include "libcasimir.h"
#include "misc.h"

int main(int argc, char *argv[])
{
    double eta = 9;
    double terms[4096] = { 0 };

    if(argc < 2)
    {
        fprintf(stderr, "Usage: %s L/R [eta]\n", argv[0]);
        exit(1);
    }

    const double LbyR = atof(argv[1]);

    if(argc > 2)
        eta = atof(argv[2]);

    const int ldim = fmax(100, ceil(eta/LbyR));

    casimir_t *casimir = casimir_init(LbyR);
    casimir_set_ldim(casimir, ldim);

    for(int m = 0; true; m++)
    {
        double v;
        casimir_logdetD0(casimir, m, 0, &v, NULL, NULL);

        if(m == 0)
        {
            printf("m=0, %.15g\n", v);
            v /= 2;
        }

        terms[m] = v;

        if(fabs(v/terms[0]) < 1e-12)
            break;
    }

    const double drude = kahan_sum(terms, sizeof(terms)/sizeof(double));

    printf("# L/R, ldim, E_HT_drude/(kB*T)\n");
    printf("%.15g, %d, %.15g\n", LbyR, ldim, drude);

    return 0;
}
