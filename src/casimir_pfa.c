#include <stdio.h>

#include "libcasimir.h"

int main(int argc, char *argv[])
{
    casimir_t *casimir;
    double LbyR, pfa, T = 0;
    double omegap = INFINITY, gamma_ = 0;
    double userdata[2] = { 0 };

    if(argc < 2)
    {
        fprintf(stderr, "Usage: %s L/R [T ωp γ]\n", argv[0]);
        return 1;
    }

    LbyR = atof(argv[1]);
    if(argc > 2)
        T = atof(argv[2]);
    if(argc > 3)
        omegap = atof(argv[3]);
    if(argc > 4)
        gamma_ = atof(argv[4]);

    casimir = casimir_init(LbyR);

    if(!isinf(omegap))
    {
        userdata[0] = omegap;
        userdata[1] = gamma_;
        casimir_set_epsilonm1(casimir, casimir_epsilonm1_drude, userdata);
    }

    pfa = casimir_pfa(casimir, T);

    casimir_free(casimir);

    printf("# L/R, T, ωp, γ, F_PFA\n");
    printf("%.15g, %.15g, %.15g, %.15g, %.15g\n", LbyR, T, omegap, gamma_, pfa);

    return 0;
}
