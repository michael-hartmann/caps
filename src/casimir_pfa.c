#include <stdio.h>

#include "libcasimir.h"

int main(int argc, char *argv[])
{
    casimir_t *casimir;
    double LbyR, pfa, T = 0;

    if(argc < 2)
    {
        fprintf(stderr, "Usage: %s L/R\n", argv[0]);
        return 1;
    }

    LbyR = atof(argv[1]);
    if(argc > 2)
        T = atof(argv[2]);

    casimir = casimir_init(LbyR);
    pfa = casimir_pfa(casimir, T);
    casimir_free(casimir);

    printf("# L/R, T, F_PFA\n");
    printf("%.15g, %.15g, %.15g\n", LbyR, T, pfa);

    return 0;
}
