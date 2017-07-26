#include <stdlib.h>
#include <stdio.h>

#include "libcasimir.h"
#include "misc.h"

void usage(const char self[], FILE *stream);

void usage(const char self[], FILE *stream)
{
    fprintf(stream,
        "%s L/R nT m lmax\n\n"
        "Options:\n"
        "    L/R:   aspect ratio L/R\n"
        "    nT:    xi*(L+R)/c\n"
        "    m:     quantum number m\n"
        "    lmax:  lmax\n", self
    );
}

int main(int argc, char *argv[])
{
    if(argc < 5)
    {
        usage(argv[0], stderr);
        return 1;
    }

    const double LbyR = atof(argv[1]);
    const double nT = atof(argv[2]);
    const int m = atoi(argv[3]);
    const int lmax = atoi(argv[4]);

    const double factor = log10(exp(1));

    printf("# LbyR = %g\n", LbyR);
    printf("# nT = %g\n", nT);
    printf("# m = %d\n", m);
    printf("# lmax = %d\n", lmax);
    printf("#\n");
    printf("# row, column, log10|M_ij|, log10|Mij| (not symmetrized)\n");

    casimir_t *casimir = casimir_init(LbyR);
    casimir_set_ldim(casimir, lmax);

    casimir_M_t *args = casimir_M_init(casimir, m, nT);

    for(int i = 0; i < 2*lmax; i++)
    {
        sign_t dummy1, dummy2;
        double al1,al2,bl1,bl2;
        int l1 = (i % lmax) + 1;
        casimir_lnab_perf(casimir, nT, l1, &al1, &bl1, &dummy1, &dummy2);

        for(int j = 0; j < 2*lmax; j++)
        {
            int l2 = (j % lmax) + 1;

            double elem1 = log(fabs(casimir_kernel_M(i, j, args)));

            casimir_lnab_perf(casimir, nT, l2, &al2, &bl2, &dummy1, &dummy2);

            /* EE: sqrt( al1*al2 )
             * EM: sqrt( al1*bl2 )
             * ME: sqrt( bl1*al2 )
             * MM: sqrt( bl1*bl2 )
             */

            double elem2;
            if(i < lmax)
            {
                if(j < lmax)
                    /* EE */
                    elem2 = elem1 + al1/2 - al2/2;
                else
                    /* EM */
                    elem2 = elem1 + al1/2 - bl2/2;
            }
            else
            {
                if(j < lmax)
                    /* ME */
                    elem2 = elem1 + bl1/2 - al2/2;
                else
                    /* MM */
                    elem2 = elem1 + bl1/2 - bl2/2;
            }

            printf("%d, %d, %g, %g\n", i, j, factor*elem1, factor*elem2);
        }
    }

    casimir_M_free(args);
    casimir_free(casimir);

    return 0;
}
